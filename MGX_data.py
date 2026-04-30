#!/usr/bin/env python3
"""
MGX_data.py
────────────
Three subcommands for SRA metagenomics data management:

  fetch-data  Fetch SRA run metadata by taxon ID(s) and/or BioProject accession(s).
  download    Download SRA runs via prefetch + fasterq-dump.
  parse       Scan downloaded FASTQ files and update metadata CSVs.

Usage:
    python MGX_data.py fetch-data --taxon 1510822 --email you@email.com
    python MGX_data.py fetch-data --bioproject PRJNA857725 PRJNA123456
    python MGX_data.py fetch-data --taxon 1510822 --bioproject PRJNA857725 --email you@email.com
    python MGX_data.py download   --accessions runs.csv --outdir ./fastq
    python MGX_data.py parse      --reference sra_taxid1510822.csv --outdir ./fastq

Run with -h for full argument list:
    python MGX_data.py <command> -h

Requirements:
    pip install biopython pandas requests tqdm
    SRA Toolkit (prefetch + fasterq-dump) on PATH for 'download'
    Install: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
"""

from __future__ import annotations

import argparse
import io
import json
import os
import re
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

# ── defaults ───────────────────────────────────────────────────────────────
FASTERQ_THREADS  = 8
PREFETCH_MAX_GB  = 50

ESEARCH_BATCH    = 5000   # UIDs per esearch page
RUNINFO_BATCH    = 200    # UIDs per efetch runinfo call
BIOSAMPLE_BATCH  = 100    # accessions per BioSample efetch call
ESUMMARY_BATCH   = 200    # UIDs per esummary call (BioProject UID → accession)

EFETCH_URL       = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESEARCH_URL      = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
BIOSAMPLE_URL    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
SRA_RUNS_URL     = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runs"

TARGET_STRATEGIES = {"WGS", "METAGENOMIC"}

ZLAB_COLS = [
    "Zlab_sort",
    "Zlab_SRA_path",
    "Zlab_metadata_path",
]


# ══════════════════════════════════════════════════════════════════════════
# Shared helpers
# ══════════════════════════════════════════════════════════════════════════

def _get(url: str, params: dict, timeout: int = 120) -> requests.Response:
    """GET with 5-attempt exponential backoff."""
    delay = 5
    for attempt in range(5):
        try:
            resp = requests.get(url, params=params, timeout=timeout)
            resp.raise_for_status()
            return resp
        except requests.RequestException as exc:
            if attempt == 4:
                raise
            print(f"\n  [retry {attempt+1}/5] {exc} — waiting {delay}s")
            time.sleep(delay)
            delay *= 2


def fetch_biosample_attributes(biosample_accs: list[str]) -> pd.DataFrame:
    """
    Fetch BioSample XML and extract all sample attributes.
    Returns a DataFrame keyed on BioSample accession.
    """
    print(f"\n  Fetching BioSample attributes for {len(biosample_accs):,} samples ...")
    all_records: list[dict] = []
    n_batches = (len(biosample_accs) + BIOSAMPLE_BATCH - 1) // BIOSAMPLE_BATCH

    for i in tqdm(range(0, len(biosample_accs), BIOSAMPLE_BATCH),
                  desc="  BioSample", total=n_batches):
        batch = biosample_accs[i : i + BIOSAMPLE_BATCH]
        try:
            resp = _get(BIOSAMPLE_URL,
                        {"db": "biosample", "id": ",".join(batch),
                         "rettype": "xml", "retmode": "xml"})
        except Exception as exc:
            tqdm.write(f"  WARNING: batch {i} failed — {exc}")
            continue

        time.sleep(0.4)

        try:
            root = ET.fromstring(resp.text)
        except ET.ParseError as exc:
            tqdm.write(f"  WARNING: XML parse error — {exc}")
            continue

        for bs in root.findall(".//BioSample"):
            acc = bs.get("accession", "")
            if not acc:
                continue
            record: dict = {"BioSample": acc}
            for attr in bs.findall(".//Attribute"):
                name = (attr.get("attribute_name") or
                        attr.get("harmonized_name") or "").strip()
                if name:
                    record[name] = (attr.text or "").strip()
            all_records.append(record)

    if not all_records:
        print("  WARNING: No BioSample attributes retrieved.")
        return pd.DataFrame({"BioSample": biosample_accs})

    df_bs = pd.DataFrame(all_records)
    print(f"  BioSample attribute columns: {len(df_bs.columns) - 1}")
    return df_bs


def load_accessions_two_col(path: str) -> pd.DataFrame:
    """
    Load a plain-text two-column accessions file (no header).
    Column 1: SRA run accession (SRR/ERR/DRR).
    Column 2: BioProject accession.
    Accepts tab or comma as separator.
    """
    for sep in ("\t", ","):
        try:
            df = pd.read_csv(path, sep=sep, dtype=str, header=None)
            if df.shape[1] >= 2:
                df = df.iloc[:, :2].copy()
                df.columns = ["Run", "BioProject"]
                return df.dropna().reset_index(drop=True)
        except Exception:
            continue
    sys.exit(f"[ERROR] Cannot parse {path} as two-column.\n"
             "  Expected: SRR_accession<tab or comma>BioProject  (no header)")


def load_accessions_any(path: str) -> list[str]:
    """Load single-column or two-column accessions file; return list of SRRs."""
    for sep in (",", "\t"):
        df = pd.read_csv(path, sep=sep, dtype=str, header=None)
        if df.shape[1] >= 2:
            first = str(df.iloc[0, 0]).strip().upper()
            if not (first.startswith("SRR") or first.startswith("ERR")
                    or first.startswith("DRR")):
                df = df.iloc[1:].reset_index(drop=True)
            return df.iloc[:, 0].dropna().str.strip().tolist()
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def load_bioproject_list(values: list[str]) -> list[str]:
    """
    Resolve --bioproject values to a flat list of accession strings.

    If exactly one value is given and it is an existing file path,
    accessions are read from that plain-text file (one accession per line;
    blank lines and comment lines starting with '#' are ignored).
    Otherwise the values are returned as-is (direct accessions on the
    command line).
    """
    if len(values) == 1 and Path(values[0]).is_file():
        path = values[0]
        with open(path) as fh:
            accs = [ln.strip() for ln in fh
                    if ln.strip() and not ln.startswith("#")]
        print(f"  Loaded {len(accs):,} BioProject accessions from {path}")
        return accs
    return values


def is_downloaded(srr: str, bioproject_dir: Path) -> bool:
    """Return True if any FASTQ file exists in <bioproject_dir>/<srr>/."""
    srr_dir = bioproject_dir / srr
    return srr_dir.is_dir() and bool(
        list(srr_dir.glob(f"{srr}*.fastq.gz")) or
        list(srr_dir.glob(f"{srr}*.fastq"))
    )


def scan_downloaded_srrs(location: Path) -> dict[str, str]:
    """
    Walk <location>/<BioProject>/<SRR>/ and return {SRR: BioProject}
    for every run that has FASTQ files present on disk.
    """
    result: dict[str, str] = {}
    if not location.is_dir():
        return result
    for bp_dir in sorted(location.iterdir()):
        if not bp_dir.is_dir():
            continue
        for srr_dir in sorted(bp_dir.iterdir()):
            if not srr_dir.is_dir():
                continue
            srr = srr_dir.name
            if is_downloaded(srr, bp_dir):
                result[srr] = bp_dir.name
    return result


def save_csv_merge(df_new: pd.DataFrame, path: str) -> int:
    """
    Write df_new to path. If path already exists, merge by outer-joining
    columns and deduplicating on 'Run' (keep latest).

    Column order: [existing columns] + [new columns only (not in existing)]
    Returns total row count.
    """
    p = Path(path)
    if p.exists():
        df_existing = pd.read_csv(p, dtype=str, skip_blank_lines=True)
        df_new_copy = df_new.astype(str)

        # Identify old vs new columns
        old_cols = list(df_existing.columns)
        new_cols = [c for c in df_new_copy.columns if c not in old_cols]

        # Merge and deduplicate
        df_combined = pd.concat([df_existing, df_new_copy], ignore_index=True)
        df_combined = (df_combined
                       .drop_duplicates(subset=["Run"], keep="last")
                       .reset_index(drop=True)
                       .dropna(how="all"))  # Drop rows that are completely empty

        # Reorder columns: existing first, then new at the end
        col_order = old_cols + new_cols
        df_combined = df_combined[col_order]

        # Fill NaN with empty string
        df_combined = df_combined.fillna("")
    else:
        df_combined = df_new.astype(str)
    df_combined.to_csv(path, index=False)
    return len(df_combined)


def drop_blank_columns(df: pd.DataFrame) -> pd.DataFrame:
    def all_blank(col):
        return df[col].isna().all() or df[col].astype(str).str.strip().eq("").all()
    return df[[c for c in df.columns if not all_blank(c)]]


def reorder_zlab_first(df: pd.DataFrame) -> pd.DataFrame:
    """Place Zlab columns first, then Run, then everything else."""
    present_zlab = [c for c in ZLAB_COLS if c in df.columns]
    has_run      = "Run" in df.columns
    other_cols   = [c for c in df.columns if c not in present_zlab and c != "Run"]
    return df[present_zlab + (["Run"] if has_run else []) + other_cols]


def check_tools(no_prefetch: bool = False, iseq: bool = False) -> None:
    if iseq:
        tools = ["iseq"]
    elif no_prefetch:
        tools = ["fasterq-dump", "gzip"]
    else:
        tools = ["prefetch", "fasterq-dump", "gzip"]
    missing = [t for t in tools
               if subprocess.run(["which", t], capture_output=True).returncode != 0]
    if missing:
        if iseq:
            sys.exit(
                f"[ERROR] Tool not found on PATH: {', '.join(missing)}\n"
                "  Install iSeq:  conda install bioconda::iseq\n"
                "  Source:        https://github.com/BioOmics/iSeq"
            )
        sys.exit(
            f"[ERROR] Tool(s) not found on PATH: {', '.join(missing)}\n"
            "  Install SRA Toolkit: "
            "https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit"
        )


# ══════════════════════════════════════════════════════════════════════════
# fetch-taxon subcommand
# ══════════════════════════════════════════════════════════════════════════

def _search_sra_uids(taxon_id: str, api_key: str | None) -> list[str]:
    """Search SRA by taxon ID; return all matching UIDs via usehistory pagination."""
    print(f"\n[Step 1] Searching SRA for taxon {taxon_id} ...")
    params = {
        "db": "sra", "term": f"txid{taxon_id}[Organism]",
        "usehistory": "y", "retmax": 0, "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key

    data   = _get(ESEARCH_URL, params).json()["esearchresult"]
    total  = int(data["count"])
    webenv = data["webenv"]
    qkey   = data["querykey"]
    print(f"  Total SRA records: {total:,}")

    if total == 0:
        sys.exit("  No SRA records found for this taxon.")

    uids: list[str] = []
    for start in tqdm(range(0, total, ESEARCH_BATCH), desc="  Fetching UIDs"):
        p = {
            "db": "sra", "term": f"txid{taxon_id}[Organism]",
            "webenv": webenv, "query_key": qkey,
            "retstart": start, "retmax": ESEARCH_BATCH, "retmode": "json",
        }
        if api_key:
            p["api_key"] = api_key
        uids.extend(_get(ESEARCH_URL, p).json()["esearchresult"]["idlist"])
        time.sleep(0.11 if api_key else 0.4)

    print(f"  UIDs collected: {len(uids):,}")
    return uids


def _search_bioproject_by_keyword(query: str, api_key: str | None) -> list[str]:
    """
    Search the BioProject database by free-text query.
    Returns a list of BioProject accessions (PRJNA…/PRJEB…/PRJDB…).

    Steps:
      1. eSearch db=bioproject → UIDs
      2. esummary on UIDs → resolve to accession strings
    """
    sleep = 0.11 if api_key else 0.4

    # ── Step 1: eSearch BioProject ────────────────────────────────────
    print(f"\n[Step 1] Searching BioProject for: {query!r} ...")
    params = {
        "db": "bioproject", "term": query,
        "usehistory": "y", "retmax": 0, "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key

    data   = _get(ESEARCH_URL, params).json()["esearchresult"]
    total  = int(data["count"])
    webenv = data["webenv"]
    qkey   = data["querykey"]
    print(f"  Total BioProject records: {total:,}")

    if total == 0:
        print(f"  No BioProject records found for query: {query!r}")
        return []

    uids: list[str] = []
    for start in tqdm(range(0, total, ESEARCH_BATCH), desc="  Fetching UIDs"):
        p = {
            "db": "bioproject", "term": query,
            "webenv": webenv, "query_key": qkey,
            "retstart": start, "retmax": ESEARCH_BATCH, "retmode": "json",
        }
        if api_key:
            p["api_key"] = api_key
        uids.extend(_get(ESEARCH_URL, p).json()["esearchresult"]["idlist"])
        time.sleep(sleep)

    print(f"  UIDs collected: {len(uids):,}")

    # ── Step 2: esummary → accession strings ─────────────────────────
    print(f"  Resolving UIDs to BioProject accessions ...")
    accessions: list[str] = []
    n_batches = (len(uids) + ESUMMARY_BATCH - 1) // ESUMMARY_BATCH

    for i in tqdm(range(0, len(uids), ESUMMARY_BATCH),
                  desc="  esummary", total=n_batches):
        batch = uids[i : i + ESUMMARY_BATCH]
        params_s: dict = {
            "db": "bioproject", "id": ",".join(batch), "retmode": "json",
        }
        if api_key:
            params_s["api_key"] = api_key
        try:
            result = _get(ESUMMARY_URL, params_s).json().get("result", {})
            for uid in batch:
                acc = result.get(uid, {}).get("project_acc", "")
                if acc:
                    accessions.append(acc)
        except Exception as exc:
            tqdm.write(f"  WARNING: esummary batch at index {i} failed — {exc}")
        time.sleep(sleep)

    print(f"  BioProject accessions resolved: {len(accessions):,}")
    return accessions


def _fetch_runinfo(uids: list[str], checkpoint_file: str,
                   api_key: str | None, sleep: float) -> pd.DataFrame:
    """
    Fetch SRA runinfo in batches, filter to WGS/METAGENOMIC, checkpoint/resume.

    Rows are written to a temp CSV on disk batch-by-batch instead of
    accumulating in memory, so large taxons (tens of thousands of UIDs)
    do not cause out-of-memory crashes.
    """
    tmp_csv   = checkpoint_file.replace("_checkpoint.json", "_partial.csv")
    start_idx = 0

    if os.path.exists(checkpoint_file):
        with open(checkpoint_file) as fh:
            ckpt = json.load(fh)
        start_idx = ckpt["next_idx"]
        # Migrate old checkpoint format: rows were stored inside the JSON
        if "rows" in ckpt and ckpt["rows"] and not os.path.exists(tmp_csv):
            pd.DataFrame(ckpt["rows"]).to_csv(tmp_csv, index=False)
            print(f"\n  Migrated {len(ckpt['rows']):,} rows from old checkpoint format.")
        rows_so_far = pd.read_csv(tmp_csv).shape[0] if os.path.exists(tmp_csv) else 0
        print(f"\n  Resuming from UID index {start_idx:,} "
              f"({rows_so_far:,} runs already saved)")

    print(f"\n[Step 2] Fetching runinfo ({len(uids):,} UIDs, "
          f"batches of {RUNINFO_BATCH}) ...")

    pbar = tqdm(total=len(uids), initial=start_idx,
                desc="  runinfo UIDs", unit="uid")

    for i in range(start_idx, len(uids), RUNINFO_BATCH):
        batch  = uids[i : i + RUNINFO_BATCH]
        params = {
            "db": "sra", "id": ",".join(batch),
            "rettype": "runinfo", "retmode": "text",
        }
        if api_key:
            params["api_key"] = api_key

        try:
            resp     = _get(EFETCH_URL, params)
            df_batch = pd.read_csv(
                io.StringIO(resp.text),
                engine="python",
                on_bad_lines="skip",
            ).dropna(how="all")
            mask     = df_batch["LibraryStrategy"].str.upper().isin(TARGET_STRATEGIES)
            df_batch = df_batch[mask]
            if not df_batch.empty:
                header = not os.path.exists(tmp_csv)
                df_batch.to_csv(tmp_csv, mode="a", header=header, index=False)
        except Exception as exc:
            tqdm.write(f"  WARNING: batch at index {i} failed — {exc}")

        time.sleep(sleep)
        pbar.update(len(batch))

        with open(checkpoint_file, "w") as fh:
            json.dump({"next_idx": i + RUNINFO_BATCH}, fh)

    pbar.close()

    for f in (checkpoint_file, ):
        if os.path.exists(f):
            os.remove(f)

    if not os.path.exists(tmp_csv):
        sys.exit("No WGS/METAGENOMIC runs found.")

    df = (pd.read_csv(tmp_csv)
            .drop_duplicates(subset=["Run"])
            .reset_index(drop=True))
    os.remove(tmp_csv)

    if df.empty:
        sys.exit("No WGS/METAGENOMIC runs found.")

    print(f"  WGS/METAGENOMIC runs: {len(df):,}")
    return df


def _fetch_sra_runinfo(accession: str) -> pd.DataFrame:
    """Fetch SRA runinfo from the SRA Run Selector backend."""
    print(f"\n[Step 1] Fetching SRA runinfo for {accession} ...")
    for attempt in range(5):
        try:
            resp = requests.get(
                SRA_RUNS_URL,
                params={"acc": accession, "rettype": "runinfo"},
                timeout=120,
            )
            resp.raise_for_status()
            break
        except requests.RequestException as exc:
            if attempt == 4:
                sys.exit(f"  ERROR: {exc}")
            wait = 5 * (2 ** attempt)
            print(f"  [retry {attempt+1}/5] {exc} — waiting {wait}s")
            time.sleep(wait)

    content = resp.text.strip()
    if not content or content.startswith("<"):
        print(f"  WARNING: No SRA data returned for {accession} — skipping.")
        return None

    df = pd.read_csv(io.StringIO(content)).dropna(how="all").reset_index(drop=True)
    print(f"  Runs: {len(df):,}")
    return df


def _print_summary(df: pd.DataFrame) -> None:
    """Print LibraryStrategy / Platform / total-data summary for a run DataFrame."""
    print("\nLibraryStrategy:")
    for val, n in df["LibraryStrategy"].value_counts().items():
        print(f"  {val:<25} {n:>6}  {'█' * (n * 30 // max(len(df), 1))}")

    print("\nPlatform:")
    for val, n in df["Platform"].value_counts().items():
        print(f"  {val:<25} {n:>6}")

    total_gb = pd.to_numeric(df["bases"], errors="coerce").sum() / 1e9
    print(f"\nTotal data: {total_gb:.1f} Gb")


def _apply_bases_filter(df: pd.DataFrame,
                        min_bases: int | None,
                        max_bases: int | None) -> pd.DataFrame:
    """Filter runs by base count. Returns filtered DataFrame."""
    if min_bases is None and max_bases is None:
        return df
    if "bases" not in df.columns:
        print("  WARNING: 'bases' column not found — size filter skipped.")
        return df
    bases = pd.to_numeric(df["bases"], errors="coerce")
    mask = pd.Series(True, index=df.index)
    if min_bases is not None:
        mask &= bases >= min_bases
    if max_bases is not None:
        mask &= bases <= max_bases
    n_removed = (~mask).sum()
    if n_removed:
        lo = f">={min_bases:,}" if min_bases else ""
        hi = f"<={max_bases:,}" if max_bases else ""
        bounds = " & ".join(filter(None, [lo, hi]))
        print(f"  Bases filter ({bounds} bases): removed {n_removed:,} run(s), "
              f"{mask.sum():,} remaining")
    return df[mask].reset_index(drop=True)


# ══════════════════════════════════════════════════════════════════════════
# fetch-data subcommand
# ══════════════════════════════════════════════════════════════════════════

def cmd_fetch_data(args: argparse.Namespace) -> None:
    taxon_ids   = args.taxon      or []
    bioprojects = load_bioproject_list(args.bioproject or [])
    sra_accs    = load_bioproject_list(args.accession  or [])
    keywords    = args.keyword    or []

    if not taxon_ids and not bioprojects and not sra_accs and not keywords:
        sys.exit("[ERROR] Provide at least one --taxon ID, --bioproject accession, "
                 "--accession SRA accession, or --keyword query.")

    try:
        from Bio import Entrez  # noqa: F401 — only needed for api_key plumbing
    except ImportError:
        sys.exit("[ERROR] biopython not installed.\n  pip install biopython")

    api_key = args.api_key or None
    sleep   = 0.11 if api_key else 0.4

    # ── taxon IDs ─────────────────────────────────────────────────────
    for taxon_id in taxon_ids:
        out_csv = args.merge_into if args.merge_into else f"sra_taxid{taxon_id}.csv"
        ckpt    = f"sra_taxid{taxon_id}_checkpoint.json"

        print("=" * 60)
        print(f"  fetch-data (taxon)  —  taxon {taxon_id}")
        print(f"  strategies          :  {', '.join(sorted(TARGET_STRATEGIES))}")
        print("=" * 60)

        uids    = _search_sra_uids(taxon_id, api_key)
        df_runs = _fetch_runinfo(uids, ckpt, api_key, sleep)
        df_runs = _apply_bases_filter(df_runs, args.min_bases, args.max_bases)
        if df_runs.empty:
            print(f"  No runs remaining after size filter for taxon {taxon_id} — skipping.")
            continue

        biosample_accs = (
            df_runs["BioSample"].dropna().unique().tolist()
            if "BioSample" in df_runs.columns else []
        )
        if biosample_accs:
            df_bs = fetch_biosample_attributes(biosample_accs)
            df    = df_runs.merge(df_bs, on="BioSample", how="left")
        else:
            print("  WARNING: No BioSample column — skipping attribute fetch.")
            df = df_runs

        merging = Path(out_csv).exists()
        print(f"\n[Step 4] {'Merging into' if merging else 'Saving →'} {out_csv} ...")
        total = save_csv_merge(df, out_csv)

        print(f"\n{'=' * 60}")
        print(f"  New runs fetched : {len(df):,}")
        print(f"  Total in file    : {total:,}  ({'merged' if merging else 'created'})  →  {out_csv}")
        print(f"{'=' * 60}")
        _print_summary(df)

    # ── BioProject accessions ─────────────────────────────────────────
    for raw_acc in bioprojects:
        accession = raw_acc.strip().upper()
        out_csv = args.merge_into if args.merge_into else f"{accession}_sra_runs.csv"

        print("=" * 60)
        print(f"  fetch-data (project)  —  {accession}")
        print("=" * 60)

        df_runs = _fetch_sra_runinfo(accession)
        if df_runs is None:
            sys.exit(f"[ERROR] No SRA data for {accession}. "
                     "Check the accession is valid and has public data.")

        if "LibraryStrategy" in df_runs.columns:
            n_before = len(df_runs)
            mask     = df_runs["LibraryStrategy"].str.upper().isin(TARGET_STRATEGIES)
            df_runs  = df_runs[mask].reset_index(drop=True)
            print(f"  Strategy filter ({'/'.join(sorted(TARGET_STRATEGIES))}): "
                  f"{len(df_runs):,} / {n_before:,} runs")
        if df_runs.empty:
            print(f"  No WGS/METAGENOMIC runs for {accession} — skipping.")
            continue

        df_runs = _apply_bases_filter(df_runs, args.min_bases, args.max_bases)
        if df_runs.empty:
            print(f"  No runs remaining after size filter for {accession} — skipping.")
            continue

        biosample_accs = (
            df_runs["BioSample"].dropna().unique().tolist()
            if "BioSample" in df_runs.columns else []
        )
        if biosample_accs:
            df_bs = fetch_biosample_attributes(biosample_accs)
            df    = df_runs.merge(df_bs, on="BioSample", how="left")
        else:
            print("  WARNING: No BioSample column — skipping attribute fetch.")
            df = df_runs

        merging = Path(out_csv).exists()
        print(f"\n[Step 3] {'Merging into' if merging else 'Saving →'} {out_csv} ...")
        total = save_csv_merge(df, out_csv)

        print(f"\n{'=' * 60}")
        print(f"  New runs fetched : {len(df):,}")
        print(f"  Total in file    : {total:,}  ({'merged' if merging else 'created'})  →  {out_csv}")
        print(f"{'=' * 60}")
        _print_summary(df)

        sra_cols = set(df_runs.columns)
        bs_cols  = [c for c in df.columns if c not in sra_cols]
        if bs_cols:
            print(f"\nBioSample attributes added ({len(bs_cols)}): {', '.join(bs_cols)}")

    # ── SRA accessions (SRR/SRS/SRP/SRX/SAMN/…) ──────────────────────
    for raw_acc in sra_accs:
        accession = raw_acc.strip().upper()
        out_csv = args.merge_into if args.merge_into else f"{accession}_sra_runs.csv"

        print("=" * 60)
        print(f"  fetch-data (accession)  —  {accession}")
        print("=" * 60)

        df_runs = _fetch_sra_runinfo(accession)
        if df_runs is None:
            print(f"  WARNING: No SRA data returned for {accession} — skipping.")
            continue

        df_runs = _apply_bases_filter(df_runs, args.min_bases, args.max_bases)
        if df_runs.empty:
            print(f"  No runs remaining after size filter for {accession} — skipping.")
            continue

        biosample_accs = (
            df_runs["BioSample"].dropna().unique().tolist()
            if "BioSample" in df_runs.columns else []
        )
        if biosample_accs:
            df_bs = fetch_biosample_attributes(biosample_accs)
            df    = df_runs.merge(df_bs, on="BioSample", how="left")
        else:
            print("  WARNING: No BioSample column — skipping attribute fetch.")
            df = df_runs

        merging = Path(out_csv).exists()
        print(f"\n[Step 3] {'Merging into' if merging else 'Saving →'} {out_csv} ...")
        total = save_csv_merge(df, out_csv)

        print(f"\n{'=' * 60}")
        print(f"  New runs fetched : {len(df):,}")
        print(f"  Total in file    : {total:,}  ({'merged' if merging else 'created'})  →  {out_csv}")
        print(f"{'=' * 60}")
        _print_summary(df)

        sra_cols = set(df_runs.columns)
        bs_cols  = [c for c in df.columns if c not in sra_cols]
        if bs_cols:
            print(f"\nBioSample attributes added ({len(bs_cols)}): {', '.join(bs_cols)}")

    # ── keyword queries (BioProject database) ────────────────────────
    for keyword in keywords:
        safe    = re.sub(r"[^\w]+", "_", keyword).strip("_").lower()
        out_csv = args.merge_into if args.merge_into else f"sra_kw_{safe}.csv"

        print("=" * 60)
        print(f"  fetch-data (keyword)  —  {keyword!r}")
        print(f"  strategies            :  {', '.join(sorted(TARGET_STRATEGIES))}")
        print("=" * 60)

        bp_accs = _search_bioproject_by_keyword(keyword, api_key)
        if not bp_accs:
            print(f"  Skipping — no BioProject results for {keyword!r}")
            continue

        print(f"\n[Step 2] Fetching SRA runs for {len(bp_accs):,} BioProjects ...")
        all_runs: list[pd.DataFrame] = []
        for acc in tqdm(bp_accs, desc="  BioProject runs", unit="project"):
            df_bp = _fetch_sra_runinfo(acc)
            if df_bp is None:
                continue
            if "LibraryStrategy" in df_bp.columns:
                mask  = df_bp["LibraryStrategy"].str.upper().isin(TARGET_STRATEGIES)
                df_bp = df_bp[mask]
            if not df_bp.empty:
                all_runs.append(df_bp)

        if not all_runs:
            print(f"  No WGS/METAGENOMIC runs found for keyword {keyword!r}")
            continue

        df_runs = (pd.concat(all_runs, ignore_index=True)
                     .drop_duplicates(subset=["Run"])
                     .reset_index(drop=True))
        print(f"  WGS/METAGENOMIC runs collected: {len(df_runs):,}")
        df_runs = _apply_bases_filter(df_runs, args.min_bases, args.max_bases)
        if df_runs.empty:
            print(f"  No runs remaining after size filter for keyword {keyword!r} — skipping.")
            continue

        biosample_accs = (
            df_runs["BioSample"].dropna().unique().tolist()
            if "BioSample" in df_runs.columns else []
        )
        if biosample_accs:
            df_bs = fetch_biosample_attributes(biosample_accs)
            df    = df_runs.merge(df_bs, on="BioSample", how="left")
        else:
            print("  WARNING: No BioSample column — skipping attribute fetch.")
            df = df_runs

        merging = Path(out_csv).exists()
        print(f"\n[Step 4] {'Merging into' if merging else 'Saving →'} {out_csv} ...")
        total = save_csv_merge(df, out_csv)

        print(f"\n{'=' * 60}")
        print(f"  New runs fetched : {len(df):,}")
        print(f"  Total in file    : {total:,}  ({'merged' if merging else 'created'})  →  {out_csv}")
        print(f"{'=' * 60}")
        _print_summary(df)


# ══════════════════════════════════════════════════════════════════════════
# download subcommand
# ══════════════════════════════════════════════════════════════════════════

def download_group_iseq(
    srr_list: list[str],
    bioproject_dir: Path,
    dry_run: bool,
    threads: int,
    parallel: int,
    database: str,
    aspera: bool,
) -> None:
    """
    Download SRA runs using iSeq.

    iSeq is called with a temp SRR-list file, writing gzip FASTQ to
    <bioproject_dir>/ flat. Afterwards each SRR's files are moved into
    per-SRR subdirectories so the structure matches is_downloaded():
      <bioproject_dir>/<SRR>/<SRR>*.fastq.gz
    """
    bioproject_dir.mkdir(parents=True, exist_ok=True)

    tmp = bioproject_dir / ".iseq_input.tmp"
    tmp.write_text("\n".join(srr_list) + "\n")

    cmd = [
        "iseq",
        "-i", str(tmp),
        "-o", str(bioproject_dir),
        "-g",
        "-t", str(threads),
        "-p", str(parallel),
        "-d", database,
    ]
    if aspera:
        cmd.append("-a")

    if dry_run:
        print(f"    [DRY RUN] {' '.join(cmd)}")
        tmp.unlink(missing_ok=True)
        return

    result = subprocess.run(cmd, text=True)
    tmp.unlink(missing_ok=True)

    if result.returncode != 0:
        print(f"  [WARNING] iSeq exited with errors for {bioproject_dir.name}")

    # Reorganise flat files into per-SRR subdirectories
    for srr in srr_list:
        matched = (
            list(bioproject_dir.glob(f"{srr}*.fastq.gz")) +
            list(bioproject_dir.glob(f"{srr}*.fastq"))
        )
        if not matched:
            continue
        srr_dir = bioproject_dir / srr
        srr_dir.mkdir(exist_ok=True)
        for f in matched:
            f.rename(srr_dir / f.name)


def download_group(
    srr_list: list[str],
    bioproject_dir: Path,
    dry_run: bool,
    threads: int,
    prefetch_max_gb: int,
    no_prefetch: bool = False,
) -> None:
    """
    Files are organised into per-SRR subdirectories:
      <bioproject_dir>/<SRR>/<SRR>_1.fastq.gz

    Default (prefetch mode):
      1. prefetch all SRRs → <bioproject_dir>/<SRR>/<SRR>.sra
      2. fasterq-dump each .sra → <bioproject_dir>/<SRR>/
      3. gzip *.fastq in each SRR subdir
      4. remove .sra file

    --no-prefetch mode (bypasses SSL issues):
      1. fasterq-dump directly on each SRR accession → <bioproject_dir>/<SRR>/
      2. gzip *.fastq in each SRR subdir
    """
    bioproject_dir.mkdir(parents=True, exist_ok=True)

    # ── no-prefetch mode ──────────────────────────────────────────────
    if no_prefetch:
        for srr in srr_list:
            srr_dir = bioproject_dir / srr
            srr_dir.mkdir(exist_ok=True)

            if dry_run:
                cmd_ex = [
                    "fasterq-dump", srr,
                    "--outdir", str(srr_dir),
                    "--split-files", "--threads", str(threads),
                ]
                print(f"    [DRY RUN] {' '.join(cmd_ex)}")
                print(f"    [DRY RUN] gzip {srr_dir}/*.fastq")
                continue

            cmd_fq = [
                "fasterq-dump", srr,
                "--outdir",    str(srr_dir),
                "--split-files",
                "--threads",   str(threads),
                "--progress",
            ]
            res = subprocess.run(cmd_fq, text=True)
            if res.returncode != 0:
                print(f"  [WARNING] fasterq-dump failed for {srr}")
                continue

            _gzip_fastq(srr_dir)

        return

    # ── prefetch mode (default) ───────────────────────────────────────
    cmd_prefetch = [
        "prefetch", *srr_list,
        "--output-directory", str(bioproject_dir),
        "--max-size", f"{prefetch_max_gb}g",
        "--progress",
    ]

    if dry_run:
        print(f"    [DRY RUN] {' '.join(cmd_prefetch)}")
        sra_ex    = bioproject_dir / srr_list[0] / f"{srr_list[0]}.sra"
        cmd_fq_ex = [
            "fasterq-dump", str(sra_ex),
            "--outdir", str(bioproject_dir / srr_list[0]),
            "--split-files", "--threads", str(threads),
        ]
        print(f"    [DRY RUN] {' '.join(cmd_fq_ex)}  (× {len(srr_list)} runs)")
        print(f"    [DRY RUN] gzip <bioproject_dir>/<SRR>/*.fastq")
        print(f"    [DRY RUN] rm <bioproject_dir>/<SRR>/<SRR>.sra")
        return

    result = subprocess.run(cmd_prefetch, text=True)
    if result.returncode != 0:
        print(f"  [WARNING] prefetch exited with errors for {bioproject_dir.name}")

    for srr in srr_list:
        srr_dir  = bioproject_dir / srr
        sra_file = srr_dir / f"{srr}.sra"
        if not sra_file.exists():
            sra_file = bioproject_dir / f"{srr}.sra"
        if not sra_file.exists():
            print(f"  [WARNING] {srr}.sra not found — skipping fasterq-dump")
            continue

        srr_dir.mkdir(exist_ok=True)

        cmd_fq = [
            "fasterq-dump", str(sra_file),
            "--outdir",    str(srr_dir),
            "--split-files",
            "--threads",   str(threads),
            "--progress",
        ]
        res = subprocess.run(cmd_fq, text=True)
        if res.returncode != 0:
            print(f"  [WARNING] fasterq-dump failed for {srr}")
            continue

        sra_file.unlink(missing_ok=True)
        _gzip_fastq(srr_dir)


def _gzip_fastq(directory: Path) -> None:
    """gzip all uncompressed *.fastq files in directory."""
    fastq_files = list(directory.glob("*.fastq"))
    if fastq_files:
        result = subprocess.run(["gzip", "-1", *[str(f) for f in fastq_files]])
        if result.returncode != 0:
            print(f"  [WARNING] gzip exited with errors in {directory}")


def cmd_download(args: argparse.Namespace) -> None:
    if not args.dry_run:
        check_tools(no_prefetch=args.no_prefetch, iseq=args.iseq)

    outdir = Path(args.outdir).resolve()

    print("\n[Step 1] Loading accessions ...")
    acc_df = load_accessions_two_col(args.accessions)
    print(f"  Runs        : {len(acc_df):,}")

    groups = acc_df.groupby("BioProject")["Run"].apply(list)
    print(f"  BioProjects : {len(groups):,}")

    if args.dry_run:
        mode = "DRY RUN"
    elif args.iseq:
        mode = f"iSeq  (db={args.database}, parallel={args.parallel})"
    elif args.no_prefetch:
        mode = "fasterq-dump direct"
    else:
        mode = "prefetch + fasterq-dump"
    print(f"\n[Step 2] Downloading ({mode}) ...")

    for bp, srrs in tqdm(groups.items(), desc="  BioProject", unit="project"):
        bp_dir = outdir / bp
        tqdm.write(f"  {bp}  ({len(srrs)} runs)  →  {bp_dir}")
        if args.iseq:
            download_group_iseq(srrs, bp_dir,
                                dry_run=args.dry_run,
                                threads=args.threads,
                                parallel=args.parallel,
                                database=args.database,
                                aspera=args.aspera)
        else:
            download_group(srrs, bp_dir,
                           dry_run=args.dry_run,
                           threads=args.threads,
                           prefetch_max_gb=args.max_size,
                           no_prefetch=args.no_prefetch)

    if args.dry_run:
        print("\n[DRY RUN] Done. No files downloaded.")
        return

    print(f"\n{'=' * 55}")
    print(f"  Download complete")
    print(f"  Runs   : {len(acc_df):,}")
    print(f"  Outdir : {outdir}")
    print(f"  Next   : python MGX_data.py parse --source-data <ref.csv> --location {outdir}")
    print(f"{'=' * 55}")


# ══════════════════════════════════════════════════════════════════════════
# parse subcommand
# ══════════════════════════════════════════════════════════════════════════

def update_parsed_file(df_new: pd.DataFrame, parsed_path: str) -> None:
    """Append df_new to parsed CSV; deduplicate by Run; Zlab columns first."""
    p = Path(parsed_path)
    if p.exists():
        df_existing = pd.read_csv(p, dtype=str)
        df_combined = pd.concat([df_existing, df_new.astype(str)], ignore_index=True)
        df_combined = (df_combined
                       .drop_duplicates(subset=["Run"], keep="last")
                       .reset_index(drop=True))
    else:
        df_combined = df_new.astype(str)

    df_combined = reorder_zlab_first(df_combined)
    df_combined.to_csv(parsed_path, index=False)
    print(f"  Parsed file updated → {parsed_path}  ({len(df_combined):,} total runs)")


def update_reference(df_ref: pd.DataFrame, df_confirmed: pd.DataFrame,
                     reference_path: str) -> None:
    """Add / update Zlab columns in reference CSV; Zlab columns first; save in-place."""
    for col in ZLAB_COLS:
        if col not in df_ref.columns:
            df_ref[col] = ""

    zlab_lookup = df_confirmed.set_index("Run")[ZLAB_COLS].to_dict("index")
    mask = df_ref["Run"].isin(zlab_lookup)
    for srr, zlab_vals in zlab_lookup.items():
        idx = df_ref.index[df_ref["Run"] == srr]
        for col, val in zlab_vals.items():
            df_ref.loc[idx, col] = val

    df_ref = reorder_zlab_first(df_ref)
    df_ref.to_csv(reference_path, index=False)
    print(f"  Source data CSV updated in-place → {reference_path}")
    print(f"  {mask.sum()} runs marked as downloaded")


def cmd_parse(args: argparse.Namespace) -> None:
    location = Path(args.location).resolve()

    print("\n[Step 1] Loading source data ...")
    df_ref = pd.read_csv(args.source_data, dtype=str)
    print(f"  Source data rows : {len(df_ref):,}")

    if "Run" not in df_ref.columns:
        sys.exit("[ERROR] Source data CSV has no 'Run' column.")
    if "BioProject" not in df_ref.columns:
        sys.exit("[ERROR] Source data CSV has no 'BioProject' column.")

    ref_runs = set(df_ref["Run"])

    print("\n[Step 2] Scanning disk for unknown SRRs ...")
    disk_srrs = scan_downloaded_srrs(location)
    unknown = {srr: bp for srr, bp in disk_srrs.items() if srr not in ref_runs}
    if unknown:
        print(f"  Found {len(unknown):,} downloaded SRR(s) not in source data — back-fetching metadata ...")
        unknown_by_bp: dict[str, list[str]] = {}
        for srr, bp in unknown.items():
            unknown_by_bp.setdefault(bp, []).append(srr)
        fetched_dfs = []
        for bp, srrs in unknown_by_bp.items():
            print(f"  Fetching runinfo for {bp} ({len(srrs)} unknown run(s)) ...")
            df_bp = _fetch_sra_runinfo(bp)
            if df_bp is None:
                continue
            df_bp = df_bp[df_bp["Run"].isin(srrs)]
            if not df_bp.empty:
                fetched_dfs.append(df_bp)
        if fetched_dfs:
            df_fetched = pd.concat(fetched_dfs, ignore_index=True)
            total = save_csv_merge(df_fetched, args.source_data)
            print(f"  Source data updated → {args.source_data}  ({total:,} total rows)")
            df_ref = pd.read_csv(args.source_data, dtype=str)
            ref_runs = set(df_ref["Run"])
            print(f"  Source data rows (after update) : {len(df_ref):,}")
        else:
            print("  No metadata retrieved for unknown SRRs.")
    else:
        print(f"  No unknown SRRs found.")

    if args.accessions:
        srr_list = load_accessions_any(args.accessions)
        missing  = [s for s in srr_list if s not in ref_runs]
        if missing:
            print(f"  [WARNING] {len(missing)} accessions not found in source data:")
            for s in missing[:10]:
                print(f"    {s}")
            if len(missing) > 10:
                print(f"    ... and {len(missing)-10} more")
        srr_list = [s for s in srr_list if s in ref_runs]
        print(f"  Accessions     : {len(srr_list):,}")
    else:
        srr_list = list(ref_runs)
        print(f"  Processing all {len(srr_list):,} source data SRRs")

    if not srr_list:
        sys.exit("[ERROR] No valid accessions to process.")

    print("\n[Step 3] Scanning disk for FASTQ files ...")
    df_targets       = df_ref[df_ref["Run"].isin(srr_list)][["Run", "BioProject"]].copy()
    confirmed, not_found = [], []

    for srr in tqdm(srr_list, desc="  Scanning", unit="run"):
        rows = df_targets.loc[df_targets["Run"] == srr, "BioProject"]
        if rows.empty:
            not_found.append(srr)
            continue
        if is_downloaded(srr, location / rows.iloc[0]):
            confirmed.append(srr)
        else:
            not_found.append(srr)

    print(f"  Found     : {len(confirmed):,}")
    print(f"  Not found : {len(not_found):,}")

    if not confirmed:
        sys.exit("[ERROR] No FASTQ files found on disk. Run 'download' first.")

    print(f"\n[Step 4] Building metadata for {len(confirmed):,} confirmed runs ...")
    df_confirmed = df_ref[df_ref["Run"].isin(confirmed)].copy()
    df_confirmed = drop_blank_columns(df_confirmed)

    df_confirmed["Zlab_sort"]          = "1"
    df_confirmed["Zlab_SRA_path"]      = df_confirmed.apply(
        lambda row: str(location / row["BioProject"] / row["Run"]), axis=1
    )
    df_confirmed["Zlab_metadata_path"] = ""
    df_confirmed = reorder_zlab_first(df_confirmed)
    print(f"  Columns in parsed output : {len(df_confirmed.columns)}")

    print("\n[Step 5] Updating parsed file ...")
    update_parsed_file(df_confirmed, args.parsed)

    print("\n[Step 6] Updating source data CSV ...")
    df_ref_fresh = pd.read_csv(args.source_data, dtype=str)
    update_reference(df_ref_fresh, df_confirmed, args.source_data)

    print(f"\n{'=' * 55}")
    print(f"  Confirmed   : {len(confirmed):,} / {len(srr_list):,} runs on disk")
    print(f"  Parsed file : {args.parsed}")
    print(f"  Source data : {args.source_data}  (updated in-place)")
    print(f"{'=' * 55}")


# ══════════════════════════════════════════════════════════════════════════
# Entry point
# ══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        prog="MGX_data.py",
        description="SRA metagenomics data management — fetch, download, parse.",
    )
    sub = parser.add_subparsers(dest="command", metavar="<command>")
    sub.required = True

    # ── fetch-data ────────────────────────────────────────────────────
    fd = sub.add_parser(
        "fetch-data",
        help="Fetch SRA run and biosample metadata",
        description=(
            "Examples:\n"
            "  python MGX_data.py fetch-data --taxon 1510822 --email you@email.com\n"
            "  python MGX_data.py fetch-data --bioproject bioproject_list.txt\n"
            "  python MGX_data.py fetch-data --accession acc_list.txt\n"
            "  python MGX_data.py fetch-data --keyword \"pig gut metagenome\" "
            "--merge-into sra_taxid1510822.csv"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    fd.add_argument("--taxon",      nargs="+", default=[], metavar="ID",
                    help="NCBI Taxon ID(s) — one or more (e.g. 1510822 9606)")
    fd.add_argument("--bioproject", nargs="+", default=[], metavar="ACC/FILE",
                    help="BioProject/SRA Study accession(s), or a path to a plain-text "
                         "file with one accession per line (e.g. PRJNA857725)")
    fd.add_argument("--accession",  nargs="+", default=[], metavar="ACC/FILE",
                    help="SRA accession(s) — any type accepted: SRR, SRS, SRP, SRX, "
                         "SAMN, ERR, ERS, DRR, etc. Or a path to a plain-text file "
                         "with one accession per line.")
    fd.add_argument("--keyword",     nargs="+", default=[], metavar="QUERY",
                    help="Free-text keyword query/queries. Supports NCBI query syntax, "
                         "e.g. --keyword \"pig gut metagenome\" ")
    fd.add_argument("--merge-into", default=None, metavar="CSV",
                    help="Output CSV file for all results. If the file already exists, "
                         "new runs are appended and deduplicated on Run ID. ")
    fd.add_argument("--min-bases", type=int, default=None, metavar="N",
                    help="Exclude runs with fewer than N bases (e.g. 1000000000 for 1 Gb)")
    fd.add_argument("--max-bases", type=int, default=None, metavar="N",
                    help="Exclude runs with more than N bases (e.g. 50000000000 for 50 Gb)")
    fd.add_argument("--email",      default=None,
                    help="Your email — recommended by NCBI policy")
    fd.add_argument("--api-key",    default=None,
                    help="NCBI API key — raises rate limit 3 → 10 req/s")

    # ── download ──────────────────────────────────────────────────────
    dl = sub.add_parser(
        "download",
        help="Download SRA runs via prefetch + fasterq-dump",
        description=(
            "Download SRA runs via prefetch + fasterq-dump.\n\n"
            "Accessions file format (Run,BioProject):\n"
            "  SRR21388550,PRJNA857725\n"
            "  SRR21388551,PRJNA857725\n"
            "  (header row optional)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    dl.add_argument("--accessions", required=True,
                    help="Plain-text two-column file, no header. "
                         "Column 1: SRA run accession. Column 2: BioProject. ")
    dl.add_argument("--outdir",     required=True,
                    help="Download root directory")
    dl.add_argument("--threads",    type=int, default=FASTERQ_THREADS,
                    help=f"Threads for conversion/compression (default: {FASTERQ_THREADS})")
    dl.add_argument("--dry-run",    action="store_true",
                    help="Print commands without executing")

    # ── iSeq mode ─────────────────────────────────────────────────────
    dl.add_argument("--iseq",       action="store_true",
                    help="Use iSeq for downloading instead of SRA Toolkit. "
                         "Recommended: faster, supports ENA mirror, resumable. "
                         "Requires: conda install bioconda::iseq")
    dl.add_argument("--parallel",   type=int, default=8, metavar="N",
                    help="iSeq: parallel download connections (default: 8). "
                         "Only used with --iseq.")
    dl.add_argument("--database",   default="ena", choices=["ena", "sra"],
                    help="iSeq: download source database (default: ena). "
                         "ena is generally faster outside the US. "
                         "Only used with --iseq.")
    dl.add_argument("--aspera",     action="store_true",
                    help="iSeq: use Aspera for faster transfers. "
                         "Only used with --iseq.")

    # ── SRA Toolkit mode ──────────────────────────────────────────────
    dl.add_argument("--max-size",   type=int, default=PREFETCH_MAX_GB,
                    help=f"Max SRA file size in GB for prefetch (default: {PREFETCH_MAX_GB}). "
                         "Only used without --iseq.")
    dl.add_argument("--no-prefetch", action="store_true",
                    help="Skip prefetch; call fasterq-dump directly on each SRR accession. "
                         "Use this if prefetch fails with SSL/TLS errors. "
                         "Only used without --iseq.")

    # ── parse ─────────────────────────────────────────────────────────
    pa = sub.add_parser(
        "parse",
        help="Scan downloaded FASTQ files and update metadata CSVs",
        description=(
            "Scan downloaded FASTQ files on disk, build curated metadata,\n"
            "and update parsed CSV + source data CSV in-place.\n\n"
            "Zlab columns added (first 3 columns of both output CSVs):\n"
            "  Zlab_sort           manual sort order (default '1')\n"
            "  Zlab_SRA_path       absolute path to BioProject FASTQ directory\n"
            "  Zlab_metadata_path  blank — fill manually if needed"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    pa.add_argument("--source-data", required=True, dest="source_data",
                    help="Source metadata CSV from fetch-data. "
                         "Updated in-place with Zlab columns and any back-fetched runs.")
    pa.add_argument("--location",   required=True,
                    help="Root directory containing BioProject subdirectories "
                         "(e.g. <location>/<BioProject>/<SRR>/)")
    pa.add_argument("--parsed",     default="parsed_metadata.csv",
                    help="Curated output CSV, appended across batches "
                         "(default: parsed_metadata.csv)")
    pa.add_argument("--accessions",
                    help="Optional: limit scan to these SRRs. Accepts single-column "
                         "SRR list or two-column Run,BioProject file. "
                         "Omit to scan all SRRs in the source data.")

    args = parser.parse_args()

    dispatch = {
        "fetch-data": cmd_fetch_data,
        "download":   cmd_download,
        "parse":      cmd_parse,
    }
    dispatch[args.command](args)


if __name__ == "__main__":
    main()
