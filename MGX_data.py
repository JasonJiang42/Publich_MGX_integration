#!/usr/bin/env python3
"""
MGX_data.py
────────────
Four subcommands for SRA metagenomics data management:

  fetch-taxon    Fetch all WGS/METAGENOMIC SRA runs for a taxon ID.
  fetch-project  Fetch full SRA run metadata for one BioProject/Study.
  download       Download SRA runs via prefetch + fasterq-dump.
  parse          Scan downloaded FASTQ files and update metadata CSVs.

Usage:
    python MGX_data.py fetch-taxon   --taxon-id 1510822 --email you@email.com
    python MGX_data.py fetch-project --accession PRJNA857725 --email you@email.com
    python MGX_data.py download      --accessions runs.csv --outdir ./fastq
    python MGX_data.py parse         --reference sra_taxid1510822.csv --outdir ./fastq

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

EFETCH_URL       = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESEARCH_URL      = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
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
    Load a two-column accessions file (Run, BioProject).
    Accepts comma or tab separation, with or without a header row.
    """
    for sep in (",", "\t"):
        df = pd.read_csv(path, sep=sep, dtype=str, header=None)
        if df.shape[1] >= 2:
            break
    else:
        sys.exit(f"[ERROR] Cannot parse {path} as two-column.\n"
                 "  Expected: Run,BioProject (comma or tab separated)")

    df = df.iloc[:, :2].copy()
    df.columns = ["Run", "BioProject"]
    first = str(df.iloc[0, 0]).strip().upper()
    if not (first.startswith("SRR") or first.startswith("ERR") or first.startswith("DRR")):
        df = df.iloc[1:].reset_index(drop=True)
    return df.dropna().reset_index(drop=True)


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


def is_downloaded(srr: str, bioproject_dir: Path) -> bool:
    """Return True if any FASTQ file exists in <bioproject_dir>/<srr>/."""
    srr_dir = bioproject_dir / srr
    return srr_dir.is_dir() and bool(
        list(srr_dir.glob(f"{srr}*.fastq.gz")) or
        list(srr_dir.glob(f"{srr}*.fastq"))
    )


def save_csv_merge(df_new: pd.DataFrame, path: str) -> int:
    """
    Write df_new to path. If path already exists, merge by outer-joining
    columns and deduplicating on 'Run' (keep latest). Returns total row count.
    """
    p = Path(path)
    if p.exists():
        df_existing = pd.read_csv(p, dtype=str)
        df_combined = pd.concat([df_existing, df_new.astype(str)], ignore_index=True)
        df_combined = (df_combined
                       .drop_duplicates(subset=["Run"], keep="last")
                       .reset_index(drop=True))
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


def check_tools(no_prefetch: bool = False) -> None:
    tools = ["fasterq-dump", "gzip"] if no_prefetch else ["prefetch", "fasterq-dump", "gzip"]
    missing = [t for t in tools
               if subprocess.run(["which", t], capture_output=True).returncode != 0]
    if missing:
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


def _fetch_runinfo(uids: list[str], checkpoint_file: str,
                   api_key: str | None, sleep: float) -> pd.DataFrame:
    """Fetch SRA runinfo in batches, filter to WGS/METAGENOMIC, checkpoint/resume."""
    rows: list[dict] = []
    start_idx = 0

    if os.path.exists(checkpoint_file):
        with open(checkpoint_file) as fh:
            ckpt = json.load(fh)
        rows      = ckpt["rows"]
        start_idx = ckpt["next_idx"]
        print(f"\n  Resuming from UID index {start_idx:,} "
              f"({len(rows):,} runs already collected)")

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
            df_batch = pd.read_csv(io.StringIO(resp.text)).dropna(how="all")
            mask     = df_batch["LibraryStrategy"].str.upper().isin(TARGET_STRATEGIES)
            df_batch = df_batch[mask]
            if not df_batch.empty:
                rows.extend(df_batch.to_dict("records"))
        except Exception as exc:
            tqdm.write(f"  WARNING: batch at index {i} failed — {exc}")

        time.sleep(sleep)
        pbar.update(len(batch))

        with open(checkpoint_file, "w") as fh:
            json.dump({"rows": rows, "next_idx": i + RUNINFO_BATCH}, fh)

    pbar.close()

    if os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)

    if not rows:
        sys.exit("No WGS/METAGENOMIC runs found.")

    df = pd.DataFrame(rows).drop_duplicates(subset=["Run"]).reset_index(drop=True)
    print(f"  WGS/METAGENOMIC runs: {len(df):,}")
    return df


def cmd_fetch_taxon(args: argparse.Namespace) -> None:
    try:
        from Bio import Entrez  # noqa: F401 — only needed for api_key plumbing
    except ImportError:
        sys.exit("[ERROR] biopython not installed.\n  pip install biopython")

    api_key = args.api_key or None
    sleep   = 0.11 if api_key else 0.4
    label   = args.label or f"sra_taxid{args.taxon_id}"
    out_csv = f"{label}.csv"
    ckpt    = f"{label}_checkpoint.json"

    print("=" * 60)
    print(f"  fetch-taxon  —  taxon {args.taxon_id}")
    print(f"  strategies   :  {', '.join(sorted(TARGET_STRATEGIES))}")
    print("=" * 60)

    uids    = _search_sra_uids(args.taxon_id, api_key)
    df_runs = _fetch_runinfo(uids, ckpt, api_key, sleep)

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

    print(f"\nLibraryStrategy:")
    for val, n in df["LibraryStrategy"].value_counts().items():
        print(f"  {val:<22} {n:>6}  {'█' * (n * 30 // max(len(df), 1))}")

    print(f"\nPlatform:")
    for val, n in df["Platform"].value_counts().items():
        print(f"  {val:<22} {n:>6}")

    total_gb = pd.to_numeric(df["bases"], errors="coerce").sum() / 1e9
    print(f"\nTotal data: {total_gb:.1f} Gb")


# ══════════════════════════════════════════════════════════════════════════
# fetch-project subcommand
# ══════════════════════════════════════════════════════════════════════════

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
        sys.exit(f"  ERROR: No data returned for {accession}. "
                 "Check the accession is valid and has public data.")

    df = pd.read_csv(io.StringIO(content)).dropna(how="all").reset_index(drop=True)
    print(f"  Runs: {len(df):,}")
    return df


def cmd_fetch_project(args: argparse.Namespace) -> None:
    accession = args.accession.strip().upper()
    out_csv   = args.out or f"{accession}_sra_runs.csv"

    print("=" * 60)
    print(f"  fetch-project  —  {accession}")
    print("=" * 60)

    df_runs = _fetch_sra_runinfo(accession)

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

    print(f"\nLibraryStrategy:")
    for val, n in df["LibraryStrategy"].value_counts().items():
        print(f"  {val:<25} {n:>5}  {'█' * (n * 30 // max(len(df), 1))}")

    print(f"\nPlatform:")
    for val, n in df["Platform"].value_counts().items():
        print(f"  {val:<25} {n:>5}")

    total_gb = pd.to_numeric(df["bases"], errors="coerce").sum() / 1e9
    print(f"\nTotal data: {total_gb:.1f} Gb")

    sra_cols = set(df_runs.columns)
    bs_cols  = [c for c in df.columns if c not in sra_cols]
    if bs_cols:
        print(f"\nBioSample attributes added ({len(bs_cols)}): {', '.join(bs_cols)}")


# ══════════════════════════════════════════════════════════════════════════
# download subcommand
# ══════════════════════════════════════════════════════════════════════════

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
        check_tools(no_prefetch=args.no_prefetch)

    outdir = Path(args.outdir).resolve()

    print("\n[Step 1] Loading accessions ...")
    acc_df = load_accessions_two_col(args.accessions)
    print(f"  Runs        : {len(acc_df):,}")

    groups = acc_df.groupby("BioProject")["Run"].apply(list)
    print(f"  BioProjects : {len(groups):,}")

    mode = "DRY RUN" if args.dry_run else ("fasterq-dump direct" if args.no_prefetch else "prefetch + fasterq-dump")
    print(f"\n[Step 2] Downloading ({mode}) ...")

    for bp, srrs in tqdm(groups.items(), desc="  BioProject", unit="project"):
        bp_dir = outdir / bp
        tqdm.write(f"  {bp}  ({len(srrs)} runs)  →  {bp_dir}")
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
    print(f"  Next   : python MGX_data.py parse --reference <ref.csv> --outdir {outdir}")
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
    print(f"  Reference CSV updated in-place → {reference_path}")
    print(f"  {mask.sum()} runs marked as downloaded")


def cmd_parse(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir).resolve()

    print("\n[Step 1] Loading reference ...")
    df_ref = pd.read_csv(args.reference, dtype=str)
    print(f"  Reference rows : {len(df_ref):,}")

    if "Run" not in df_ref.columns:
        sys.exit("[ERROR] Reference CSV has no 'Run' column.")
    if "BioProject" not in df_ref.columns:
        sys.exit("[ERROR] Reference CSV has no 'BioProject' column.")

    ref_runs = set(df_ref["Run"])

    if args.accessions:
        srr_list = load_accessions_any(args.accessions)
        missing  = [s for s in srr_list if s not in ref_runs]
        if missing:
            print(f"  [WARNING] {len(missing)} accessions not found in reference:")
            for s in missing[:10]:
                print(f"    {s}")
            if len(missing) > 10:
                print(f"    ... and {len(missing)-10} more")
        srr_list = [s for s in srr_list if s in ref_runs]
        print(f"  Accessions     : {len(srr_list):,}")
    else:
        srr_list = list(ref_runs)
        print(f"  Scanning all {len(srr_list):,} reference SRRs")

    if not srr_list:
        sys.exit("[ERROR] No valid accessions to process.")

    print("\n[Step 2] Scanning disk for FASTQ files ...")
    df_targets       = df_ref[df_ref["Run"].isin(srr_list)][["Run", "BioProject"]].copy()
    confirmed, not_found = [], []

    for srr in tqdm(srr_list, desc="  Scanning", unit="run"):
        rows = df_targets.loc[df_targets["Run"] == srr, "BioProject"]
        if rows.empty:
            not_found.append(srr)
            continue
        if is_downloaded(srr, outdir / rows.iloc[0]):
            confirmed.append(srr)
        else:
            not_found.append(srr)

    print(f"  Found     : {len(confirmed):,}")
    print(f"  Not found : {len(not_found):,}")

    if not confirmed:
        sys.exit("[ERROR] No FASTQ files found on disk. Run 'download' first.")

    print(f"\n[Step 3] Building metadata for {len(confirmed):,} confirmed runs ...")
    df_confirmed = df_ref[df_ref["Run"].isin(confirmed)].copy()
    df_confirmed = drop_blank_columns(df_confirmed)

    df_confirmed["Zlab_sort"]          = "1"
    df_confirmed["Zlab_SRA_path"]      = df_confirmed.apply(
        lambda row: str(outdir / row["BioProject"] / row["Run"]), axis=1
    )
    df_confirmed["Zlab_metadata_path"] = ""
    df_confirmed = reorder_zlab_first(df_confirmed)
    print(f"  Columns in parsed output : {len(df_confirmed.columns)}")

    print("\n[Step 4] Updating parsed file ...")
    update_parsed_file(df_confirmed, args.parsed)

    print("\n[Step 5] Updating reference CSV ...")
    df_ref_fresh = pd.read_csv(args.reference, dtype=str)
    update_reference(df_ref_fresh, df_confirmed, args.reference)

    print(f"\n{'=' * 55}")
    print(f"  Confirmed   : {len(confirmed):,} / {len(srr_list):,} runs on disk")
    print(f"  Parsed file : {args.parsed}")
    print(f"  Reference   : {args.reference}  (updated in-place)")
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

    # ── fetch-taxon ───────────────────────────────────────────────────
    ft = sub.add_parser(
        "fetch-taxon",
        help="Fetch all WGS/METAGENOMIC SRA runs for a taxon ID",
        description=(
            "Search SRA by NCBI Taxon ID and fetch all WGS/METAGENOMIC run metadata,\n"
            "merged with BioSample attributes. Supports checkpoint/resume."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ft.add_argument("--taxon-id", required=True,
                    help="NCBI Taxon ID (e.g. 1510822 for pig gut metagenome)")
    ft.add_argument("--email",    required=True,
                    help="Your email — required by NCBI policy")
    ft.add_argument("--api-key",  default=None,
                    help="NCBI API key — raises rate limit 3 → 10 req/s")
    ft.add_argument("--label",    default=None,
                    help="Output filename prefix (default: sra_taxid<ID>)")

    # ── fetch-project ─────────────────────────────────────────────────
    fp = sub.add_parser(
        "fetch-project",
        help="Fetch full SRA run metadata for one BioProject or SRA Study",
        description=(
            "Fetch SRA Run Selector-style metadata for a single BioProject\n"
            "or SRA Study accession, including all BioSample attributes."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    fp.add_argument("--accession", required=True,
                    help="BioProject (PRJNA…) or SRA Study (SRP…) accession")
    fp.add_argument("--email",     default=None,
                    help="Your email (recommended for NCBI requests)")
    fp.add_argument("--out",       default=None,
                    help="Output CSV filename (default: <accession>_sra_runs.csv)")

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
                    help="Two-column CSV/TSV: Run, BioProject (header optional)")
    dl.add_argument("--outdir",     required=True,
                    help="Download root directory")
    dl.add_argument("--threads",    type=int, default=FASTERQ_THREADS,
                    help=f"fasterq-dump threads (default: {FASTERQ_THREADS})")
    dl.add_argument("--max-size",   type=int, default=PREFETCH_MAX_GB,
                    help=f"Max SRA file size in GB for prefetch (default: {PREFETCH_MAX_GB})")
    dl.add_argument("--dry-run",    action="store_true",
                    help="Print commands without executing")
    dl.add_argument("--no-prefetch", action="store_true",
                    help="Skip prefetch; call fasterq-dump directly on each SRR accession. "
                         "Use this if prefetch fails with SSL/TLS errors (common in mainland China)")

    # ── parse ─────────────────────────────────────────────────────────
    pa = sub.add_parser(
        "parse",
        help="Scan downloaded FASTQ files and update metadata CSVs",
        description=(
            "Scan downloaded FASTQ files on disk, build curated metadata,\n"
            "and update parsed CSV + reference CSV in-place.\n\n"
            "Zlab columns added (first 3 columns of both output CSVs):\n"
            "  Zlab_sort           manual sort order (default '1')\n"
            "  Zlab_SRA_path       absolute path to BioProject FASTQ directory\n"
            "  Zlab_metadata_path  blank — fill manually if needed"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    pa.add_argument("--reference",  required=True,
                    help="Reference metadata CSV from fetch-taxon or fetch-project")
    pa.add_argument("--outdir",     required=True,
                    help="Root directory containing BioProject subdirectories")
    pa.add_argument("--parsed",     default="parsed_metadata.csv",
                    help="Curated output CSV, appended across batches "
                         "(default: parsed_metadata.csv)")
    pa.add_argument("--accessions",
                    help="Optional: limit scan to these SRRs. Accepts single-column "
                         "SRR list or two-column Run,BioProject file. "
                         "Omit to scan all SRRs in the reference.")

    args = parser.parse_args()

    dispatch = {
        "fetch-taxon":   cmd_fetch_taxon,
        "fetch-project": cmd_fetch_project,
        "download":      cmd_download,
        "parse":         cmd_parse,
    }
    dispatch[args.command](args)


if __name__ == "__main__":
    main()
