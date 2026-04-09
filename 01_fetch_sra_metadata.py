#!/usr/bin/env python3
"""
01_fetch_sra_metadata.py
─────────────────────────
Fetch all SRA run metadata for a BioProject accession,
matching the full metadata table from the NCBI SRA Run Selector page
(https://www.ncbi.nlm.nih.gov/Traces/study/).

Includes BioSample attributes: geo_loc_name, host, collection_date, etc.

Steps:
  1. Fetch SRA runinfo (Run, BioSample, LibraryStrategy, spots, bases ...)
  2. Fetch BioSample XML for each unique BioSample accession
  3. Merge → output full table matching SRA Run Selector

Usage:
    python 01_fetch_sra_metadata.py --accession PRJNA857725 --email your@email.com

Requirements:
    pip install requests pandas
"""

from __future__ import annotations

import argparse
import io
import sys
import time
import xml.etree.ElementTree as ET

import pandas as pd
import requests

# ── endpoints ──────────────────────────────────────────────────────────────
SRA_RUNS_URL   = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runs"
BIOSAMPLE_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

SLEEP          = 0.4   # seconds between requests (NCBI free tier = 3 req/s)
BATCH_SIZE     = 100   # BioSamples per efetch call


# ══════════════════════════════════════════════════════════════════════════
# Step 1 — SRA runinfo
# ══════════════════════════════════════════════════════════════════════════

def fetch_sra_runinfo(accession: str) -> pd.DataFrame:
    """
    Fetch SRA runinfo from the SRA Run Selector backend.
    Returns a DataFrame with ~46 columns per run.
    """
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
                 f"Check the accession is valid and has public data.")

    df = pd.read_csv(io.StringIO(content)).dropna(how="all").reset_index(drop=True)
    print(f"  Runs: {len(df)}")
    return df


# ══════════════════════════════════════════════════════════════════════════
# Step 2 — BioSample attributes
# ══════════════════════════════════════════════════════════════════════════

def fetch_biosample_attributes(biosample_accs: list[str]) -> pd.DataFrame:
    """
    Fetch BioSample XML for each accession and extract all sample attributes
    (geo_loc_name, host, collection_date, isolation_source, etc.).
    Returns a DataFrame keyed on BioSample accession.
    """
    print(f"\n[Step 2] Fetching BioSample attributes for "
          f"{len(biosample_accs)} unique samples ...")

    all_records: list[dict] = []
    n_batches = (len(biosample_accs) + BATCH_SIZE - 1) // BATCH_SIZE

    for i in range(0, len(biosample_accs), BATCH_SIZE):
        batch     = biosample_accs[i : i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        print(f"  Batch {batch_num}/{n_batches} ({len(batch)} samples) ...", end=" ")

        for attempt in range(5):
            try:
                resp = requests.get(
                    BIOSAMPLE_URL,
                    params={"db": "biosample", "id": ",".join(batch),
                            "rettype": "xml", "retmode": "xml"},
                    timeout=60,
                )
                resp.raise_for_status()
                break
            except requests.RequestException as exc:
                if attempt == 4:
                    print(f"\n  ERROR on batch {batch_num}: {exc}")
                    break
                wait = 5 * (2 ** attempt)
                print(f"\n  [retry {attempt+1}/5] {exc} — waiting {wait}s")
                time.sleep(wait)

        time.sleep(SLEEP)

        # Parse XML — each <BioSample> has <Attributes><Attribute .../></Attributes>
        try:
            root = ET.fromstring(resp.text)
        except ET.ParseError as exc:
            print(f"XML parse error: {exc}")
            continue

        n = 0
        for bs in root.findall(".//BioSample"):
            acc = bs.get("accession", "")
            if not acc:
                continue
            record: dict = {"BioSample": acc}

            # Extract all attributes
            for attr in bs.findall(".//Attribute"):
                name = (attr.get("attribute_name") or
                        attr.get("harmonized_name") or "").strip()
                if name:
                    record[name] = (attr.text or "").strip()

            all_records.append(record)
            n += 1

        print(f"{n} parsed")

    if not all_records:
        print("  WARNING: No BioSample attributes retrieved.")
        return pd.DataFrame({"BioSample": biosample_accs})

    df_bs = pd.DataFrame(all_records)
    print(f"  Total BioSample attributes columns: {len(df_bs.columns) - 1}")
    return df_bs


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Fetch full SRA Run Selector metadata (incl. BioSample attributes)."
    )
    parser.add_argument("--accession", required=True,
                        help="BioProject accession (e.g. PRJNA857725) "
                             "or SRA Study accession (e.g. SRP395053)")
    parser.add_argument("--email",     required=False, default=None,
                        help="Your email (recommended for NCBI requests)")
    parser.add_argument("--out",       default=None,
                        help="Output CSV filename "
                             "(default: <accession>_sra_runs.csv)")
    args = parser.parse_args()

    accession = args.accession.strip().upper()
    out_csv   = args.out or f"{accession}_sra_runs.csv"

    print("=" * 60)
    print(f"  01_fetch_sra_metadata.py  —  {accession}")
    print("=" * 60)

    # ── Step 1: SRA runinfo ───────────────────────────────────────────
    df_runs = fetch_sra_runinfo(accession)

    # ── Step 2: BioSample attributes ─────────────────────────────────
    biosample_accs = (
        df_runs["BioSample"].dropna().unique().tolist()
        if "BioSample" in df_runs.columns else []
    )

    if biosample_accs:
        df_bs = fetch_biosample_attributes(biosample_accs)
        # Merge BioSample attributes into run table
        df = df_runs.merge(df_bs, on="BioSample", how="left")
    else:
        print("  WARNING: No BioSample column found — skipping attribute fetch.")
        df = df_runs

    # ── Save ──────────────────────────────────────────────────────────
    print(f"\n[Step 3] Saving to {out_csv} ...")
    df.to_csv(out_csv, index=False)

    # ── Summary ───────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print(f"  Done — {len(df):,} runs  x  {len(df.columns)} columns → {out_csv}")
    print(f"{'=' * 60}")

    # BioSample attribute columns added
    sra_cols = set(df_runs.columns)
    new_cols  = [c for c in df.columns if c not in sra_cols]
    if new_cols:
        print(f"\nBioSample attributes added ({len(new_cols)}):")
        print(f"  {', '.join(new_cols)}")

    print(f"\nLibraryStrategy:")
    for val, n in df["LibraryStrategy"].value_counts().items():
        bar = "█" * (n * 30 // max(len(df), 1))
        print(f"  {val:<25} {n:>5}  {bar}")

    print(f"\nPlatform:")
    for val, n in df["Platform"].value_counts().items():
        print(f"  {val:<25} {n:>5}")

    total_gb = pd.to_numeric(df["bases"], errors="coerce").sum() / 1e9
    print(f"\nTotal data size : {total_gb:.1f} Gb")


if __name__ == "__main__":
    main()
