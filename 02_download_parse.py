#!/usr/bin/env python3
"""
02_download_parse.py
─────────────────────
Download SRA runs and parse metadata. Designed for multi-batch use —
the parsed output file is safely appended across multiple runs.

Input:
  --accessions   text file with one SRR accession per line
  --reference    reference metadata CSV from sra_fetch.py
  --outdir       download root directory
  --parsed       parsed output CSV (appended in-place across batches)
  --curator      curator name recorded in Zlab_curator column

Download structure:
  <outdir>/<BioProject>/SRR*.fastq.gz

Steps:
  1  Load reference CSV + accession list
  2  Group SRRs by BioProject, download via iSeq
  3  Verify downloaded files exist on disk
  4  Extract metadata for confirmed downloads
  5  Drop all-blank columns, add 5 Zlab columns
  6  Append to parsed CSV (multi-batch safe, deduplicated by Run)
  7  Update reference CSV in-place with Zlab columns

Usage:
    python 02_download_parse.py \\
        --accessions runs.txt \\
        --reference  sra_taxid1510822.csv \\
        --outdir     ./fastq \\
        --parsed     parsed_metadata.csv \\
        --curator    "Jason"

    # Preview without downloading
    python 02_download_parse.py ... --dry-run

Requirements:
    pip install pandas tqdm
    iSeq must be on PATH (https://github.com/BioOmics/iSeq)
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import date
from pathlib import Path

import pandas as pd
from tqdm import tqdm

# ── iSeq defaults ──────────────────────────────────────────────────────────
ISEQ_THREADS  = 8
ISEQ_PARALLEL = 4
ISEQ_DATABASE = "ena"   # ENA mirror — faster in mainland China


# ══════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════

def load_accession_list(path: str) -> list[str]:
    with open(path) as fh:
        accs = [line.strip() for line in fh if line.strip()]
    return accs


def check_iseq() -> None:
    if subprocess.run(["which", "iseq"], capture_output=True).returncode != 0:
        sys.exit(
            "[ERROR] iSeq not found on PATH.\n"
            "  Install: wget https://github.com/BioOmics/iSeq/releases/download/"
            "v1.9.8/iSeq-v1.9.8.tar.gz && tar -zxf iSeq-v1.9.8.tar.gz\n"
            "  Then add iSeq-v1.9.8/bin to PATH."
        )


def is_downloaded(srr: str, bioproject_dir: Path) -> bool:
    """Return True if any FASTQ file for this SRR exists in the directory."""
    return bool(list(bioproject_dir.glob(f"{srr}*.fastq*")))


def drop_blank_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Drop columns where ALL values are NaN or empty string."""
    def all_blank(col):
        return df[col].isna().all() or df[col].astype(str).str.strip().eq("").all()
    keep = [c for c in df.columns if not all_blank(c)]
    return df[keep]


# ══════════════════════════════════════════════════════════════════════════
# Step 2 — Download via iSeq
# ══════════════════════════════════════════════════════════════════════════

def download_group(
    srr_list: list[str],
    bioproject_dir: Path,
    dry_run: bool,
    aspera: bool,
) -> None:
    """Write a temp accession file and call iSeq for one BioProject group."""
    bioproject_dir.mkdir(parents=True, exist_ok=True)
    acc_file = bioproject_dir / "_accessions.txt"
    acc_file.write_text("\n".join(srr_list) + "\n")

    cmd = [
        "iseq",
        "-i", str(acc_file),
        "-o", str(bioproject_dir),
        "-d", ISEQ_DATABASE,
        "-t", str(ISEQ_THREADS),
        "-p", str(ISEQ_PARALLEL),
        "-g",   # gzip output
    ]
    if aspera:
        cmd.append("-a")

    if dry_run:
        print(f"    [DRY RUN] {' '.join(cmd)}")
        return

    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        print(f"  [WARNING] iSeq exited with errors for {bioproject_dir.name}")


# ══════════════════════════════════════════════════════════════════════════
# Step 6 — Append to parsed CSV (multi-batch safe)
# ══════════════════════════════════════════════════════════════════════════

def update_parsed_file(df_new: pd.DataFrame, parsed_path: str) -> None:
    """
    Append df_new to the existing parsed CSV.
    - If file exists: outer-join columns, deduplicate by Run (keep latest).
    - If not exists: create fresh.
    """
    p = Path(parsed_path)
    if p.exists():
        df_existing = pd.read_csv(p, dtype=str)
        df_combined = pd.concat([df_existing, df_new.astype(str)],
                                 ignore_index=True)
        # Keep latest record for each Run
        df_combined = (df_combined
                       .drop_duplicates(subset=["Run"], keep="last")
                       .reset_index(drop=True))
    else:
        df_combined = df_new.astype(str)

    df_combined.to_csv(parsed_path, index=False)
    print(f"  Parsed file updated → {parsed_path}  ({len(df_combined):,} total runs)")


# ══════════════════════════════════════════════════════════════════════════
# Step 7 — Update reference CSV in-place
# ══════════════════════════════════════════════════════════════════════════

ZLAB_COLS = [
    "Zlab_sort",
    "Zlab_SRA_path",
    "Zlab_metadata_path",
    "Zlab_curator",
    "Zlab_curated_date",
]

def update_reference(df_ref: pd.DataFrame, df_confirmed: pd.DataFrame,
                     reference_path: str) -> None:
    """
    Add Zlab columns to reference CSV if absent, then fill values
    for confirmed downloaded runs. Saves in-place.
    """
    for col in ZLAB_COLS:
        if col not in df_ref.columns:
            df_ref[col] = ""

    # Build lookup: Run → Zlab values from confirmed batch
    zlab_lookup = df_confirmed.set_index("Run")[ZLAB_COLS].to_dict("index")

    mask = df_ref["Run"].isin(zlab_lookup)
    for srr, zlab_vals in zlab_lookup.items():
        idx = df_ref.index[df_ref["Run"] == srr]
        for col, val in zlab_vals.items():
            df_ref.loc[idx, col] = val

    df_ref.to_csv(reference_path, index=False)
    print(f"  Reference CSV updated in-place → {reference_path}")
    print(f"  {mask.sum()} runs marked as downloaded")


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Download SRA runs via iSeq and parse metadata."
    )
    parser.add_argument("--accessions", required=True,
                        help="Text file with one SRR accession per line")
    parser.add_argument("--reference",  required=True,
                        help="Reference metadata CSV from sra_fetch.py")
    parser.add_argument("--outdir",     required=True,
                        help="Download root directory")
    parser.add_argument("--parsed",     default="parsed_metadata.csv",
                        help="Parsed output CSV, appended across batches "
                             "(default: parsed_metadata.csv)")
    parser.add_argument("--curator",    required=True,
                        help="Curator name for Zlab_curator column")
    parser.add_argument("--aspera",     action="store_true",
                        help="Use Aspera for faster download (requires ascp)")
    parser.add_argument("--dry-run",    action="store_true",
                        help="Show iSeq commands without downloading")
    args = parser.parse_args()

    if not args.dry_run:
        check_iseq()

    outdir = Path(args.outdir)
    today  = date.today().isoformat()   # e.g. 2026-04-09

    # ── Step 1: Load inputs ───────────────────────────────────────────
    print("\n[Step 1] Loading inputs ...")

    srr_list = load_accession_list(args.accessions)
    print(f"  Accessions loaded : {len(srr_list):,}")

    df_ref = pd.read_csv(args.reference, dtype=str)
    print(f"  Reference rows    : {len(df_ref):,}")

    if "Run" not in df_ref.columns:
        sys.exit("[ERROR] Reference CSV has no 'Run' column.")
    if "BioProject" not in df_ref.columns:
        sys.exit("[ERROR] Reference CSV has no 'BioProject' column.")

    # Check all SRRs exist in reference
    ref_runs  = set(df_ref["Run"])
    missing   = [s for s in srr_list if s not in ref_runs]
    if missing:
        print(f"  [WARNING] {len(missing)} accessions not found in reference CSV:")
        for s in missing[:10]:
            print(f"    {s}")
        if len(missing) > 10:
            print(f"    ... and {len(missing)-10} more")

    # Only proceed with accessions present in reference
    srr_list = [s for s in srr_list if s in ref_runs]
    if not srr_list:
        sys.exit("[ERROR] No valid accessions to process.")

    # Group by BioProject
    df_targets = df_ref[df_ref["Run"].isin(srr_list)][["Run", "BioProject"]].copy()
    groups     = df_targets.groupby("BioProject")["Run"].apply(list)
    print(f"  BioProjects       : {len(groups):,}")

    # ── Step 2: Download ──────────────────────────────────────────────
    print(f"\n[Step 2] Downloading via iSeq "
          f"({'DRY RUN' if args.dry_run else ISEQ_DATABASE.upper()}) ...")

    for bp, srrs in tqdm(groups.items(), desc="  BioProject", unit="project"):
        bp_dir = outdir / bp
        tqdm.write(f"  {bp}  ({len(srrs)} runs)  →  {bp_dir}")
        download_group(srrs, bp_dir, dry_run=args.dry_run, aspera=args.aspera)

    if args.dry_run:
        print("\n[DRY RUN] Skipping verification and metadata steps.")
        return

    # ── Step 3: Verify downloads ──────────────────────────────────────
    print("\n[Step 3] Verifying downloads ...")

    confirmed, failed = [], []
    for srr in srr_list:
        bp    = df_targets.loc[df_targets["Run"] == srr, "BioProject"].iloc[0]
        bp_dir = outdir / bp
        if is_downloaded(srr, bp_dir):
            confirmed.append(srr)
        else:
            failed.append(srr)

    print(f"  Confirmed : {len(confirmed):,}")
    print(f"  Failed    : {len(failed):,}")
    if failed:
        failed_path = outdir / "failed_runs.txt"
        failed_path.write_text("\n".join(failed) + "\n")
        print(f"  Failed list → {failed_path}")

    if not confirmed:
        sys.exit("[ERROR] No runs confirmed downloaded. Check iSeq output.")

    # ── Step 4: Extract metadata for confirmed runs ───────────────────
    print("\n[Step 4] Building metadata for confirmed downloads ...")

    df_confirmed = df_ref[df_ref["Run"].isin(confirmed)].copy()

    # ── Step 5: Drop all-blank columns + add Zlab columns ────────────
    df_confirmed = drop_blank_columns(df_confirmed)

    df_confirmed["Zlab_sort"]          = "1"
    df_confirmed["Zlab_SRA_path"]      = df_confirmed["BioProject"].apply(
        lambda bp: str(outdir / bp)
    )
    df_confirmed["Zlab_metadata_path"] = ""
    df_confirmed["Zlab_curator"]       = args.curator
    df_confirmed["Zlab_curated_date"]  = today

    print(f"  Columns in parsed output : {len(df_confirmed.columns)}")

    # ── Step 6: Update parsed file ────────────────────────────────────
    print("\n[Step 5] Updating parsed file ...")
    update_parsed_file(df_confirmed, args.parsed)

    # ── Step 7: Update reference CSV in-place ────────────────────────
    print("\n[Step 6] Updating reference CSV ...")

    # Reload reference to get fresh copy (avoids dtype issues)
    df_ref_fresh = pd.read_csv(args.reference, dtype=str)
    update_reference(df_ref_fresh, df_confirmed, args.reference)

    # ── Summary ───────────────────────────────────────────────────────
    print(f"\n{'=' * 55}")
    print(f"  Batch complete")
    print(f"  Downloaded  : {len(confirmed):,} / {len(srr_list):,} runs")
    print(f"  Parsed file : {args.parsed}")
    print(f"  Reference   : {args.reference}  (updated in-place)")
    if failed:
        print(f"  Failed runs : {outdir}/failed_runs.txt")
    print(f"{'=' * 55}")


if __name__ == "__main__":
    main()
