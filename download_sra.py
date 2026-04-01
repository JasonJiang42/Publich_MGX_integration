#!/usr/bin/env python3
"""
download_sra.py — Download METAGENOMIC and WGS SRA FASTQ files for BioProject accessions.

Usage:
    python download_sra.py --accessions PRJNA123456 PRJNA789012 --email you@example.com
    python download_sra.py --accessions PRJNA123456 --outdir ./fastq --threads 6
    python download_sra.py --accessions PRJNA123456 --dry-run
"""

import argparse
import csv
import http.client
import io
import subprocess
import sys
import time
import urllib.error
import xml.etree.ElementTree as ET
from pathlib import Path

import requests
from Bio import Entrez

# ---------------------------------------------------------------------------
# Entrez helpers
# ---------------------------------------------------------------------------

def entrez_call(fn, *args, max_attempts: int = 5, **kwargs):
    """Call a Bio.Entrez function with exponential backoff."""
    delay = 5
    for attempt in range(1, max_attempts + 1):
        try:
            return fn(*args, **kwargs)
        except (RuntimeError,
                urllib.error.HTTPError,
                urllib.error.URLError,
                http.client.IncompleteRead,
                ConnectionResetError) as exc:
            if attempt == max_attempts:
                raise
            print(f"  [retry {attempt}/{max_attempts}] {exc} — waiting {delay}s")
            time.sleep(delay)
            delay *= 2


def accession_to_uid(accession: str) -> str | None:
    """Convert a BioProject accession (e.g. PRJNA123456) to its numeric UID."""
    handle = entrez_call(
        Entrez.esearch, db="bioproject", term=f"{accession}[PRJNA]", retmax=1
    )
    result = Entrez.read(handle)
    handle.close()
    ids = result.get("IdList", [])
    return ids[0] if ids else None


def get_sra_runs(bioproject_uid: str, bioproject_acc: str) -> list[dict]:
    """
    Return runinfo rows for all SRA runs linked to a BioProject UID.
    Each row is a dict of all runinfo CSV columns plus 'bioproject_accession'.
    Only rows with LibraryStrategy in {METAGENOMIC, WGS} are returned.
    """
    # Step 1: BioProject UID → SRA experiment UIDs via elink
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
        f"?dbfrom=bioproject&db=sra&id={bioproject_uid}&retmode=xml"
    )
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()

    root = ET.fromstring(resp.text)
    sra_uids = [
        el.text
        for el in root.findall(".//LinkSetDb[LinkName='bioproject_sra']/Link/Id")
    ]
    if not sra_uids:
        return []

    # Step 2: SRA UIDs → full runinfo rows via efetch
    all_rows: list[dict] = []
    batch_size = 200

    for i in range(0, len(sra_uids), batch_size):
        batch = sra_uids[i : i + batch_size]
        handle = entrez_call(
            Entrez.efetch,
            db="sra",
            id=",".join(batch),
            rettype="runinfo",
            retmode="text",
        )
        raw = handle.read()
        handle.close()
        csv_text = raw.decode("utf-8") if isinstance(raw, bytes) else raw

        # Parse runinfo CSV (may contain multiple blank-line-separated blocks)
        reader = csv.DictReader(io.StringIO(csv_text))
        for row in reader:
            if not row.get("Run"):          # skip blank / malformed rows
                continue
            if row.get("LibraryStrategy", "").upper() not in {"METAGENOMIC", "WGS"}:
                continue
            row["bioproject_accession"] = bioproject_acc
            all_rows.append(row)

    return all_rows


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------

def check_tool(name: str) -> bool:
    return subprocess.run(["which", name], capture_output=True).returncode == 0


def download_run(run: str, outdir: Path, threads: int, use_prefetch: bool) -> bool:
    """Download a single SRA run as FASTQ. Returns True on success."""
    # Skip if already downloaded
    existing = list(outdir.glob(f"{run}*.fastq*")) + list((outdir / run).glob("*.fastq*"))
    if existing:
        print(f"  [skip] {run} — already downloaded")
        return True

    if use_prefetch:
        print(f"  [prefetch] {run}")
        result = subprocess.run(
            ["prefetch", run, "--output-directory", str(outdir)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"  [ERROR] prefetch failed for {run}:\n{result.stderr.strip()}")
            return False

    print(f"  [fasterq-dump] {run}")
    fasterq_cmd = [
        "fasterq-dump", run,
        "--outdir", str(outdir),
        "--threads", str(threads),
        "--split-files",
        "--skip-technical",
    ]
    # Use local .sra file if prefetch placed it there
    if use_prefetch:
        sra_file = outdir / run / f"{run}.sra"
        if sra_file.exists():
            fasterq_cmd[1] = str(sra_file)

    result = subprocess.run(fasterq_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] fasterq-dump failed for {run}:\n{result.stderr.strip()}")
        return False

    return True


def save_metadata(runs: list[dict], path: Path) -> None:
    """Write run metadata to a CSV file."""
    if not runs:
        return
    # Put key columns first, then the rest alphabetically
    priority = [
        "Run", "bioproject_accession", "BioProject", "BioSample",
        "SampleName", "LibraryStrategy", "LibraryLayout",
        "Platform", "Model", "spots", "bases", "size_MB",
        "Organism", "SRAStudy", "Experiment",
    ]
    all_keys = list(runs[0].keys())
    ordered = priority + [k for k in all_keys if k not in priority]
    fieldnames = [k for k in ordered if k in all_keys]

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(runs)
    print(f"  Metadata saved → {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Download METAGENOMIC and WGS SRA runs for given BioProject accessions."
    )
    parser.add_argument(
        "--accessions", nargs="+", required=True, metavar="PRJNA",
        help="One or more BioProject accessions (e.g. PRJNA123456)",
    )
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument(
        "--outdir", default="sra_downloads",
        help="Output directory for FASTQ files (default: sra_downloads)",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Threads for fasterq-dump (default: 4)",
    )
    parser.add_argument(
        "--max-runs", type=int, default=None,
        help="Limit total number of runs downloaded (for testing)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="List runs that would be downloaded without downloading",
    )
    parser.add_argument(
        "--no-prefetch", action="store_true",
        help="Skip prefetch, use fasterq-dump directly",
    )
    args = parser.parse_args()

    Entrez.email = args.email
    Entrez.api_key = None

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Check tools
    if not args.dry_run:
        if not check_tool("fasterq-dump"):
            sys.exit("[ERROR] fasterq-dump not found. Install SRA Toolkit.")
        if not args.no_prefetch and not check_tool("prefetch"):
            print("[WARNING] prefetch not found — using fasterq-dump directly.")
            args.no_prefetch = True

    use_prefetch = not args.no_prefetch

    # -----------------------------------------------------------------------
    # Collect METAGENOMIC runs across all BioProjects
    # -----------------------------------------------------------------------
    all_run_meta: list[dict] = []

    for accession in args.accessions:
        accession = accession.strip().upper()
        print(f"\n[{accession}] Resolving BioProject UID...")
        uid = accession_to_uid(accession)
        if not uid:
            print(f"  [WARNING] Could not find UID for {accession} — skipping")
            continue
        print(f"  UID: {uid}")
        print(f"  Fetching SRA run list (METAGENOMIC + WGS)...")

        rows = get_sra_runs(uid, accession)
        print(f"  METAGENOMIC/WGS runs found: {len(rows)}")
        all_run_meta.extend(rows)

    if not all_run_meta:
        print("\nNo METAGENOMIC runs found. Exiting.")
        return

    print(f"\nTotal METAGENOMIC/WGS runs: {len(all_run_meta)}")

    if args.max_runs:
        all_run_meta = all_run_meta[: args.max_runs]
        print(f"  (limited to {args.max_runs} by --max-runs)")

    # -----------------------------------------------------------------------
    # Dry run: show table and save metadata, skip download
    # -----------------------------------------------------------------------
    if args.dry_run:
        print(f"\n{'Run':<15} {'BioProject':<14} {'LibraryLayout':<14} "
              f"{'spots':>10}  {'Organism'}")
        print("-" * 75)
        for r in all_run_meta:
            print(f"{r.get('Run',''):<15} {r.get('bioproject_accession',''):<14} "
                  f"{r.get('LibraryLayout',''):<14} "
                  f"{r.get('spots',''):>10}  {r.get('Organism','')}")
        meta_path = outdir / "sra_metadata.csv"
        save_metadata(all_run_meta, meta_path)
        return

    # -----------------------------------------------------------------------
    # Save metadata before downloading
    # -----------------------------------------------------------------------
    meta_path = outdir / "sra_metadata.csv"
    save_metadata(all_run_meta, meta_path)

    # -----------------------------------------------------------------------
    # Download
    # -----------------------------------------------------------------------
    success, failed = 0, []
    for i, row in enumerate(all_run_meta, 1):
        run = row["Run"]
        proj = row["bioproject_accession"]
        print(f"\n[{i}/{len(all_run_meta)}] {run}  ({proj})")
        ok = download_run(run, outdir, args.threads, use_prefetch)
        if ok:
            success += 1
        else:
            failed.append(run)

    print(f"\n--- Summary ---")
    print(f"  Downloaded : {success}/{len(all_run_meta)}")
    print(f"  Metadata   : {meta_path}")
    if failed:
        print(f"  Failed ({len(failed)}):")
        for r in failed:
            print(f"    {r}")
        failed_file = outdir / "failed_runs.txt"
        failed_file.write_text("\n".join(failed) + "\n")
        print(f"  Failed list: {failed_file}")


if __name__ == "__main__":
    main()
