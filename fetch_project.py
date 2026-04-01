#!/usr/bin/env python3
"""
fetch_project.py  —  Step 1: Fetch BioProject Metadata
═══════════════════════════════════════════════════════
Queries NCBI for all BioProjects linked to a taxon ID and saves a
metadata table (CSV + JSON).

Each project is automatically grouped into one of four sequence-type
categories based on keywords in the project title and description:

    16S       → 16S rRNA amplicon studies
    Amplicon  → amplicon sequencing (non-16S)
    RNA       → RNA-seq / metatranscriptomics
    Other     → WGS, metagenomics, or unclassified

Collected fields per project
─────────────────────────────
  bioproject_id          PRJNA / PRJEB / PRJDB accession
  title                  Project title
  description            Full description
  project_data_type      metagenome / raw sequence reads / targeted loci ...
  organism_name          Organism name from NCBI taxonomy
  taxon_id               NCBI Taxon ID of the target organism
  biosample_count        Number of BioSamples linked to the project
  sra_library_strategies Pipe-separated SRA library types (WGS, AMPLICON ...)
  seq_type               Auto-assigned group: 16S | Amplicon | RNA | Other
  registration_date      Date project was registered with NCBI
  submission_date        Date project was submitted
  pubmed_ids             Associated PubMed IDs (pipe-separated)
  relevance              NCBI domain tag (Agricultural, Medical ...)

Usage
─────
  python fetch_project.py --taxon-id 1510822 --email your@email.com
  python fetch_project.py --taxon-id 408170  --email your@email.com --label human_gut

Requirements
────────────
  pip install biopython pandas requests tqdm
"""

import argparse
import re
import http.client
import json
import os
import time
import urllib.error
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

import pandas as pd
import requests
from Bio import Entrez
from tqdm import tqdm

# ── constants ──────────────────────────────────────────────────────────────────
BATCH_SIZE  = 200   # IDs per efetch / esearch page
ELINK_BATCH = 50    # UIDs per elink call (SRA)
SLEEP       = 0.35  # seconds between requests  (~3 req/s, NCBI free tier)

ELINK_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"

# Keywords that define each seq_type group
SEQ_TYPE_RULES = {
    "16S":      r"16S",
    "Amplicon": r"[Aa]mplicon",
    "RNA":      r"\bRNA\b|\bRNA-[Ss]eq\b|\bmetatranscriptom",
}


# ══════════════════════════════════════════════════════════════════════════════
# NCBI helpers
# ══════════════════════════════════════════════════════════════════════════════

def entrez_call(fn, *args, retries: int = 5, **kwargs):
    """
    Call any Bio.Entrez function with automatic exponential-backoff retry.
    Handles connection drops and HTTP errors from NCBI.

    Backoff schedule: 5 → 10 → 20 → 40 → 80 seconds
    """
    delay = 5
    for attempt in range(retries):
        try:
            handle = fn(*args, **kwargs)
            result = Entrez.read(handle)
            handle.close()
            return result
        except (RuntimeError, urllib.error.HTTPError, urllib.error.URLError) as exc:
            if attempt == retries - 1:
                raise
            print(f"\n  [retry {attempt + 1}/{retries}] {str(exc)[:120]} "
                  f"— waiting {delay}s")
            time.sleep(delay)
            delay *= 2


def esearch_ids(db: str, term: str) -> list[str]:
    """
    Search an NCBI database and return all matching IDs.

    Uses usehistory to page through large result sets without
    repeating the search query on each page.
    """
    # First call: get total count + server-side session bookmark
    handle = Entrez.esearch(db=db, term=term, usehistory="y", retmax=0)
    record = Entrez.read(handle)
    handle.close()

    total     = int(record["Count"])
    webenv    = record["WebEnv"]
    query_key = record["QueryKey"]
    print(f"  Found {total:,} records in {db}")

    ids = []
    for start in range(0, total, BATCH_SIZE):
        handle = Entrez.esearch(
            db=db, term=term,
            webenv=webenv, query_key=query_key,
            retstart=start, retmax=BATCH_SIZE,
        )
        batch = Entrez.read(handle)
        handle.close()
        ids.extend(batch["IdList"])
        time.sleep(SLEEP)

    return ids


# ══════════════════════════════════════════════════════════════════════════════
# Step 2 — BioProject XML fetch + parse
# ══════════════════════════════════════════════════════════════════════════════

def fetch_bioproject_records(bioproject_ids: list[str],
                             checkpoint_file: str | None = None) -> list[dict]:
    """
    Download full BioProject XML in batches of BATCH_SIZE and parse
    each record into a flat dict.

    Checkpoint / resume
    ───────────────────
    If checkpoint_file is given, progress is saved after every batch.
    On restart the already-fetched records are loaded and the download
    continues from the next unfinished batch — no work is repeated.
    The checkpoint file is deleted automatically on successful completion.
    """
    # Load existing checkpoint if available
    records: list[dict] = []
    start_batch = 0
    if checkpoint_file and os.path.exists(checkpoint_file):
        with open(checkpoint_file) as fh:
            ckpt = json.load(fh)
        records     = ckpt["records"]
        start_batch = ckpt["next_batch"]
        print(f"  Resuming from batch {start_batch} "
              f"({len(records)} records already fetched)")

    batches = list(range(0, len(bioproject_ids), BATCH_SIZE))

    for batch_idx in tqdm(range(start_batch, len(batches)),
                          desc="Fetching BioProject records",
                          initial=start_batch, total=len(batches)):
        start = batches[batch_idx]
        batch = bioproject_ids[start: start + BATCH_SIZE]
        delay = 5
        for attempt in range(5):
            try:
                handle   = Entrez.efetch(db="bioproject",
                                         id=",".join(batch), rettype="xml")
                xml_text = handle.read()
                handle.close()
                break
            except (RuntimeError,
                    urllib.error.HTTPError,
                    urllib.error.URLError,
                    http.client.IncompleteRead,
                    ConnectionResetError) as exc:
                if attempt == 4:
                    raise
                print(f"\n  [retry {attempt + 1}/5] {str(exc)[:120]} "
                      f"— waiting {delay}s")
                time.sleep(delay)
                delay *= 2

        time.sleep(SLEEP)
        records.extend(_parse_bioproject_xml(xml_text))

        # Save checkpoint after every successful batch
        if checkpoint_file:
            with open(checkpoint_file, "w") as fh:
                json.dump({"records": records,
                           "next_batch": batch_idx + 1}, fh)

    # Clean up checkpoint on success
    if checkpoint_file and os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)

    return records


def _parse_bioproject_xml(xml_text: str | bytes) -> list[dict]:
    """
    Parse raw BioProject XML.

    NCBI efetch returns one XML document per record. All declarations
    are stripped before wrapping into a single root so Python's
    ElementTree can parse the whole batch at once.
    """
    if isinstance(xml_text, bytes):
        xml_text = xml_text.decode("utf-8", errors="replace")
    xml_text = re.sub(r"<\?xml[^?]*\?>", "", xml_text)
    xml_text = re.sub(r"<!DOCTYPE[^>]*>",  "", xml_text)
    xml_text = f"<RecordSet>{xml_text.strip()}</RecordSet>"

    root    = ET.fromstring(xml_text)
    results = []

    for pkg in root.iter("DocumentSummary"):
        rec = {}

        # IDs
        rec["bioproject_id"]  = _attr(pkg, "Project/ProjectID/ArchiveID", "accession")
        rec["bioproject_uid"] = pkg.get("uid", "")

        # Title / description
        rec["title"]       = _text(pkg, "Project/ProjectDescr/Title")
        rec["description"] = _text(pkg, "Project/ProjectDescr/Description")

        # Project type
        rec["project_data_type"] = _text(
            pkg, "Project/ProjectType/ProjectTypeSubmission"
                 "/ProjectDataTypeSet/DataType")

        # Organism
        rec["organism_name"] = _text(
            pkg, "Project/ProjectType/ProjectTypeSubmission"
                 "/Target/Organism/OrganismName")
        rec["taxon_id"] = _text(
            pkg, "Project/ProjectType/ProjectTypeSubmission"
                 "/Target/Organism/TaxID")

        # Dates
        rec["registration_date"] = (
            _text(pkg, "Project/ProjectDescr/ProjectReleaseDate")
            or _text(pkg, "Submission/Description/Access")
        )
        rec["submission_date"] = _attr(pkg, "Submission", "submitted")

        # Publications
        pmids = [el.text for el in
                 pkg.findall(".//Publication/DbType[.='ePubmed']/../ID")]
        rec["pubmed_ids"] = "|".join(pmids) if pmids else ""

        # Relevance domain tag
        rel = pkg.find("Project/ProjectDescr/Relevance")
        rec["relevance"] = (
            "|".join(child.tag for child in rel)
            if rel is not None else ""
        )

        results.append(rec)

    return results


def _text(element, xpath: str) -> str:
    node = element.find(xpath)
    return (node.text or "").strip() if node is not None else ""


def _attr(element, xpath: str, attr: str) -> str:
    node = element.find(xpath)
    return node.get(attr, "") if node is not None else ""


# ══════════════════════════════════════════════════════════════════════════════
# Step 3 — SRA library strategies
# ══════════════════════════════════════════════════════════════════════════════

def fetch_sra_strategies(uids: list[str]) -> dict[str, set[str]]:
    """
    For each BioProject UID, collect the set of SRA library strategies
    (AMPLICON, WGS, METAGENOMICS, RNA-Seq …) by linking to SRA and
    reading the ExpXml field from esummary.
    """
    uid_to_strategies: dict[str, set[str]] = defaultdict(set)

    print(f"\nLinking {len(uids):,} BioProjects → SRA ...")
    for start in tqdm(range(0, len(uids), ELINK_BATCH),
                      desc="BioProject→SRA"):
        batch = uids[start: start + ELINK_BATCH]
        try:
            link_records = entrez_call(Entrez.elink,
                                       dbfrom="bioproject", db="sra",
                                       id=",".join(batch))
        except Exception:
            continue
        time.sleep(SLEEP)

        for link_set in link_records:
            bp_uid  = link_set.get("IdList", [""])[0]
            sra_ids = [lnk["Id"] for lnk in
                       link_set.get("LinkSetDb", [{}])[0].get("Link", [])
                       if link_set.get("LinkSetDb")]
            if not sra_ids:
                continue

            try:
                summaries = entrez_call(Entrez.esummary,
                                        db="sra",
                                        id=",".join(sra_ids[:200]))
                time.sleep(SLEEP)
            except Exception:
                continue

            for s in summaries:
                exp_xml = s.get("ExpXml", "")
                if not exp_xml:
                    continue
                try:
                    eroot = ET.fromstring(f"<r>{exp_xml}</r>")
                    for ls in eroot.iter("LIBRARY_STRATEGY"):
                        if ls.text:
                            uid_to_strategies[bp_uid].add(
                                ls.text.strip().upper())
                except ET.ParseError:
                    pass

    return uid_to_strategies


# ══════════════════════════════════════════════════════════════════════════════
# Step 4 — BioSample counts
# ══════════════════════════════════════════════════════════════════════════════

def fetch_biosample_counts(uids: list[str]) -> dict[str, int]:
    """
    For each BioProject UID, count linked BioSamples via the elink REST API.

    Each UID is queried individually — batching merges all BioSample IDs
    into one undifferentiated pool with no per-project breakdown.
    Uses linkname=bioproject_biosample_all to include sub-project samples.
    """
    uid_to_count: dict[str, int] = {}

    print(f"\nFetching BioSample counts for {len(uids):,} projects ...")
    for uid in tqdm(uids, desc="BioProject→BioSample"):
        params = {
            "dbfrom":   "bioproject",
            "db":       "biosample",
            "id":       uid,
            "linkname": "bioproject_biosample_all",
            "retmode":  "xml",
            "email":    Entrez.email,
        }
        delay = 5
        for attempt in range(5):
            try:
                resp = requests.get(ELINK_URL, params=params, timeout=30)
                resp.raise_for_status()
                xml  = re.sub(r"<\?xml[^?]*\?>",  "", resp.text)
                xml  = re.sub(r"<!DOCTYPE[^>]*>", "", xml)
                root = ET.fromstring(xml.strip())
                uid_to_count[uid] = len(root.findall(".//LinkSetDb/Link/Id"))
                break
            except (requests.exceptions.SSLError,
                    requests.exceptions.ConnectionError,
                    requests.exceptions.Timeout,
                    requests.exceptions.HTTPError) as exc:
                if attempt == 4:
                    uid_to_count[uid] = 0
                    break
                print(f"\n  [retry {attempt + 1}/5] {str(exc)[:120]} "
                      f"— waiting {delay}s")
                time.sleep(delay)
                delay *= 2
        time.sleep(SLEEP)

    return uid_to_count


# ══════════════════════════════════════════════════════════════════════════════
# Seq-type grouping
# ══════════════════════════════════════════════════════════════════════════════

def assign_seq_type(title: str, description: str) -> str:
    """
    Assign a sequence-type label from project title + description.

    Rules (applied in order, multiple labels are pipe-joined):
        16S       → contains "16S"
        Amplicon  → contains "Amplicon" / "amplicon"
        RNA       → contains "RNA", "RNA-Seq", or "metatranscriptom"
        Other     → none of the above matched
    """
    text    = f"{title} {description}"
    matched = [k for k, pat in SEQ_TYPE_RULES.items()
               if re.search(pat, text)]
    return "|".join(matched) if matched else "Other"


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        prog="fetch_project.py",
        description="Step 1 — Fetch BioProject metadata from NCBI for a given taxon.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python fetch_project.py --taxon-id 1510822 --email your@email.com
  python fetch_project.py --taxon-id 408170  --email your@email.com --label human_gut
        """,
    )
    parser.add_argument("--taxon-id", required=True,
                        help="NCBI Taxon ID  (pig gut: 1510822 | human gut: 408170)")
    parser.add_argument("--email",    required=True,
                        help="Your email address — required by NCBI Entrez policy")
    parser.add_argument("--label",    default=None,
                        help="Output filename prefix  (default: projects_taxid<ID>)")
    args = parser.parse_args()

    Entrez.email   = args.email
    Entrez.api_key = None
    taxon_id       = args.taxon_id
    label          = args.label or f"projects_taxid{taxon_id}"

    print("=" * 62)
    print(f"  fetch_project.py  —  taxon {taxon_id}")
    print("=" * 62)

    # ── Step 1: search ─────────────────────────────────────────────
    print("\n[Step 1/4] Searching BioProject ...")
    bp_ids = esearch_ids("bioproject", f"txid{taxon_id}[ORGN]")

    # ── Step 2: fetch full records ──────────────────────────────────
    print("\n[Step 2/4] Fetching BioProject records ...")
    checkpoint_file = f"{label}_checkpoint.json"
    records = fetch_bioproject_records(bp_ids, checkpoint_file=checkpoint_file)

    uid_to_acc = {r["bioproject_uid"]: r["bioproject_id"] for r in records}
    uids       = [r["bioproject_uid"] for r in records if r["bioproject_uid"]]

    # ── Step 3: SRA library strategies ─────────────────────────────
    print("\n[Step 3/4] Fetching SRA library strategies ...")
    uid_to_strategies = fetch_sra_strategies(uids)
    acc_to_strategies = {
        uid_to_acc.get(uid, uid): "|".join(sorted(s))
        for uid, s in uid_to_strategies.items()
    }

    # ── Step 4: BioSample counts ────────────────────────────────────
    print("\n[Step 4/4] Fetching BioSample counts ...")
    uid_to_count = fetch_biosample_counts(uids)

    # ── Merge ───────────────────────────────────────────────────────
    for rec in records:
        rec["sra_library_strategies"] = acc_to_strategies.get(
            rec["bioproject_id"],
            acc_to_strategies.get(rec["bioproject_uid"], ""),
        )
        rec["biosample_count"] = uid_to_count.get(rec["bioproject_uid"], 0)

    df = pd.DataFrame(records)

    # Auto-assign seq_type group
    df["seq_type"] = df.apply(
        lambda r: assign_seq_type(r["title"], r["description"]), axis=1
    )

    # Flag projects whose biosample_count exceeds the mean of the "Other" group.
    # The average is derived from "Other" only — 16S / Amplicon / RNA projects
    # are excluded from the baseline because they are not the target.
    df["biosample_count"] = pd.to_numeric(df["biosample_count"], errors="coerce").fillna(0)
    other_mean = df.loc[df["seq_type"] == "Other", "biosample_count"].mean()
    df["above_avg_biosample"] = df["biosample_count"].apply(
        lambda x: "above_avg" if x > other_mean else "below_avg"
    )

    # Column order
    col_order = [
        "bioproject_id", "bioproject_uid",
        "title", "description",
        "project_data_type",
        "organism_name", "taxon_id",
        "biosample_count",
        "above_avg_biosample",
        "sra_library_strategies",
        "seq_type",
        "registration_date", "submission_date",
        "pubmed_ids", "relevance",
    ]
    df = df[[c for c in col_order if c in df.columns]]

    # ── Save ────────────────────────────────────────────────────────
    out_csv  = f"{label}.csv"
    out_json = f"{label}.json"
    df.to_csv(out_csv,  index=False)
    df.to_json(out_json, orient="records", indent=2)

    # ── Summary ─────────────────────────────────────────────────────
    print(f"\n{'=' * 62}")
    print(f"  Done — {len(df):,} BioProjects saved")
    print(f"  {out_csv}")
    print(f"  {out_json}")
    print(f"{'=' * 62}")

    print(f"\nseq_type breakdown:")
    for tag, n in Counter(df["seq_type"]).most_common():
        bar = "█" * (n * 40 // len(df))
        print(f"  {tag:<12} {n:>5}  {bar}")

    print(f"\nBioSample counts:  "
          f"total={df['biosample_count'].sum():,}  "
          f"median={df['biosample_count'].median():.0f}  "
          f"max={df['biosample_count'].max()}")

    above_n = (df["above_avg_biosample"] == "above_avg").sum()
    below_n = (df["above_avg_biosample"] == "below_avg").sum()
    print(f"\nabove_avg_biosample threshold (mean of 'Other' group): "
          f"{other_mean:.1f}")
    print(f"  above_avg: {above_n} / {len(df)} "
          f"({100 * above_n / len(df):.1f}%)")
    print(f"  below_avg: {below_n} / {len(df)} "
          f"({100 * below_n / len(df):.1f}%)")

    print(f"\nTop SRA library strategies:")
    strats = [s for row in df["sra_library_strategies"].dropna()
              for s in row.split("|") if s]
    for strat, n in Counter(strats).most_common(10):
        print(f"  {strat:<30} {n}")


if __name__ == "__main__":
    main()
