# Gut Microbiome — NCBI Data Download Pipeline

A two-step pipeline to collect BioProject metadata and download SRA sequencing
data from NCBI for any gut microbiome taxon.

---

## Table of Contents

1. [Pipeline Overview](#1-pipeline-overview)
2. [Prerequisites](#2-prerequisites)
3. [Quick Start](#3-quick-start)
4. [Step 1 — Fetch Project Metadata](#4-step-1--fetch-project-metadata-fetch_projectpy)
5. [Step 2 — Download SRA FASTQ](#5-step-2--download-sra-fastq-download_srapy)
6. [Output Files](#6-output-files)
7. [seq_type Groups](#7-seq_type-groups)
8. [SRA Library Strategy Reference](#8-sra-library-strategy-reference)
9. [Taxon ID Reference](#9-taxon-id-reference)
10. [Troubleshooting](#10-troubleshooting)

---

## 1. Pipeline Overview

```
┌──────────────────────────────────────────────────────────────┐
│                     INPUT: NCBI Taxon ID                      │
│          pig gut = 1510822  |  human gut = 408170             │
└────────────────────────────┬─────────────────────────────────┘
                             │
                             ▼
┌──────────────────────────────────────────────────────────────┐
│                  STEP 1 — fetch_project.py                    │
│                                                               │
│  esearch (BioProject)                                         │
│    → efetch XML (title, description, dates, pubmed)           │
│    → elink BioProject→SRA (library strategies)                │
│    → elink BioProject→BioSample (sample counts)               │
│    → keyword grouping (seq_type)                              │
│    → above/below average BioSample labelling                  │
│                                                               │
│  Output:  <label>.csv   /   <label>.json                      │
└────────────────────────────┬─────────────────────────────────┘
                             │
             ┌───────────────┼──────────────┐
             ▼               ▼              ▼
          16S            Amplicon       Other / WGS
        projects         projects       projects  ◄── target
             │               │              │
             └───────────────┴──────────────┘
                             │
                  (pick BioProject accessions)
                             │
                             ▼
┌──────────────────────────────────────────────────────────────┐
│                  STEP 2 — download_sra.py                     │
│                                                               │
│  BioProject accession → NCBI UID                              │
│    → elink BioProject→SRA (experiment UIDs)                   │
│    → efetch runinfo CSV (filter METAGENOMIC + WGS)            │
│    → save sra_metadata.csv                                    │
│    → prefetch + fasterq-dump → FASTQ files                    │
│                                                               │
│  Output:  <outdir>/SRR*.fastq                                 │
│           <outdir>/sra_metadata.csv                           │
└──────────────────────────────────────────────────────────────┘
```

---

## 2. Prerequisites

### Python packages

```bash
pip install biopython pandas requests tqdm
```

### SRA Toolkit — required for Step 2 only

```bash
# Option A: conda (recommended)
conda install -c bioconda -c conda-forge sra-tools

# Option B: Homebrew (macOS)
brew install sratoolkit

# Option C: direct download from NCBI
#   https://github.com/ncbi/sra-tools/releases
#   Download the binary for your platform, extract, and add bin/ to PATH

# Verify installation
fasterq-dump --version
prefetch --version
```

> **Important:** Use SRA Toolkit **v3.x or newer**. Older versions (e.g. v2.9.x)
> have SSL certificate issues that prevent `prefetch` from connecting to NCBI.
> If upgrading is not possible, use the `--no-prefetch` flag (see Step 2).

### Initial SRA Toolkit configuration

Run once after installation:

```bash
vdb-config -i
# Accept defaults and save — this creates ~/.ncbi/user-settings.mkfg
```

### NCBI account (optional but recommended)

A free API key from https://www.ncbi.nlm.nih.gov/account/ raises the request
rate from 3 to 10 per second, cutting Step 1 runtime by ~3×.

---

## 3. Quick Start

```bash
# Step 1 — collect metadata for pig gut microbiome
python fetch_project.py \
    --taxon-id 1510822 \
    --email    your@email.com

# Step 2 — preview METAGENOMIC + WGS runs (no download)
python download_sra.py \
    --accessions PRJNA857725 PRJNA123456 \
    --email      your@email.com \
    --dry-run

# Step 2 — download
python download_sra.py \
    --accessions PRJNA857725 \
    --email      your@email.com \
    --outdir     ./fastq \
    --threads    8
```

---

## 4. Step 1 — Fetch Project Metadata (`fetch_project.py`)

### What it does

`fetch_project.py` queries four NCBI databases and builds a metadata table for
every BioProject linked to a taxon of interest.

```
fetch_project.py
│
├── [1/4]  esearch("bioproject", "txid<ID>[ORGN]")
│           └─ returns all BioProject UIDs for the taxon
│
├── [2/4]  efetch("bioproject", ids, rettype="xml")
│           └─ downloads full XML in batches of 200
│           └─ parses: accession, title, description, data_type,
│                      organism, dates, pubmed_ids, relevance
│           └─ checkpoint/resume: saves progress after every batch
│                      so a crash can be continued from where it stopped
│
├── [3/4]  elink(bioproject → sra)  +  esummary("sra")
│           └─ finds all SRA experiments per project
│           └─ reads LIBRARY_STRATEGY from ExpXml
│              (e.g. AMPLICON, WGS, METAGENOMIC, RNA-Seq)
│
└── [4/4]  elink REST (bioproject → biosample, one UID at a time)
            └─ counts linked BioSamples per project
            └─ uses linkname=bioproject_biosample_all
               (includes samples from sub-projects)
```

> **Why one UID at a time for BioSample?** Batching multiple UIDs in one
> elink call causes NCBI to merge all BioSample IDs into a single list with
> no per-project breakdown. Querying individually gives accurate counts.

### Usage

```bash
# Pig gut microbiome
python fetch_project.py \
    --taxon-id 1510822 \
    --email    your@email.com

# Human gut microbiome (large dataset — uses checkpoint/resume)
python fetch_project.py \
    --taxon-id 408170 \
    --email    your@email.com \
    --label    human_gut
```

### Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--taxon-id` | Yes | — | NCBI Taxon ID |
| `--email` | Yes | — | Your email (required by NCBI policy) |
| `--label` | No | `projects_taxid<ID>` | Output filename prefix |

### Checkpoint / Resume

For large datasets (e.g. human gut ~5,900 projects), the script saves a
checkpoint file `<label>_checkpoint.json` after every batch of 200 records.
If the script crashes or is interrupted, simply re-run the same command —
it will resume from the last saved batch automatically. The checkpoint file
is deleted on successful completion.

### Expected terminal output

```
══════════════════════════════════════════════════════════════
  fetch_project.py  —  taxon 1510822
══════════════════════════════════════════════════════════════

[Step 1/4] Searching BioProject ...
  Found 1,133 records in bioproject

[Step 2/4] Fetching BioProject records ...
Fetching BioProject records: 100%|██████| 6/6 [02:40<00:00]

[Step 3/4] Fetching SRA library strategies ...
BioProject→SRA: 100%|███████████████████| 23/23 [04:15<00:00]

[Step 4/4] Fetching BioSample counts ...
BioProject→BioSample: 100%|████████████| 1133/1133 [06:35<00:00]

══════════════════════════════════════════════════════════════
  Done — 1,133 BioProjects saved
  projects_taxid1510822.csv
  projects_taxid1510822.json
══════════════════════════════════════════════════════════════

seq_type breakdown:
  Other          842  ████████████████████████████████
  16S            187  ███████
  Amplicon        72  ██
  RNA             32  █

BioSample counts:  total=128,451  median=24  max=9,079

above_avg_biosample threshold (mean of 'Other' group): 87.3
  above_avg:  312 / 1133 (27.5%)
  below_avg:  821 / 1133 (72.5%)

Top SRA library strategies:
  AMPLICON                         412
  WGS                              289
  RNA-SEQ                           87
```

### Output columns

| Column | Example | Description |
|--------|---------|-------------|
| `bioproject_id` | PRJNA1248931 | NCBI accession |
| `bioproject_uid` | 1248931 | Internal NCBI numeric UID |
| `title` | Colonic Microbiota… | Project title |
| `description` | This study… | Full description |
| `project_data_type` | metagenome | NCBI data type |
| `organism_name` | pig gut metagenome | Target organism |
| `taxon_id` | 1510822 | NCBI Taxon ID |
| `biosample_count` | 119 | Number of linked BioSamples |
| `above_avg_biosample` | above_avg | Whether biosample_count exceeds mean of the "Other" group |
| `sra_library_strategies` | AMPLICON\|WGS | Pipe-separated SRA library types |
| `seq_type` | 16S\|Amplicon | Auto-assigned group (see §7) |
| `registration_date` | 2024-01-15 | NCBI registration date |
| `submission_date` | 2024-01-10 | Submission date |
| `pubmed_ids` | 38123456 | Linked PubMed IDs (pipe-separated) |
| `relevance` | Agricultural | NCBI domain tag |

### `above_avg_biosample` explained

The threshold is the **mean biosample_count of the "Other" seq_type group only**
(shotgun metagenomics / WGS projects). 16S, Amplicon, and RNA projects are
excluded from the baseline because they are not the primary download target.
All projects — including 16S/Amplicon/RNA — are then labelled against this
threshold:

- `above_avg` — biosample_count > mean of Other group
- `below_avg` — biosample_count ≤ mean of Other group

### Expected runtime

| Taxon | Projects | Approx. time (no API key) |
|-------|----------|--------------------------|
| Pig gut (1510822) | ~1,100 | 15–20 min |
| Human gut (408170) | ~5,900 | 90–120 min |

---

## 5. Step 2 — Download SRA FASTQ (`download_sra.py`)

### What it does

`download_sra.py` accepts one or more BioProject accessions, finds all
METAGENOMIC and WGS SRA runs linked to them, saves run metadata to a CSV,
then downloads FASTQ files using `prefetch` + `fasterq-dump`.

```
download_sra.py
│
├── For each BioProject accession:
│   │
│   ├── esearch("bioproject", "PRJNA…[PRJNA]")
│   │    └─ resolves accession → numeric UID
│   │
│   ├── elink REST (bioproject → sra, retmode=xml)
│   │    └─ returns SRA experiment UIDs
│   │
│   └── efetch("sra", ids, rettype="runinfo", retmode="text")
│        └─ downloads runinfo CSV (200 UIDs per batch)
│        └─ filters: LibraryStrategy in {METAGENOMIC, WGS}
│        └─ adds bioproject_accession column
│
├── Save sra_metadata.csv  (all runs, all runinfo columns)
│
└── For each SRR/ERR/DRR run:
     ├── skip if FASTQ already exists
     ├── prefetch <run> --output-directory <outdir>
     └── fasterq-dump <run> --split-files --skip-technical
```

### Always dry-run first

```bash
python download_sra.py \
    --accessions PRJNA857725 \
    --email      your@email.com \
    --dry-run
```

Dry-run output example:
```
[PRJNA857725] Resolving BioProject UID...
  UID: 857725
  Fetching SRA run list (METAGENOMIC + WGS)...
  METAGENOMIC/WGS runs found: 121

Total METAGENOMIC/WGS runs: 121

Run             BioProject      LibraryLayout       spots  Organism
---------------------------------------------------------------------------
SRR21388550     PRJNA857725     PAIRED             8523112  pig gut metagenome
SRR21388551     PRJNA857725     PAIRED             7841034  pig gut metagenome
...

  Metadata saved → sra_downloads/sra_metadata.csv
```

### Download

```bash
python download_sra.py \
    --accessions PRJNA857725 PRJNA123456 \
    --email      your@email.com \
    --outdir     ./fastq \
    --threads    8
```

### Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--accessions` | Yes | — | One or more BioProject accessions |
| `--email` | Yes | — | Your email (required by NCBI policy) |
| `--outdir` | No | `sra_downloads` | Output directory for FASTQ files |
| `--threads` | No | `4` | Threads passed to `fasterq-dump` |
| `--max-runs` | No | None | Cap total runs (useful for testing) |
| `--dry-run` | No | False | List runs and save metadata without downloading |
| `--no-prefetch` | No | False | Skip `prefetch`, use `fasterq-dump` directly |

### Resume interrupted downloads

`download_sra.py` automatically skips runs whose FASTQ files already exist in
`--outdir`. To resume after a crash, simply re-run the same command.

Failed runs are saved to `<outdir>/failed_runs.txt`. To retry them:

```bash
while read srr; do
    fasterq-dump "$srr" --outdir ./fastq --split-files --threads 4
done < fastq/failed_runs.txt
```

### Run metadata CSV columns

The `sra_metadata.csv` produced before downloading contains all runinfo
columns from NCBI, including:

| Column | Description |
|--------|-------------|
| `Run` | SRR/ERR/DRR accession |
| `bioproject_accession` | Input BioProject accession |
| `BioProject` | NCBI BioProject ID |
| `BioSample` | Linked BioSample accession |
| `SampleName` | Sample name |
| `LibraryStrategy` | `METAGENOMIC` or `WGS` |
| `LibraryLayout` | `PAIRED` or `SINGLE` |
| `Platform` | `ILLUMINA`, `OXFORD_NANOPORE`, etc. |
| `Model` | Instrument model |
| `spots` | Total read pairs / spots |
| `bases` | Total bases sequenced |
| `size_MB` | Approximate file size |
| `Organism` | Organism name |
| `SRAStudy` | SRA study accession (SRP…) |
| `Experiment` | SRA experiment accession (SRX…) |

---

## 6. Output Files

```
pig_dataset/
│
├── fetch_project.py               ← Step 1 script
├── download_sra.py                ← Step 2 script
├── README.md                      ← this document
│
├── projects_taxid1510822.csv      ← Step 1 output: BioProject metadata
├── projects_taxid1510822.json
│
└── fastq/                         ← Step 2 output
    ├── sra_metadata.csv           ← run-level metadata for all downloaded runs
    ├── SRR21388550_1.fastq        ← forward reads (PAIRED)
    ├── SRR21388550_2.fastq        ← reverse reads (PAIRED)
    ├── SRR21388551.fastq          ← single-end reads
    └── failed_runs.txt            ← SRR accessions that failed (if any)
```

---

## 7. seq_type Groups

`fetch_project.py` automatically assigns each project to a group using keyword
matching on the project title and description.

```
Title + Description text
        │
        ├── contains "16S"                      ──► seq_type: "16S"
        │
        ├── contains "Amplicon" / "amplicon"     ──► seq_type: "Amplicon"
        │
        ├── contains "RNA", "RNA-Seq",           ──► seq_type: "RNA"
        │   or "metatranscriptom"
        │
        └── none of the above                   ──► seq_type: "Other"
```

A project can carry **multiple labels** (pipe-joined) if its text matches more
than one rule — e.g. a 16S amplicon study gets `16S|Amplicon`.

| seq_type | Typical studies |
|----------|----------------|
| `16S` | 16S rRNA gene amplicon sequencing |
| `Amplicon` | Amplicon sequencing (ITS, 18S, etc.) |
| `RNA` | RNA-seq, metatranscriptomics |
| `Other` | Shotgun metagenomics, WGS, unclassified |

---

## 8. SRA Library Strategy Reference

The `LibraryStrategy` field in the SRA runinfo CSV uses NCBI controlled vocabulary.
`download_sra.py` filters for **METAGENOMIC** and **WGS**.

| Strategy | Meaning | Relevant? |
|----------|---------|-----------|
| **METAGENOMIC** | Shotgun metagenomics | **Yes — downloaded** |
| **WGS** | Whole Genome Shotgun | **Yes — downloaded** |
| AMPLICON | Targeted amplicon (e.g. 16S rRNA) | No (use separate pipeline) |
| RNA-Seq | Transcriptome sequencing | No |
| miRNA-Seq | Small/micro RNA | No |
| ncRNA-Seq | Non-coding RNA | No |
| ChIP-Seq | Chromatin immunoprecipitation | No |
| Bisulfite-Seq | DNA methylation | No |
| WXS | Whole Exome Sequencing | No |
| OTHER | Does not fit other categories | No |

To include additional strategies, edit line in `get_sra_runs()`:

```python
# Current
if row.get("LibraryStrategy", "").upper() not in {"METAGENOMIC", "WGS"}:

# To also include AMPLICON
if row.get("LibraryStrategy", "").upper() not in {"METAGENOMIC", "WGS", "AMPLICON"}:
```

---

## 9. Taxon ID Reference

| Organism | Taxon ID | NCBI link |
|----------|----------|-----------|
| Pig gut metagenome | 1510822 | https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1510822 |
| Human gut metagenome | 408170 | https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=408170 |
| Mouse gut metagenome | 1268534 | https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1268534 |
| Chicken gut metagenome | 1385655 | https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1385655 |

Find any taxon: https://www.ncbi.nlm.nih.gov/taxonomy

---

## 10. Troubleshooting

### Common errors

| Error | Cause | Fix |
|-------|-------|-----|
| `HTTP Error 400` | `api_key=""` (empty string) | Set to `None` or a valid key |
| `RuntimeError: EOF` | NCBI dropped connection | Script retries automatically (5×). Wait and re-run if persistent. |
| `biosample_count = 0` | Wrong elink linkname | Must use `bioproject_biosample_all`, not `bioproject_biosample` |
| `fasterq-dump not found` | SRA Toolkit not installed | `conda install -c bioconda sra-tools` |
| `IncompleteRead` error | Connection reset mid-transfer | Script retries automatically. For large runs, re-run to resume. |
| `PackagesNotFoundError: sra-tools=3.1.1` | Exact version unavailable for your platform | Use `conda install -c bioconda sra-tools` without a version pin |

### prefetch SSL certificate error

```
mbedtls_ssl_handshake returned -9984
Certificate verification failed
The certificate is signed with an unacceptable hash
```

**Cause:** SRA Toolkit v2.9.x uses an outdated TLS library that rejects modern
NCBI certificates.

**Fixes (in order of preference):**

```bash
# Option 1: upgrade SRA Toolkit
conda install -c bioconda -c conda-forge sra-tools   # no version pin

# Option 2: Homebrew (macOS)
brew install sratoolkit

# Option 3: bypass prefetch entirely (immediate workaround)
python download_sra.py \
    --accessions PRJNA857725 \
    --email      your@email.com \
    --no-prefetch
```

### Resume a large Step 1 run

`fetch_project.py` saves a checkpoint file `<label>_checkpoint.json` after each
batch. If the script is interrupted, re-run the exact same command:

```bash
python fetch_project.py --taxon-id 408170 --email your@email.com --label human_gut
# → "Resuming from batch 12 (2400 records already fetched)"
```

### Resume failed Step 2 downloads

```bash
# Retry all SRR accessions listed in failed_runs.txt
while read srr; do
    fasterq-dump "$srr" --outdir ./fastq --split-files --threads 4
done < fastq/failed_runs.txt
```

### Test with a small number of runs

```bash
# Download only the first 3 runs to verify the pipeline
python download_sra.py \
    --accessions PRJNA857725 \
    --email      your@email.com \
    --max-runs   3 \
    --outdir     ./test_fastq
```
