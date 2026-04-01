# Gut Microbiome — NCBI Data Download Pipeline

Two scripts to collect BioProject metadata and download SRA sequencing data from NCBI.

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Quick Start](#2-quick-start)
3. [Step 1 — Fetch Project Metadata](#3-step-1--fetch-project-metadata)
4. [Step 2 — Download SRA FASTQ](#4-step-2--download-sra-fastq)
5. [Output Files](#5-output-files)
6. [seq_type Groups](#6-seq_type-groups)
7. [Taxon ID Reference](#7-taxon-id-reference)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Prerequisites

### Python packages

```bash
pip install biopython pandas requests tqdm
```

### SRA Toolkit (Step 2 only)

```bash
# conda
conda install -c bioconda -c conda-forge sra-tools

# macOS Homebrew
brew install sratoolkit
```

> Use SRA Toolkit **v3.x or newer**. Older versions have SSL issues with `prefetch`.
> If you cannot upgrade, use the `--no-prefetch` flag.

Run once after installation:

```bash
vdb-config -i    # accept defaults and save
```

---

## 2. Quick Start

```bash
# Step 1 — collect metadata
python fetch_project.py \
    --taxon-id 1510822 \
    --email    your@email.com

# Step 2 — preview runs (no download)
python download_sra.py \
    --accessions PRJNA857725 \
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

## 3. Step 1 — Fetch Project Metadata

**Script:** `fetch_project.py`

Searches NCBI BioProject by taxon ID and builds a metadata table (CSV + JSON).
Each project is automatically grouped by sequence type and labelled by BioSample count.

### Usage

```bash
# Pig gut
python fetch_project.py --taxon-id 1510822 --email your@email.com

# Human gut (large dataset — checkpoint/resume enabled)
python fetch_project.py --taxon-id 408170 --email your@email.com --label human_gut
```

### Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--taxon-id` | Yes | — | NCBI Taxon ID |
| `--email` | Yes | — | Your email (NCBI policy) |
| `--label` | No | `projects_taxid<ID>` | Output filename prefix |

### Output columns

| Column | Description |
|--------|-------------|
| `bioproject_id` | NCBI accession (e.g. PRJNA123456) |
| `title` | Project title |
| `description` | Full description |
| `project_data_type` | metagenome / raw sequence reads / etc. |
| `organism_name` | Target organism |
| `taxon_id` | NCBI Taxon ID |
| `biosample_count` | Number of linked BioSamples |
| `above_avg_biosample` | `above_avg` or `below_avg` (threshold = mean of "Other" group) |
| `sra_library_strategies` | Pipe-separated SRA library types (e.g. `AMPLICON\|WGS`) |
| `seq_type` | Auto-assigned group: `16S`, `Amplicon`, `RNA`, or `Other` |
| `registration_date` | Date registered with NCBI |
| `pubmed_ids` | Linked PubMed IDs |
| `relevance` | NCBI domain tag (e.g. Agricultural, Medical) |

### Notes

- **Checkpoint/resume:** progress is saved to `<label>_checkpoint.json` after every
  batch. If the script crashes, re-run the same command to continue from where it stopped.
- **`above_avg_biosample`:** threshold is the mean biosample_count of the `Other`
  group only (WGS / metagenomics). All projects are then labelled `above_avg` or
  `below_avg` against this threshold.

### Expected runtime

| Taxon | Projects | Approx. time |
|-------|----------|-------------|
| Pig gut (1510822) | ~1,100 | 15–20 min |
| Human gut (408170) | ~5,900 | 90–120 min |

---

## 4. Step 2 — Download SRA FASTQ

**Script:** `download_sra.py`

Takes one or more BioProject accessions, finds all **METAGENOMIC** and **WGS**
runs, saves run metadata to a CSV, then downloads FASTQ files.

### Usage

```bash
# Dry run — list runs and save metadata, no download
python download_sra.py \
    --accessions PRJNA857725 \
    --email      your@email.com \
    --dry-run

# Download multiple BioProjects
python download_sra.py \
    --accessions PRJNA857725 PRJNA123456 \
    --email      your@email.com \
    --outdir     ./fastq \
    --threads    8

# If prefetch has SSL errors, bypass it
python download_sra.py \
    --accessions PRJNA857725 \
    --email      your@email.com \
    --no-prefetch
```

### Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--accessions` | Yes | — | One or more BioProject accessions |
| `--email` | Yes | — | Your email (NCBI policy) |
| `--outdir` | No | `sra_downloads` | Output directory for FASTQ files |
| `--threads` | No | `4` | Threads for `fasterq-dump` |
| `--max-runs` | No | None | Limit runs downloaded (for testing) |
| `--dry-run` | No | False | Preview without downloading |
| `--no-prefetch` | No | False | Skip `prefetch`, stream directly via `fasterq-dump` |

### Notes

- Only runs with `LibraryStrategy` = `METAGENOMIC` or `WGS` are downloaded.
- Already-downloaded runs are skipped automatically — safe to re-run.
- `sra_metadata.csv` is saved before any download starts.
- Failed runs are written to `failed_runs.txt` for easy retry.

---

## 5. Output Files

```
pig_dataset/
├── fetch_project.py
├── download_sra.py
├── README.md
│
├── projects_taxid1510822.csv      # BioProject metadata table
├── projects_taxid1510822.json
│
└── fastq/
    ├── sra_metadata.csv           # Run-level metadata for all downloaded runs
    ├── SRR21388550_1.fastq        # Forward reads (paired)
    ├── SRR21388550_2.fastq        # Reverse reads (paired)
    ├── SRR21388551.fastq          # Single-end reads
    └── failed_runs.txt            # Failed accessions (if any)
```

---

## 6. seq_type Groups

Keyword matching on project title + description:

| seq_type | Keyword matched | Typical studies |
|----------|----------------|----------------|
| `16S` | "16S" | 16S rRNA amplicon |
| `Amplicon` | "Amplicon" / "amplicon" | ITS, 18S, other amplicons |
| `RNA` | "RNA", "RNA-Seq", "metatranscriptom" | Metatranscriptomics |
| `Other` | none of the above | Shotgun metagenomics, WGS |

A project gets multiple labels (pipe-joined) if it matches more than one rule.

---

## 7. Taxon ID Reference

| Organism | Taxon ID |
|----------|----------|
| Pig gut metagenome | 1510822 |
| Human gut metagenome | 408170 |
| Mouse gut metagenome | 1268534 |
| Chicken gut metagenome | 1385655 |

Find any taxon: https://www.ncbi.nlm.nih.gov/taxonomy

---

## 8. Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| `HTTP Error 400` | Empty `api_key` string | Set to `None` or a valid key |
| `RuntimeError: EOF` | NCBI dropped connection | Script retries automatically (5×) |
| `biosample_count = 0` | Wrong elink linkname | Must use `bioproject_biosample_all` |
| `fasterq-dump not found` | SRA Toolkit not installed | `conda install -c bioconda sra-tools` |
| `IncompleteRead` | Connection reset mid-transfer | Re-run — checkpoint resumes progress |
| `PackagesNotFoundError` | Exact version not available | Remove version pin: `conda install -c bioconda sra-tools` |

### prefetch SSL error

```
mbedtls_ssl_handshake returned -9984 — Certificate verification failed
```

```bash
# Option 1: upgrade SRA Toolkit
conda install -c bioconda -c conda-forge sra-tools

# Option 2: bypass prefetch (immediate workaround)
python download_sra.py --accessions PRJNA857725 --email you@email.com --no-prefetch
```

### Retry failed downloads

```bash
while read srr; do
    fasterq-dump "$srr" --outdir ./fastq --split-files --threads 4
done < fastq/failed_runs.txt
```
