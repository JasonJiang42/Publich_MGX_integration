# Gut Microbiome — NCBI Data Download Pipeline

Scripts to fetch SRA run metadata from NCBI and download FASTQ files.

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Quick Start](#2-quick-start)
3. [Workflow A — By Taxon ID (large-scale)](#3-workflow-a--by-taxon-id-large-scale)
4. [Workflow B — By BioProject (single project)](#4-workflow-b--by-bioproject-single-project)
5. [Step 3 — Download & Parse](#5-step-3--download--parse)
6. [Output Files](#6-output-files)
7. [Taxon ID Reference](#7-taxon-id-reference)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Prerequisites

### Python packages

```bash
pip install biopython pandas requests tqdm
```

### SRA Toolkit

Required for `03_download_parse_prefetch.py`.

```bash
# conda (recommended)
conda install -c bioconda sra-tools

# macOS Homebrew
brew install sratoolkit
```

> Use SRA Toolkit **v3.x or newer**. Older versions have SSL issues with `prefetch`.

Run once after installation to accept defaults:

```bash
vdb-config -i
```

---

## 2. Quick Start

```bash
# Workflow A — fetch metadata by taxon ID, then download
python sra_fetch.py --taxon-id 1510822 --email your@email.com
python 03_download_parse_prefetch.py \
    --accessions runs.txt \
    --reference  sra_taxid1510822.csv \
    --outdir     ./fastq \
    --parsed     parsed_metadata.csv \
    --curator    "Jason"

# Workflow B — fetch metadata for one BioProject, then download
python 01_fetch_sra_metadata.py --accession PRJNA857725 --email your@email.com
python 03_download_parse_prefetch.py \
    --accessions runs.txt \
    --reference  PRJNA857725_sra_runs.csv \
    --outdir     ./fastq \
    --parsed     parsed_metadata.csv \
    --curator    "Jason"
```

---

## 3. Workflow A — By Taxon ID (large-scale)

**Script:** `sra_fetch.py`

Searches the entire SRA database by NCBI Taxon ID. Returns all runs with
`LibraryStrategy = WGS` or `METAGENOMIC`, merged with BioSample attributes.

### Usage

```bash
python sra_fetch.py \
    --taxon-id 1510822 \
    --email    your@email.com

# With NCBI API key (raises rate limit 3 → 10 req/s, recommended for large taxons)
python sra_fetch.py \
    --taxon-id 1510822 \
    --email    your@email.com \
    --api-key  YOUR_NCBI_API_KEY
```

### Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--taxon-id` | Yes | — | NCBI Taxon ID |
| `--email` | Yes | — | Your email (NCBI policy) |
| `--label` | No | `sra_taxid<ID>` | Output filename prefix |
| `--api-key` | No | None | NCBI API key for higher rate limits |

### Output

| File | Description |
|------|-------------|
| `sra_taxid1510822.csv` | All WGS/METAGENOMIC runs with BioSample attributes |

### Notes

- **Checkpoint/resume:** saved to `<label>_checkpoint.json` after every batch.
  Re-run the same command to resume after a crash.
- Filtered to `WGS` and `METAGENOMIC` strategies only.
- BioSample attributes (geo_loc_name, host, collection_date, etc.) are merged automatically.

---

## 4. Workflow B — By BioProject (single project)

**Script:** `01_fetch_sra_metadata.py`

Fetches the full SRA run metadata for a single BioProject or SRA Study accession,
matching the output of the NCBI SRA Run Selector page.

### Usage

```bash
python 01_fetch_sra_metadata.py \
    --accession PRJNA857725 \
    --email     your@email.com

# Custom output filename
python 01_fetch_sra_metadata.py \
    --accession PRJNA857725 \
    --email     your@email.com \
    --out       pig_gut_runs.csv
```

### Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--accession` | Yes | — | BioProject (PRJNA…) or SRA Study (SRP…) accession |
| `--email` | No | None | Your email (recommended for NCBI) |
| `--out` | No | `<accession>_sra_runs.csv` | Output CSV filename |

### Output columns (selected)

| Column | Description |
|--------|-------------|
| `Run` | SRR accession |
| `BioProject` | BioProject accession |
| `BioSample` | BioSample accession |
| `LibraryStrategy` | WGS, METAGENOMIC, AMPLICON, etc. |
| `LibraryLayout` | PAIRED or SINGLE |
| `Platform` | ILLUMINA, PACBIO, etc. |
| `spots` | Number of reads |
| `bases` | Total base pairs |
| `geo_loc_name` | Collection location (from BioSample) |
| `host` | Host organism (from BioSample) |
| `collection_date` | Sample collection date (from BioSample) |

### Notes

- Equivalent to downloading the metadata table from
  `https://www.ncbi.nlm.nih.gov/Traces/study/`
- Includes all BioSample attributes (geo_loc_name, host, collection_date, etc.).

### Prepare accession list for download

After fetching metadata, create the `runs.txt` input for the download step:

```bash
# All runs in the project
tail -n +2 PRJNA857725_sra_runs.csv | cut -d',' -f1 > runs.txt

# WGS/METAGENOMIC runs only
python -c "
import pandas as pd
df = pd.read_csv('PRJNA857725_sra_runs.csv')
mask = df['LibraryStrategy'].str.upper().isin({'WGS','METAGENOMIC'})
df.loc[mask, 'Run'].to_csv('runs.txt', index=False, header=False)
"
```

---

## 5. Step 3 — Download & Parse

**Script:** `03_download_parse_prefetch.py`

Downloads SRA runs using `prefetch` + `fasterq-dump`, verifies files on disk,
and builds a curated metadata CSV. Designed for multi-batch use — the output
parsed file is safely updated across multiple runs.

### Usage

```bash
# Preview commands without downloading
python 03_download_parse_prefetch.py \
    --accessions runs.txt \
    --reference  sra_taxid1510822.csv \
    --outdir     ./fastq \
    --parsed     parsed_metadata.csv \
    --curator    "Jason" \
    --dry-run

# Download
python 03_download_parse_prefetch.py \
    --accessions runs.txt \
    --reference  sra_taxid1510822.csv \
    --outdir     ./fastq \
    --parsed     parsed_metadata.csv \
    --curator    "Jason" \
    --threads    8 \
    --max-size   100
```

### Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--accessions` | Yes | — | Text file with one SRR accession per line |
| `--reference` | Yes | — | Reference CSV from `sra_fetch.py` or `01_fetch_sra_metadata.py` |
| `--outdir` | Yes | — | Root output directory |
| `--parsed` | No | `parsed_metadata.csv` | Curated output CSV (appended across batches) |
| `--curator` | Yes | — | Your name, recorded in `Zlab_curator` column |
| `--threads` | No | `8` | Threads for `fasterq-dump` |
| `--max-size` | No | `50` | Max SRA file size in GB for `prefetch` |
| `--dry-run` | No | False | Print commands without downloading |

### Download structure

```
fastq/
└── PRJNA857725/
    ├── SRR21388550_1.fastq.gz    # paired forward
    ├── SRR21388550_2.fastq.gz    # paired reverse
    └── SRR21388551.fastq.gz      # single-end
```

### Pipeline steps

| Step | Description |
|------|-------------|
| 1 | Load reference CSV + accession list |
| 2 | Group SRRs by BioProject; run `prefetch` → `fasterq-dump` → `gzip` |
| 3 | Verify FASTQ files exist on disk; write `failed_runs.txt` if any failed |
| 4 | Extract metadata for confirmed downloads from reference CSV |
| 5 | Drop all-blank columns; add 5 Zlab curation columns |
| 6 | Append to parsed CSV (deduplication by Run, multi-batch safe) |
| 7 | Update reference CSV in-place with Zlab columns |

### Zlab columns added

| Column | Description |
|--------|-------------|
| `Zlab_sort` | Manual sort order (default `1`) |
| `Zlab_SRA_path` | Local path to BioProject FASTQ directory |
| `Zlab_metadata_path` | Path to additional metadata (fill manually if needed) |
| `Zlab_curator` | Curator name (from `--curator`) |
| `Zlab_curated_date` | Date of download (ISO format, auto-filled) |

### Multi-batch usage

Run the script multiple times with different `runs.txt` files.
The parsed file accumulates — each batch's records are merged by `Run` accession.

```bash
python 03_download_parse_prefetch.py --accessions batch1.txt ...
python 03_download_parse_prefetch.py --accessions batch2.txt ...
# parsed_metadata.csv now contains all confirmed runs from both batches
```

### Retry failed downloads

```bash
# Re-run using the failed list as the next accession input
python 03_download_parse_prefetch.py \
    --accessions fastq/failed_runs.txt \
    --reference  sra_taxid1510822.csv \
    --outdir     ./fastq \
    --parsed     parsed_metadata.csv \
    --curator    "Jason"
```

---

## 6. Output Files

```
pig_dataset/
├── sra_fetch.py
├── 01_fetch_sra_metadata.py
├── 03_download_parse_prefetch.py
│
├── sra_taxid1510822.csv          # All WGS/METAGENOMIC runs by taxon ID
├── PRJNA857725_sra_runs.csv      # Runs for a single BioProject (Workflow B)
│
├── parsed_metadata.csv           # Curated metadata for all downloaded runs
│                                 # (appended across batches, includes Zlab columns)
│
└── fastq/
    ├── failed_runs.txt           # Failed SRR accessions (if any)
    ├── PRJNA857725/
    │   ├── SRR21388550_1.fastq.gz
    │   ├── SRR21388550_2.fastq.gz
    │   └── SRR21388551.fastq.gz
    └── PRJNA123456/
        └── SRR99999999.fastq.gz
```

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
| `HTTP Error 400` in sra_fetch.py | Empty `api_key` string passed | Set to `None` or provide a valid key |
| `RuntimeError: EOF` | NCBI dropped connection mid-request | Script retries automatically (5×) |
| `IncompleteRead` | Connection reset during transfer | Re-run — checkpoint resumes progress |
| `fasterq-dump not found` | SRA Toolkit not installed | `conda install -c bioconda sra-tools` |
| `PackagesNotFoundError` | Exact version unavailable | Remove version pin: `conda install -c bioconda sra-tools` |
| `No data returned` in 01 script | Invalid or private accession | Check accession exists and has public data |

### prefetch SSL error

```
mbedtls_ssl_handshake returned -9984 — Certificate verification failed
```

Cause: SRA Toolkit older than v3.x uses an outdated TLS library.

```bash
# Upgrade SRA Toolkit
conda install -c bioconda sra-tools
```

### prefetch file too large

```
This file is larger than the maximum allowed: 20,971,520,000
```

```bash
# Increase the limit with --max-size (in GB)
python 03_download_parse_prefetch.py ... --max-size 100
```

### fasterq-dump temp space

`fasterq-dump` writes large temporary files during conversion. If your disk is full:

```bash
# Point temp files to a larger disk
fasterq-dump SRR.sra --outdir ./fastq --temp /path/to/large/disk/tmp
```

For batch use, set the `--temp` flag by editing `FASTERQ_THREADS` / `cmd_fq` in
`03_download_parse_prefetch.py` line ~115.
