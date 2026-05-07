# MGX_data — SRA Metagenomics Data Pipeline

A single-script pipeline for fetching, downloading, and organising SRA metagenomic sequencing data.

---

## Latest updates

The current `MGX_data.py` is updated from the former `MGX_data copy.py` version. Major changes:

- Replaced the old `fetch-taxon` and `fetch-project` commands with one unified `fetch-data` command.
- `fetch-data` can now collect metadata by:
  - NCBI taxon ID (`--taxon`)
  - BioProject / SRA Study accession (`--bioproject`)
  - SRA accession such as SRR, SRP, SRS, SRX, SAMN, ERR, DRR (`--accession`)
  - BioProject keyword search (`--keyword`)
- Added `--merge-into` to combine results into one CSV and deduplicate by `Run`.
- Added run-size filters: `--min-bases` and `--max-bases`.
- Improved large taxon fetching by writing partial results to disk instead of keeping all rows in memory.
- Added optional iSeq download mode with `--iseq`, `--database`, `--parallel`, and `--aspera`.
- Updated `parse` arguments from `--reference` / `--outdir` to `--source-data` / `--location`.
- `parse` can now scan downloaded folders for unknown SRRs and back-fetch missing metadata.
- The expected download layout is now `<location>/<BioProject>/<SRR>/<SRR>*.fastq.gz`.

---

## Repository layout

```text
.
├── MGX_data.py                 # Main command-line pipeline
├── README.md                   # Project overview and workflow
├── requirements.txt            # Python package dependencies
├── environment.yml             # Conda environment, including SRA Toolkit
├── pyproject.toml              # Installable Python project metadata
├── examples/
│   └── accession_example.txt   # Example two-column run/BioProject file
├── docs/
│   └── USAGE.md                # Additional usage notes
├── tests/
│   └── test_cli.py             # CLI smoke tests
└── .github/workflows/ci.yml    # GitHub Actions CI
```

## Installation

### Pip / virtualenv

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Then run either:

```bash
python MGX_data.py <command> -h
mgx-data <command> -h
```

### Conda

```bash
conda env create -f environment.yml
conda activate mgx-data
```

## Requirements

Python packages:

```bash
pip install biopython pandas requests tqdm
```

For standard SRA Toolkit downloads:

```bash
conda install -c bioconda sra-tools
```

Optional, recommended faster downloader:

```bash
conda install -c bioconda iseq
```

---

## Commands

| Command | Description |
|---------|-------------|
| `fetch-data` | Fetch SRA run metadata and BioSample attributes |
| `download` | Download FASTQ files with SRA Toolkit or iSeq |
| `parse` | Scan downloaded files and update metadata CSVs |

```bash
python MGX_data.py <command> -h
```

---

## Workflow

### Step 1 — Fetch metadata

#### Fetch by taxon ID

```bash
python MGX_data.py fetch-data \
    --taxon 1510822 \
    --email you@email.com \
    --api-key YOUR_NCBI_KEY
```

Default output:

```text
sra_taxid1510822.csv
```

#### Fetch by BioProject

```bash
python MGX_data.py fetch-data \
    --bioproject PRJNA857725 \
    --email you@email.com
```

Default output:

```text
PRJNA857725_sra_runs.csv
```

You can also pass a plain-text file with one BioProject accession per line:

```bash
python MGX_data.py fetch-data \
    --bioproject bioproject_list.txt \
    --email you@email.com
```

#### Fetch by SRA accession

```bash
python MGX_data.py fetch-data \
    --accession SRR21388550 \
    --email you@email.com
```

The `--accession` option can accept SRR, SRP, SRS, SRX, SAMN, ERR, ERS, DRR, and related SRA identifiers. It can also read a one-accession-per-line text file.

#### Fetch by keyword

```bash
python MGX_data.py fetch-data \
    --keyword "pig gut metagenome" \
    --email you@email.com
```

This searches NCBI BioProject by keyword, resolves BioProject accessions, and fetches matching WGS/METAGENOMIC SRA runs.

#### Merge multiple searches into one CSV

```bash
python MGX_data.py fetch-data \
    --taxon 1510822 \
    --bioproject PRJNA857725 \
    --merge-into source_data.csv \
    --email you@email.com
```

If the output CSV already exists, new records are appended and deduplicated by `Run`.

#### Filter by sequencing size

```bash
python MGX_data.py fetch-data \
    --taxon 1510822 \
    --min-bases 1000000000 \
    --max-bases 50000000000 \
    --email you@email.com
```

---

### Step 2 — Download

Prepare a two-column accessions file with **no header**. See `examples/accession_example.txt` for an example:

```text
SRR21388550,PRJNA857725
SRR21388551,PRJNA857725
SRR21388552,PRJNA123456
```

Tab-separated format is also accepted.

#### Standard SRA Toolkit download

```bash
python MGX_data.py download \
    --accessions runs.csv \
    --outdir ./fastq
```

#### iSeq download mode

```bash
python MGX_data.py download \
    --accessions runs.csv \
    --outdir ./fastq \
    --iseq \
    --database ena \
    --parallel 8
```

Optional download flags:

| Flag | Default | Description |
|------|---------|-------------|
| `--threads` | 8 | Threads for conversion/compression |
| `--dry-run` | — | Preview commands without downloading |
| `--iseq` | — | Use iSeq instead of SRA Toolkit |
| `--database` | `ena` | iSeq source database: `ena` or `sra` |
| `--parallel` | 8 | iSeq parallel download connections |
| `--aspera` | — | Use Aspera with iSeq |
| `--max-size` | 50 | SRA Toolkit prefetch max file size in GB |
| `--no-prefetch` | — | Skip prefetch and run fasterq-dump directly |

Downloaded files are organised as:

```text
fastq/<BioProject>/<SRR>/<SRR>*.fastq.gz
```

---

### Step 3 — Parse downloaded files

Scan downloaded files and update both the parsed metadata CSV and source-data CSV:

```bash
python MGX_data.py parse \
    --source-data source_data.csv \
    --location ./fastq \
    --parsed parsed_metadata.csv
```

To limit parsing to specific accessions:

```bash
python MGX_data.py parse \
    --accessions runs.csv \
    --source-data source_data.csv \
    --location ./fastq \
    --parsed parsed_metadata.csv
```

The parse step adds three curation columns to both output CSVs:

| Column | Description |
|--------|-------------|
| `Zlab_sort` | Manual sort order, default `1` |
| `Zlab_SRA_path` | Absolute path to the downloaded SRR folder |
| `Zlab_metadata_path` | Blank field for manual metadata path if needed |

`parsed_metadata.csv` is multi-batch safe: re-running appends new records and deduplicates by `Run`.

---

## Output file example

```text
project/
├── source_data.csv              # source metadata; updated in-place by parse
├── parsed_metadata.csv          # curated metadata for downloaded runs
└── fastq/
    ├── PRJNA857725/
    │   └── SRR21388550/
    │       ├── SRR21388550_1.fastq.gz
    │       └── SRR21388550_2.fastq.gz
    └── PRJNA123456/
        └── SRR99999999/
            └── SRR99999999.fastq.gz
```

---

## Taxon ID reference

| Organism | Taxon ID |
|----------|----------|
| Pig gut metagenome | 1510822 |
| Human gut metagenome | 408170 |
| Mouse gut metagenome | 1268534 |
| Chicken gut metagenome | 1385655 |

Find any taxon: <https://www.ncbi.nlm.nih.gov/taxonomy>
