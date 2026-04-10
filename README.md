# MGX_data — SRA Metagenomics Data Pipeline

A single-script pipeline for fetching, downloading, and organising SRA metagenomic sequencing data.

---

## Requirements

```bash
pip install biopython pandas requests tqdm
```

SRA Toolkit required for `download`:

```bash
conda install -c bioconda sra-tools   # v3.x or newer
```

---

## Commands

| Command | Description |
|---------|-------------|
| `fetch-taxon` | Fetch all WGS/METAGENOMIC run metadata for a taxon ID |
| `fetch-project` | Fetch full run metadata for a single BioProject |
| `download` | Download FASTQ files via prefetch + fasterq-dump |
| `parse` | Scan downloaded files and update metadata CSVs |

```bash
python MGX_data.py <command> -h   # full argument list
```

---

## Workflow

### Step 1 — Fetch metadata

**By taxon ID** (large scale, supports checkpoint/resume):

```bash
python MGX_data.py fetch-taxon \
    --taxon-id 1510822 \
    --email    you@email.com \
    --api-key  YOUR_NCBI_KEY      # optional — raises rate limit 3 → 10 req/s
```

Output: `sra_taxid1510822.csv` — all WGS/METAGENOMIC runs with BioSample attributes.

**By BioProject** (single project):

```bash
python MGX_data.py fetch-project \
    --accession PRJNA857725 \
    --email     you@email.com
```

Output: `PRJNA857725_sra_runs.csv` — all runs with BioSample attributes.

---

### Step 2 — Download

Prepare a two-column accessions file (`Run,BioProject`):

```
SRR21388550,PRJNA857725
SRR21388551,PRJNA857725
SRR21388552,PRJNA123456
```

Then download:

```bash
python MGX_data.py download \
    --accessions runs.csv \
    --outdir     ./fastq
```

Files are saved to `./fastq/<BioProject>/<SRR>_1.fastq.gz`.  
`.sra` files are removed automatically after conversion.

Optional flags:

| Flag | Default | Description |
|------|---------|-------------|
| `--threads` | 8 | fasterq-dump threads |
| `--max-size` | 50 | Max SRA file size in GB |
| `--dry-run` | — | Preview commands without downloading |

---

### Step 3 — Parse

Scan downloaded files and update both the parsed metadata CSV and the reference CSV:

```bash
python MGX_data.py parse \
    --reference sra_taxid1510822.csv \
    --outdir    ./fastq \
    --parsed    parsed_metadata.csv
```

To limit the scan to a specific accession list:

```bash
python MGX_data.py parse \
    --accessions runs.csv \
    --reference  sra_taxid1510822.csv \
    --outdir     ./fastq \
    --parsed     parsed_metadata.csv
```

The parse step adds three curation columns to both output CSVs:

| Column | Description |
|--------|-------------|
| `Zlab_sort` | Manual sort order (default `1`) |
| `Zlab_SRA_path` | Absolute path to BioProject FASTQ directory |
| `Zlab_metadata_path` | Blank — fill manually if needed |

`parsed_metadata.csv` is **multi-batch safe** — re-running appends new records and deduplicates by `Run`.

---

## Output Files

```
project/
├── sra_taxid1510822.csv        # reference metadata (updated in-place by parse)
├── parsed_metadata.csv         # curated metadata for downloaded runs
└── fastq/
    ├── PRJNA857725/
    │   ├── SRR21388550_1.fastq.gz
    │   └── SRR21388550_2.fastq.gz
    └── PRJNA123456/
        └── SRR99999999.fastq.gz
```

---

## Taxon ID Reference

| Organism | Taxon ID |
|----------|----------|
| Pig gut metagenome | 1510822 |
| Human gut metagenome | 408170 |
| Mouse gut metagenome | 1268534 |
| Chicken gut metagenome | 1385655 |

Find any taxon: <https://www.ncbi.nlm.nih.gov/taxonomy>
