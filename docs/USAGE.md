# Usage notes

This project is intentionally usable as a single script:

```bash
python MGX_data.py <command> -h
```

After editable installation, the same CLI is also available as:

```bash
mgx-data <command> -h
```

Large sequencing files and generated CSV outputs should not be committed to Git. Keep them under ignored folders such as `fastq/`, `data/`, `outputs/`, or `results/`.
