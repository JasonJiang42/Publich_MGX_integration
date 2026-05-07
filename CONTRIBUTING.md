# Contributing

Thank you for helping improve this project.

## Local setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .[dev]
```

## Checks

```bash
python -m py_compile MGX_data.py
pytest -q
ruff check .
```

## Pull requests

Please keep changes focused and include a short description of the workflow or bug that was tested.
