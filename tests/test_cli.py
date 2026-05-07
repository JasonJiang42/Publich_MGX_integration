import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "MGX_data.py"


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def test_top_level_help():
    result = run_cli("-h")
    assert result.returncode == 0
    assert "fetch-data" in result.stdout
    assert "download" in result.stdout
    assert "parse" in result.stdout


def test_subcommand_help():
    for command in ["fetch-data", "download", "parse"]:
        result = run_cli(command, "-h")
        assert result.returncode == 0
        assert "usage:" in result.stdout
