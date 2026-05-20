import subprocess
import sys
from pathlib import Path

from MGX_data import parse_strategy_values


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


def test_fetch_data_strategy_help():
    result = run_cli("fetch-data", "-h")
    assert result.returncode == 0
    assert "--strategy" in result.stdout
    assert "AMPLICON" in result.stdout


def test_parse_strategy_values():
    assert parse_strategy_values(["16S"]) == {"AMPLICON"}
    assert parse_strategy_values(["WGS,AMPLICON"]) == {"WGS", "AMPLICON"}
    assert parse_strategy_values(["all"]) is None
