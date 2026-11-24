"""Shared helpers: make `ng_amr` importable without installation and run the CLI."""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
FIXTURES = ROOT / "tests" / "fixtures"
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))


def run_cli(*args, cwd):
    """Run the command line through `run_amr.py`; return (returncode, stdout, stderr)."""
    result = subprocess.run([sys.executable, str(ROOT / "run_amr.py"), *map(str, args)],
                            cwd=cwd, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr


def run_folder(stdout: str) -> Path:
    """The run folder is the last line the CLI prints."""
    return Path(stdout.strip().splitlines()[-1])


def read_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))
