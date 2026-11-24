"""Run ng-amr from a source checkout without installing it (same as `ng-amr` / `python -m ng_amr`)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from ng_amr.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
