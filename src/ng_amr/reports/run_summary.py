"""Run manifest: arguments, timings, stage outcomes and source checksums."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from ng_amr import __version__


def report_run(root, args, results, audit, started):
    from datetime import datetime
    now = datetime.now().astimezone()
    package = Path(__file__).resolve().parents[1]
    manifest = dict(arguments=vars(args), started=started.isoformat(), finished=now.isoformat(),
                    python=sys.version, version=__version__,
                    source_sha256={p.relative_to(package).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sorted(package.rglob("*.py"))},
                    workbook_sha256=hashlib.sha256(args.workbook.read_bytes()).hexdigest() if args.workbook is not None else None,
                    stages=audit.rows)
    other = root / "05_Other_Outputs"
    (other / "NG_AMR_Run_Manifest.json").write_text(json.dumps(manifest, default=str, indent=2), encoding="utf-8")
