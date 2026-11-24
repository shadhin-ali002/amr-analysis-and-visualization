"""Run audit trail, logging and CSV helpers."""
from __future__ import annotations

import csv
import logging

logger = logging.getLogger("ng_analysis")


def csv_write(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as out:
        writer = csv.writer(out)
        writer.writerow(header)
        writer.writerows(rows)


class Audit:
    def __init__(self):
        self.rows = []
    def add(self, stage, status, detail):
        self.rows.append((stage, status, str(detail)))
        logger.log(logging.ERROR if status == "failed" else logging.INFO, "%s | %s | %s", stage, status, detail)


def stage(audit,name,action):
    try:
        start=len(audit.rows)
        value=action()
        child=audit.rows[start:]
        if any(status=="failed" for _,status,_ in child):
            audit.add(name,"partial","see failed step details")
        elif child and all(status=="skipped" for _,status,_ in child):
            audit.add(name,"skipped",child[0][2])
        else:
            audit.add(name,"ok","completed")
        return value
    except Exception as exc:
        audit.add(name,"failed",f"{type(exc).__name__}: {exc}")
        logger.exception("Error in %s",name)
        return None
