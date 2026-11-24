"""Isolate identifier normalisation and alias handling."""
from __future__ import annotations

import csv
import re

from ng_amr.core.constants import (
    QC,
)

ALIASES = {}


def canonical(value):
    name = str(value).strip().upper()
    if re.fullmatch(r"(QC[ _-]*)?ATCC[ _-]*49226([ _-]*(QC|CONTROL))?", name):
        return QC
    return ALIASES.get(name, name)


def assembly_id(stem):
    return canonical(re.sub(r"^\d+_(?=.)", "", stem, count=1))


def load_aliases(path):
    """Read a one-step private mapping before workbook and filename matching."""
    aliases = {}
    if path is None:
        return aliases
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle)
        if rows.fieldnames != ["alias", "canonical_id"]:
            raise ValueError("Expected alias,canonical_id in the alias CSV")
        for number, row in enumerate(rows, 2):
            alias = (row["alias"] or "").strip().upper()
            target = (row["canonical_id"] or "").strip().upper()
            if not alias or not target or alias in aliases or None in row:
                raise ValueError(f"Empty or duplicate alias, or extra field, on row {number}")
            if alias == QC or target == QC:
                raise ValueError("The quality-control identity is reserved")
            aliases[alias] = target
    if set(aliases) & set(aliases.values()):
        raise ValueError("Map aliases directly to final canonical IDs")
    return aliases


def short_isolate_label(iso: str, max_len: int = 18) -> str:
    m = re.search(r"IsolateID[_-]?(.+?)(?:_GC[AF]\d|$)", iso or "", re.I)
    if m:
        return m.group(1).strip("_")[:max_len]
    m = re.match(r"^(\d+)", iso or "")
    return m.group(1) if m else (iso or "")[:max_len]


def normalize_sample_id(value):
    return assembly_id(value).lower()


def _parse_lab_id(value):
    return canonical(value)


def _parse_accession(full_name: str) -> str:
    m = re.search(r"(GC[AF]\d[\d.]*)", full_name, re.IGNORECASE)
    return m.group(1) if m else ""


def _parse_assembly_name(full_name: str) -> str:
    m = re.search(r"(ASM\S+)", full_name, re.IGNORECASE)
    return m.group(1) if m else ""
