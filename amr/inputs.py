"""Inputs routines."""
from __future__ import annotations



import csv

import logging

import math

import re

from collections import Counter, defaultdict

from dataclasses import dataclass, field

import openpyxl

logger = logging.getLogger("ng_analysis")

QC = "ATCC49226"

ALIASES = {}

SEQUENCE_EXT = {".fa", ".fasta", ".fna"}

def canonical(value):
    name = str(value).strip().upper()
    if re.fullmatch(r"(QC[ _-]*)?ATCC[ _-]*49226([ _-]*(QC|CONTROL))?", name):
        return QC
    return ALIASES.get(name, name)

def assembly_id(stem):
    return canonical(re.sub(r"^\d+_(?=.)", "", stem, count=1))

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

def discover_assemblies(base, audit, mapping=None):
    result = defaultdict(dict)
    nt = base / "Assembly" / "NT"
    for kind, directory in (("annotated", nt), ("plane", nt/"Plane_Fasta")):
        groups = defaultdict(list)
        if directory.is_dir():
            for path in sorted(directory.iterdir()):
                if path.is_file() and path.suffix.lower() in SEQUENCE_EXT:
                    groups[(mapping or {}).get(path.name,assembly_id(path.stem))].append(path)
        for key, paths in groups.items():
            result[key][kind + "_files"] = paths
            if len(paths) == 1:
                result[key][kind] = str(paths[0])
            else:
                audit.add("assemblies", "collision", f"{kind} {key}: {paths}")
    return result

ANTIBIOTIC_ORDER = [
    "Spectinomycin", "Penicillin G", "Ciprofloxacin", "Azithromycin",
    "Tetracycline", "Ceftriaxone", "Cefixime", "Cefepime",
]

ANTIBIOTIC_SHORT = {
    "Spectinomycin": "SPC",
    "Penicillin G": "PG",
    "Ciprofloxacin": "CIP",
    "Azithromycin": "AZM",
    "Tetracycline": "TET",
    "Ceftriaxone": "CRO",
    "Cefixime": "CFM",
    "Cefepime": "FEP",
}

ORIGINAL, ROUNDED, AD_MIC = MIC_FIELDS = ("E-strip original MIC", "E-strip rounded MIC", "AD MIC")

EST_CAT, AD_CAT, DD_CAT, WGS_CAT = CATEGORIES = (
    "E-strip interpretation", "AD interpretation", "Disk diffusion category", "Recorded WGS prediction")

EST_MCF, AD_MCF = "E-strip McFarland", "AD McFarland"

DATE_FIELDS = ("E-strip applied", "E-strip read", "AD inoculated", "AD read", "E-strip McFarland date", "AD McFarland date")

MEASURE = re.compile(r"^\s*([<>≤≥]=?|=)?\s*(\d+(?:\.\d*)?|\.\d+)(?:[eE]([+-]?\d+))?\s*$")

def measurement(value):
    match = MEASURE.fullmatch(str(value))
    if match:
        operator, digits, exponent = match.groups()
        amount = float(digits + ("e" + exponent if exponent else ""))
        if math.isfinite(amount) and amount > 0:
            return operator or "", amount
    return "", None

def drug_name(header):
    text = re.sub("[^a-z]", "", str(header).lower())
    for drug in ANTIBIOTIC_ORDER:
        if re.sub("[^a-z]", "", drug.lower()) in text:
            return drug
    return "Penicillin G" if "penicillin" in text else None

@dataclass
class Phenotypes:
    values: dict = field(default_factory=lambda: defaultdict(dict))
    ids: dict = field(default_factory=lambda: defaultdict(set))
    trace: list = field(default_factory=list)
    breakpoints: list = field(default_factory=list)
    metadata: dict = field(default_factory=lambda: defaultdict(dict))

def is_data_id(value):
    if value is None:
        return False
    token = str(value).strip()
    if re.search(r"(?i)^(?:isolate|total|subtotal|mean|median|average|SD|variance|pearson|spearman|correlation|notes?|overall|statistics|summary|match|not.match|p.value|MIC[59]0|count|percentage|statistic|minimum|maximum|range)(?:\b|_)", token):
        return False
    return bool(re.fullmatch(r"[\w.\-]+", token) or canonical(token) == QC)

def read_master(path, audit):
    data = Phenotypes()
    formulas = openpyxl.load_workbook(path)
    cached = openpyxl.load_workbook(path, data_only=True)
    try:
        if "Breakpoints" in cached:
            bp=cached["Breakpoints"]
            # The reviewed sheet stores current-edition S/I/R in E:G and QC in H:I.
            if "36" in str(bp["E2"].value) and str(bp["A3"].value).lower()=="antibiotic":
                for cells in bp.iter_rows(min_row=4,values_only=True):
                    drug=drug_name(cells[0])
                    if drug not in ANTIBIOTIC_ORDER or str(cells[0]).strip()!=drug:continue
                    lo,mid,hi,qlo,qhi=cells[4:9]
                    def bound(v,sign):
                        return "—" if v in (None,"—") else sign+str(v)
                    data.breakpoints.append((str(len(data.breakpoints)+1),drug,ANTIBIOTIC_SHORT[drug],
                                             bound(lo,"≤"),str(mid or "—"),bound(hi,"≥"),
                                             f"{qlo}–{qhi}" if None not in (qlo,qhi) else "—"))
            else:audit.add("Breakpoints","review","unrecognized reference header; breakpoint table omitted")
        for sheet in ("EST_MIC", "AD_MIC", "McFarland", "Disc_Diffusion", "WGS"):
            if sheet not in formulas.sheetnames:
                audit.add(sheet, "failed", "required sheet missing")
                continue
            ws, cs = formulas[sheet], cached[sheet]
            head = next(((cell.row, cell.column) for row in ws.iter_rows(min_row=1,max_row=min(ws.max_row,100))
                         for cell in row if re.search(r"(?i)isolate[ _-]*id",str(cell.value))), None)
            if not head:
                audit.add(sheet, "failed", "isolate header missing")
                continue
            h, id_col = head
            merged = {(r,c):(area.min_row, area.min_col) for area in ws.merged_cells.ranges
                      if area.min_row <= h for r in range(area.min_row,min(area.max_row,h)+1)
                      for c in range(area.min_col,area.max_col+1)}
            flagged_headers=set()
            def header(row,col):
                location=merged.get((row,col),(row,col))
                value=cs.cell(*location).value
                if value is None and ws.cell(*location).data_type=="f" and location not in flagged_headers:
                    flagged_headers.add(location)
                    audit.add(sheet,"review",f"{ws.cell(*location).coordinate}: header/date formula has no cached value")
                return value
            groups = defaultdict(list)
            for row in range(h+1, ws.max_row+1):
                cell = cs.cell(row,id_col)
                if ws.cell(row,id_col).data_type=="f" and cell.value is None:
                    audit.add(sheet,"review",f"{cell.coordinate}: isolate ID formula has no cached value")
                if is_data_id(cell.value):
                    raw = str(cell.value).strip()
                    if isinstance(cell.value,int) and re.fullmatch(r"0+",cell.number_format):
                        raw = raw.zfill(len(cell.number_format))
                    groups[canonical(raw)].append((row,raw))
            columns = []
            for c in range(1,ws.max_column+1):
                title = str(header(h,c) or "")
                drug = drug_name(title)
                if sheet == "McFarland":
                    drug = next((drug_name(header(r,c)) for r in range(h-1,0,-1) if drug_name(header(r,c))),None)
                    name = {"e-strip":EST_MCF,"agar dilution":AD_MCF}.get(title.lower().strip())
                elif sheet in ("EST_MIC","AD_MIC"):
                    name = (EST_CAT if sheet=="EST_MIC" else AD_CAT) if "mic" not in title.lower() else (
                        AD_MIC if sheet=="AD_MIC" else ROUNDED if "rounded" in title.lower() else ORIGINAL)
                else:
                    name = DD_CAT if sheet=="Disc_Diffusion" else WGS_CAT
                if drug and name:
                    unit = None
                    if name in MIC_FIELDS:
                        unitmatch = re.search(r"\(([^()]*?/[^()]*)\)",title)
                        unit = unitmatch.group(1) if unitmatch else None
                        if (unit or "").lower().replace("μ","u").replace("µ","u").replace(" ","") != "ug/ml":
                            audit.add(sheet,"failed",f"unrecognized MIC unit in {ws.cell(h,c).coordinate}: {title}")
                            continue
                    columns.append((c,drug,name,unit))
            if len({(d,n) for _,d,n,_ in columns}) != len(columns):
                audit.add(sheet,"failed","ambiguous duplicate field headers")
                continue
            for key, group in groups.items():
                data.ids[key].update(raw for _,raw in group)
                if len(group) != 1:
                    audit.add(sheet, "collision", f"{key}: {group}; all rows in this sheet excluded")
                    continue
                r,raw = group[0]
                for c,drug,name,unit in columns:
                    src = ws.cell(r,c)
                    value = cs.cell(r,c).value
                    status = "recorded"
                    if src.data_type=="f" and value is None:
                        status = "formula cache unavailable or empty"
                    if cs.cell(r,c).data_type=="e":
                        status,value = "Excel error",None
                    if isinstance(value,str) and not value.strip():
                        value=None
                    if name in MIC_FIELDS and value is not None and measurement(value)[1] is None:
                        status,value="invalid MIC",None
                    if status != "recorded":
                        audit.add(sheet,"review",f"{src.coordinate}: {status}")
                    row = data.values[key,drug]
                    row[name] = value
                    data.trace.append((key,raw,drug,name,sheet,src.coordinate,value,unit,status if value is not None or status!="recorded" else "blank"))
                    for dr in range(1,h):
                        label = str(header(dr,id_col)).lower()
                        if "date" not in label: continue
                        if sheet=="McFarland":
                            datefield = DATE_FIELDS[4 if name==EST_MCF else 5]
                        elif name in (ORIGINAL,AD_MIC):
                            datefield = DATE_FIELDS[(0 if sheet=="EST_MIC" else 2) + int("reading" in label or "observation" in label)]
                        else:
                            continue
                        row[datefield] = header(dr,c)
                if sheet in ("WGS","Disc_Diffusion"):
                    for c in range(1,ws.max_column+1):
                        title=str(header(h,c) or "").lower().strip()
                        if title in ("population","study time","place","specimen type"):
                            val=cs.cell(r,c).value
                            if val is not None: data.metadata[key][title]=val
    finally:
        formulas.close()
        cached.close()
    return data

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
