"""Mutation x isolate master matrices and the combined master workbook."""
from __future__ import annotations

import os
import re
from collections import Counter
from typing import List

import openpyxl
import pandas as pd

from ng_amr.core.audit import (
    logger,
)
from ng_amr.core.constants import (
    ANTIBIOTIC_ORDER,
    GENE_SHEET_ORDER,
    MIC_FIELDS,
    QC,
)
from ng_amr.core.ids import (
    _parse_accession,
    _parse_assembly_name,
    _parse_lab_id,
)


def _mutation_sort_key(notation: str):
    m = re.search(r"(\d+)", notation or "")
    pos = int(m.group(1)) if m else 0
    return (pos, notation or "")

def write_master_matrix(gene_results: List[dict], output_root: str, thresholds: List[int], book):
    if pd is None:
        logger.info("  pandas missing – skip master matrix")
        return None

    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font, PatternFill
    except ImportError:
        logger.info("  openpyxl missing – skip multi-sheet master matrix")
        return None

    isolate_set = set()
    for gr in gene_results:
        for d in gr.get("all_data") or []:
            isolate_set.add(d["isolate"])
        for d in gr.get("all_mutation_data") or []:
            isolate_set.add(d["isolate"])
    def _iso_order(name: str):
        m = re.match(r"^(\d+)", (name or "").strip())
        idx = int(m.group(1)) if m else 10**9
        return (idx, _parse_lab_id(name), name)
    isolates = sorted(isolate_set, key=_iso_order)
    n_isolates = len(isolates)
    if n_isolates == 0:
        logger.info("  No isolates – skip master matrix")
        return None

    lab_ids = {iso: _parse_lab_id(iso) for iso in isolates}

    gene_iso = {}
    gene_muts = {}
    gene_mut_meta = {}

    for gr in gene_results:
        gene = gr["gene"]
        curated_set = gr.get("curated_set") or set()
        iso_info = {}
        mut_meta = {}
        present_counts = {}

        data_by_iso = {d["isolate"]: d for d in (gr.get("all_data") or [])}
        mutdata_by_iso = {d["isolate"]: set(d.get("mutations") or [])
                          for d in (gr.get("all_mutation_data") or [])}

        for iso in isolates:
            ad = data_by_iso.get(iso)
            if ad is None:
                iso_info[iso] = {"found": False, "pident": None, "muts": set()}
                continue
            found = bool(ad.get("gene_found", ad.get("status") == "OK"))
            muts = mutdata_by_iso.get(iso, set())
            if found:
                pass
            iso_info[iso] = {
                "found": found,
                "pident": ad.get("pident") if found else None,
                "muts": muts if found else set(),
            }
            if found:
                for m in muts:
                    present_counts[m] = present_counts.get(m, 0) + 1

            for r in ad.get("mutations") or []:
                notation = r.get("mutation") or r.get("notation")
                if not notation:
                    continue
                if notation not in mut_meta:
                    mut_meta[notation] = {
                        "type": r.get("type", "substitution"),
                        "position": r.get("position", 0),
                        "curated": bool(r.get("curated", notation in curated_set)),
                    }

        ordered = sorted(present_counts.keys(), key=_mutation_sort_key)
        for m in ordered:
            if m not in mut_meta:
                mut_meta[m] = {
                    "type": "substitution",
                    "position": _mutation_sort_key(m)[0],
                    "curated": m in curated_set,
                }
            mut_meta[m]["count"] = present_counts.get(m, 0)

        gene_iso[gene] = iso_info
        gene_muts[gene] = ordered
        gene_mut_meta[gene] = mut_meta

    genes_present = [g for g in GENE_SHEET_ORDER if g in gene_iso]
    for g in gene_iso:
        if g not in genes_present:
            genes_present.append(g)

    flat_cols = []
    for g in genes_present:
        for m in gene_muts[g]:
            flat_cols.append((g, m))

    wb = Workbook()

    header_fill = PatternFill("solid", fgColor="1F4E79")
    header_font = Font(color="FFFFFF", bold=True, size=10)
    meta_fill = PatternFill("solid", fgColor="D6EAF8")
    yes_fill = PatternFill("solid", fgColor="C6EFCE")
    no_fill = PatternFill("solid", fgColor="FFFFFF")
    na_fill = PatternFill("solid", fgColor="FCE4D6")
    mut_present_fill = PatternFill("solid", fgColor="F4B183")
    white_fill = PatternFill("solid", fgColor="FFFFFF")
    center = Alignment(horizontal="center", vertical="center", wrap_text=True)

    def _style_header_row(ws, row, start_col, end_col):
        for c in range(start_col, end_col + 1):
            cell = ws.cell(row, c)
            cell.fill = header_fill
            cell.font = header_font
            cell.alignment = center

    mic_columns = [(d, f) for d in ANTIBIOTIC_ORDER for f in MIC_FIELDS]
    mic_headers = [f"{d} | {f} (µg/mL)" for d, f in mic_columns]

    ws = wb.active
    ws.title = "MIC_Mutation_Matrix"
    meta_labels = [
        "Gene", "Mutation", "Curated", "Type", "Position",
        "Isolates with Mutation", "Isolates Tested", "Prevalence %",
    ]
    n_mic = len(mic_headers)
    first_mut_col = 2 + n_mic

    for i, lab in enumerate(meta_labels, start=1):
        ws.cell(i, 1, lab)

    ws.cell(9, 1, "Isolate ID")
    for j, drug in enumerate(mic_headers):
        ws.cell(9, 2 + j, drug)

    for j, (g, m) in enumerate(flat_cols):
        col = first_mut_col + j
        meta = gene_mut_meta[g][m]
        cnt = meta.get("count", 0)
        prev = round(100.0 * cnt / n_isolates, 2) if n_isolates else 0.0
        ws.cell(1, col, g)
        ws.cell(2, col, m)
        curated_cell = ws.cell(3, col, "TRUE" if meta.get("curated") else "FALSE")
        if meta.get("curated"):
            curated_cell.fill = yes_fill
        ws.cell(4, col, meta.get("type", "substitution"))
        ws.cell(5, col, meta.get("position", 0))
        ws.cell(6, col, cnt)
        ws.cell(7, col, n_isolates)
        prev_str = f"{int(prev)}%" if prev == int(prev) else f"{prev}%"
        ws.cell(8, col, prev_str)
        ws.cell(9, col, f"{g}_{m}")

    for i in range(1, 9):
        for c in range(1, first_mut_col + len(flat_cols)):
            cell = ws.cell(i, c)
            cell.fill = meta_fill
            cell.alignment = center
    for j, (g, m) in enumerate(flat_cols):
        col = first_mut_col + j
        meta = gene_mut_meta[g][m]
        if meta.get("curated"):
            cell = ws.cell(3, col)
            cell.fill = yes_fill
            cell.font = Font(color="006100", bold=True, size=10)
    _style_header_row(ws, 9, 1, first_mut_col + len(flat_cols) - 1)

    for ri, iso in enumerate(isolates):
        row = 10 + ri
        cell = ws.cell(row, 1, lab_ids[iso])
        cell.alignment = center
        cell.fill = white_fill
        for j, drug in enumerate(mic_headers):
            mic_val = book.values.get((iso, mic_columns[j][0]), {}).get(mic_columns[j][1])
            cell = ws.cell(row, 2 + j, mic_val if mic_val is not None else None)
            cell.alignment = center
            cell.fill = white_fill
        for j, (g, m) in enumerate(flat_cols):
            col = first_mut_col + j
            info = gene_iso[g][iso]
            val = int(m in info["muts"]) if info["found"] else None
            cell = ws.cell(row, col, val)
            cell.alignment = center
            if val == 1:
                cell.fill = mut_present_fill
            else:
                cell.fill = white_fill

    ws = wb.create_sheet("Mutation_Summary")
    for i, lab in enumerate(meta_labels, start=1):
        ws.cell(i, 1, lab)

    first_sum_mut_col = 3
    for j, (g, m) in enumerate(flat_cols):
        col = first_sum_mut_col + j
        meta = gene_mut_meta[g][m]
        cnt = meta.get("count", 0)
        prev = round(100.0 * cnt / n_isolates, 2) if n_isolates else 0.0
        prev_str = f"{int(prev)}%" if prev == int(prev) else f"{prev}%"
        ws.cell(1, col, g)
        ws.cell(2, col, m)
        curated_cell = ws.cell(3, col, "TRUE" if meta.get("curated") else "FALSE")
        if meta.get("curated"):
            curated_cell.fill = yes_fill
            curated_cell.font = Font(color="006100", bold=True, size=10)
        ws.cell(4, col, meta.get("type", "substitution"))
        ws.cell(5, col, meta.get("position", 0))
        ws.cell(6, col, cnt)
        ws.cell(7, col, n_isolates)
        ws.cell(8, col, prev_str)

    for i in range(1, 9):
        for c in range(1, first_sum_mut_col + len(flat_cols)):
            cell = ws.cell(i, c)
            cell.fill = meta_fill
            cell.alignment = center
    for j, (g, m) in enumerate(flat_cols):
        col = first_sum_mut_col + j
        if gene_mut_meta[g][m].get("curated"):
            cell = ws.cell(3, col)
            cell.fill = yes_fill
            cell.font = Font(color="006100", bold=True, size=10)

    ws.cell(9, 1, "Isolate_ID")
    ws.cell(9, 2, "Isolate_Full_Name")
    for j, (g, m) in enumerate(flat_cols):
        ws.cell(9, first_sum_mut_col + j, f"{g}_{m}")
    _style_header_row(ws, 9, 1, first_sum_mut_col + len(flat_cols) - 1)

    for ri, iso in enumerate(isolates):
        row = 10 + ri
        ws.cell(row, 1, lab_ids[iso])
        ws.cell(row, 2, iso)
        for j, (g, m) in enumerate(flat_cols):
            info = gene_iso[g][iso]
            if not info["found"]:
                val = "N/A"
            elif m in info["muts"]:
                val = "YES"
            else:
                val = "NO"
            cell = ws.cell(row, first_sum_mut_col + j, val)
            if val == "YES":
                cell.fill = yes_fill
            elif val == "N/A":
                cell.fill = na_fill

    ws = wb.create_sheet("Mutation_Matrix")
    for i, lab in enumerate(meta_labels, start=1):
        ws.cell(i, 1, lab)

    first_mm_mut_col = 2 + n_mic
    for j, (g, m) in enumerate(flat_cols):
        col = first_mm_mut_col + j
        meta = gene_mut_meta[g][m]
        cnt = meta.get("count", 0)
        prev = round(100.0 * cnt / n_isolates, 2) if n_isolates else 0.0
        prev_str = f"{int(prev)}%" if prev == int(prev) else f"{prev}%"
        ws.cell(1, col, g)
        ws.cell(2, col, m)
        curated_cell = ws.cell(3, col, "TRUE" if meta.get("curated") else "FALSE")
        if meta.get("curated"):
            curated_cell.fill = yes_fill
            curated_cell.font = Font(color="006100", bold=True, size=10)
        ws.cell(4, col, meta.get("type", "substitution"))
        ws.cell(5, col, meta.get("position", 0))
        ws.cell(6, col, cnt)
        ws.cell(7, col, n_isolates)
        ws.cell(8, col, prev_str)

    for i in range(1, 9):
        for c in range(1, first_mm_mut_col + len(flat_cols)):
            cell = ws.cell(i, c)
            cell.fill = meta_fill
            cell.alignment = center
    for j, (g, m) in enumerate(flat_cols):
        col = first_mm_mut_col + j
        if gene_mut_meta[g][m].get("curated"):
            cell = ws.cell(3, col)
            cell.fill = yes_fill
            cell.font = Font(color="006100", bold=True, size=10)

    ws.cell(9, 1, "Isolate_ID")
    for j, drug in enumerate(mic_headers):
        ws.cell(9, 2 + j, drug)
    for j, (g, m) in enumerate(flat_cols):
        ws.cell(9, first_mm_mut_col + j, f"{g}_{m}")
    _style_header_row(ws, 9, 1, first_mm_mut_col + len(flat_cols) - 1)

    for ri, iso in enumerate(isolates):
        row = 10 + ri
        ws.cell(row, 1, lab_ids[iso])
        for j in range(n_mic):
            ws.cell(row, 2 + j, book.values.get((iso, mic_columns[j][0]), {}).get(mic_columns[j][1]))
        for j, (g, m) in enumerate(flat_cols):
            info = gene_iso[g][iso]
            if not info["found"]:
                val = "N/A"
            elif m in info["muts"]:
                val = "YES"
            else:
                val = "NO"
            cell = ws.cell(row, first_mm_mut_col + j, val)
            if val == "YES":
                cell.fill = yes_fill
            elif val == "N/A":
                cell.fill = na_fill

    for g in genes_present:
        muts = gene_muts[g]
        meta = gene_mut_meta[g]
        ws = wb.create_sheet(g[:31])

        meta_row_labels = [
            (1, "Mutation"),
            (2, "Type"),
            (3, "Position"),
            (4, "Curated"),
            (5, "Isolates with Mutation"),
            (6, "Isolates Tested"),
            (7, "Prevalence %"),
        ]
        for r, lab in meta_row_labels:
            ws.cell(r, 1, lab)

        for j, m in enumerate(muts):
            col = 3 + j
            cnt = meta[m].get("count", 0)
            prev = round(100.0 * cnt / n_isolates, 2) if n_isolates else 0.0
            prev_str = f"{int(prev)}%" if prev == int(prev) else f"{prev}%"
            ws.cell(1, col, m)
            ws.cell(2, col, meta[m].get("type", "substitution"))
            ws.cell(3, col, meta[m].get("position", 0))
            ws.cell(4, col, "TRUE" if meta[m].get("curated") else "FALSE")
            ws.cell(5, col, cnt)
            ws.cell(6, col, n_isolates)
            ws.cell(7, col, prev_str)

        for r in range(1, 8):
            for c in range(1, 3 + len(muts)):
                cell = ws.cell(r, c)
                cell.fill = meta_fill
                cell.alignment = center
        for j, m in enumerate(muts):
            if meta[m].get("curated"):
                cell = ws.cell(4, 3 + j)
                cell.fill = yes_fill
                cell.font = Font(color="006100", bold=True, size=10)

        ws.cell(8, 1, "Isolate")
        ws.cell(8, 2, "% Identity")
        for j, m in enumerate(muts):
            ws.cell(8, 3 + j, m)
        _style_header_row(ws, 8, 1, 2 + len(muts))

        for ri, iso in enumerate(isolates):
            row = 9 + ri
            info = gene_iso[g][iso]
            ws.cell(row, 1, iso)
            if not info["found"]:
                ws.cell(row, 2, "N/A")
                for j in range(len(muts)):
                    cell = ws.cell(row, 3 + j, "N/A")
                    cell.fill = na_fill
            else:
                pident = info["pident"]
                has_any = bool(info["muts"] & set(muts))
                if pident in ("", None) and not has_any:
                    ws.cell(row, 2, "N/R")
                else:
                    ws.cell(row, 2, pident if pident not in ("", None) else "N/R")
                for j, m in enumerate(muts):
                    if m in info["muts"]:
                        cell = ws.cell(row, 3 + j, "YES")
                        cell.fill = yes_fill
                    else:
                        cell = ws.cell(row, 3 + j, "NO")
                        cell.fill = no_fill

    ws = wb.create_sheet("Gene_Call_Status")
    headers = ["Isolate_ID", "Isolate_Full_Name"] + genes_present
    for c, h in enumerate(headers, 1):
        ws.cell(1, c, h)
    _style_header_row(ws, 1, 1, len(headers))
    for ri, iso in enumerate(isolates):
        row = 2 + ri
        ws.cell(row, 1, lab_ids[iso])
        ws.cell(row, 2, iso)
        for j, g in enumerate(genes_present):
            found = gene_iso[g][iso]["found"]
            cell = ws.cell(row, 3 + j, "Yes" if found else "No")
            if found:
                cell.fill = yes_fill
            else:
                cell.fill = na_fill

    ws = wb.create_sheet("Isolate_Metadata")
    headers = ["Isolate_Index", "Isolate_Full_Name", "Isolate_Lab_ID",
               "Assembly_Accession", "Assembly_Name"]
    for c, h in enumerate(headers, 1):
        ws.cell(1, c, h)
    _style_header_row(ws, 1, 1, len(headers))
    for ri, iso in enumerate(isolates):
        row = 2 + ri
        ws.cell(row, 1, ri + 1)
        ws.cell(row, 2, iso)
        ws.cell(row, 3, lab_ids[iso])
        ws.cell(row, 4, _parse_accession(iso))
        ws.cell(row, 5, _parse_assembly_name(iso))

    ws = wb.create_sheet("Key_Notes")
    notes = [
        ("Isolate_Metadata",
         f"One row per isolate (n={n_isolates}). Parsed identifiers for joining to external metadata or MIC tables."),
        ("Mutation_Summary",
         "Rows 1-8: mutation metadata (Gene, Mutation, Curated, Type, Position, Isolates with Mutation, Isolates Tested, Prevalence %); row 9 headers; row 10+: isolate x mutation with YES/NO/N/A."),
        ("Mutation_Matrix",
         "Same 8 metadata rows as MIC_Mutation_Matrix, then Isolate_ID + measured MICs + Gene_Mutation columns with YES/NO/N/A."),
        ("MIC_Mutation_Matrix",
         "Transposed matrix: mutation metadata in header rows, isolates as rows, 0/1 presence for called loci; method-specific MIC columns."),
        ("Gene_Call_Status",
         "One row per isolate flagging whether each gene returned a usable call (Yes/No)."),
        ("Per-gene sheets",
         "Per gene: metadata rows (Mutation, Type, Position, Curated, Isolates with Mutation, Isolates Tested, Prevalence %) then isolate x mutation with % Identity; YES/NO/N/A."),
        ("", ""),
        ("Column / Value Definitions", ""),
        ("Gene", "AMR determinant gene name."),
        ("Mutation", "Mutation code, e.g. 'F504L' (substitution) or 'ins346D' (insertion)."),
        ("Type", "substitution, insertion, deletion, stop_gain, frameshift, etc."),
        ("Position", "Amino acid (or nucleotide) position in the reference."),
        ("Curated", "TRUE = established resistance marker in curated_variants.json; FALSE = other detected change."),
        ("YES / NO / N/A",
         "YES = mutation present; NO = gene called and wildtype at that site; N/A = gene not called."),
        ("% Identity = N/R",
         "N/R = gene called but no mutations reported (wildtype); N/A = gene not called."),
        ("", ""),
        ("Prevalence Definitions", ""),
        ("Prevalence % (MIC sheet)",
         "N_Isolates_With_Mutation / N_Total_Isolates (fixed panel denominator)."),
        ("", ""),
        ("Notes", ""),
        ("MIC columns",
         "Cefepime, Cefixime, Ceftriaxone, Ciprofloxacin, Penicillin G, Tetracycline, Azithromycin, Spectinomycin – original and rounded E-strip and agar dilution in µg/mL."),
        ("No prevalence filter",
         "Every mutation detected by the pipeline is included (no 25/30/35/40/50% cutoff in this workbook)."),
        ("Genes included",
         ", ".join(genes_present)),
    ]
    for r, (a, b) in enumerate(notes, 1):
        ws.cell(r, 1, a)
        ws.cell(r, 2, b)
        if a and not b:
            ws.cell(r, 1).font = Font(bold=True, size=11)

    for ws in wb.worksheets:
        for row in ws.iter_rows(min_row=1, max_row=ws.max_row,
                                min_col=1, max_col=ws.max_column):
            for cell in row:
                cell.alignment = center

    out_path = os.path.join(output_root, "NG_AMR_Mutation_MIC_Matrix.xlsx")
    wb.save(out_path)
    logger.info(f"\nMaster mutation matrix (multi-sheet) → {out_path}")
    logger.info(f"  Isolates={n_isolates}  Genes={len(genes_present)}  Mutations={len(flat_cols)}")

    return out_path


def combined_master(root):
    """Combine the two workbooks produced in this run. Inputs never come from a search."""
    from copy import copy
    mutation=root/"NG_AMR_Mutation_MIC_Matrix.xlsx"
    phenotype=root/"NG_AMR_MIC.xlsx"
    if not mutation.is_file() or not phenotype.is_file():
        raise FileNotFoundError("the current run needs both the mutation matrix and phenotype workbook")
    combined=openpyxl.load_workbook(mutation)
    source=openpyxl.load_workbook(phenotype)
    try:
        for sheet in source:
            if sheet.title in combined.sheetnames:raise ValueError(f"duplicate output sheet: {sheet.title}")
            dest=combined.create_sheet(sheet.title)
            for row in sheet:
                dest.append([cell.value for cell in row])
            for cell in sheet[1]:
                out=dest[cell.coordinate]
                out.font=copy(cell.font)
                out.fill=copy(cell.fill)
                out.alignment=copy(cell.alignment)
            dest.freeze_panes=sheet.freeze_panes
            dest.auto_filter.ref=sheet.auto_filter.ref
            dest.row_dimensions[1].height=sheet.row_dimensions[1].height
            for name,dim in sheet.column_dimensions.items():dest.column_dimensions[name].width=dim.width
        combined.save(root/"NG_AMR_MIC_MasterMatrix.xlsx")
    finally:
        combined.close()
        source.close()


def clinical_view(results, selected=None):
    view=[]
    for gr in results:
        rows=[d for d in gr["all_data"] if d["isolate"]!=QC and (selected is None or d["isolate"] in selected)]
        called=[d for d in rows if d.get("gene_found")]
        mutations=[{"isolate":d["isolate"],"mutations":[m["mutation"] for m in d["mutations"]]} for d in called]
        counts=Counter(m for r in mutations for m in set(r["mutations"]))
        view.append(dict(gr,all_data=rows,all_mutation_data=mutations,processed=len(called),
                         skipped=len(rows)-len(called),n_found=len(called),total_isolates=len(rows),
                         n_collection=len(rows),mutation_counts=dict(counts)))
    return view
