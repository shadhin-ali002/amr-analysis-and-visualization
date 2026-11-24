"""Comparison of called proteins with supplied original protein FASTA files."""
from __future__ import annotations

import csv
import os
import re
from typing import Dict, List, Optional

import pandas as pd

try:
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
    from Bio.Seq import Seq
except ImportError:
    SeqIO = Seq = PairwiseAligner = None

from ng_amr.core.audit import (
    logger,
)
from ng_amr.core.constants import (
    AA_COMPARISON_FOLDER_NAME,
    DEFAULT_GENES,
    GENE_AA_PRODUCT_SYNONYMS,
    GENE_CONFIG,
    ORIGINAL_AA_FOLDER,
)
from ng_amr.core.ids import (
    assembly_id,
    canonical,
)
from ng_amr.genomics.alignment import (
    A4,
    REPORTLAB_OK,
    HexColor,
    _gene_token_in_text,
    black,
    mm,
    rl_canvas,
    white,
)
from ng_amr.reports.pdf import (
    _PDF_CURATED,
    _PDF_DISRUPT,
    _PDF_HEADER,
    _PDF_MARGIN,
    _PDF_MUTED,
    _PDF_NONCUR,
    _pdf_draw_kv_table,
    _pdf_write_seq,
    _PDFPage,
)


def _aa_clean_seq(seq: str) -> str:
    if not seq:
        return ""
    s = re.sub(r"[^A-Za-z*]", "", str(seq).upper())
    return s.rstrip("*")


def _aa_isolate_keys(name):
    return [assembly_id(name)]


def _gene_from_aa_header(description: str, rid: str = "") -> Optional[str]:
    text = f"{rid} {description}".lower()
    text = text.replace("_", " ")
    if "subunit beta'" in text or "subunit beta’" in text:
        text_for_rpob = text.replace("subunit beta'", "SUBUNIT_BETA_PRIME").replace(
            "subunit beta’", "SUBUNIT_BETA_PRIME"
        )
    else:
        text_for_rpob = text
    best_gene = None
    best_len = 0
    for gene, syns in GENE_AA_PRODUCT_SYNONYMS.items():
        probe = text_for_rpob if gene == "rpoB" else text
        for syn in syns:
            s = (syn or "").lower().strip()
            if not s:
                continue
            if s in probe and len(s) > best_len:
                best_gene = gene
                best_len = len(s)
    if best_gene:
        return best_gene
    known = list(dict.fromkeys(list(DEFAULT_GENES) + list(GENE_CONFIG.keys())))
    for g in known:
        if _gene_token_in_text(description, g) or _gene_token_in_text(rid, g):
            return g
    return None


def load_original_aa_index(aa_folder: str) -> Dict[str, Dict[str, str]]:
    index: Dict[str, Dict[str, str]] = {}
    if not aa_folder or not os.path.isdir(aa_folder):
        logger.warning("  [AA_Comparison] Original AA folder not found: %s", aa_folder)
        return index

    fasta_ext = (".fasta", ".fa", ".fna", ".faa", ".aa", ".txt")
    files: List[str] = []
    for root, _dirs, names in os.walk(aa_folder):
        for n in names:
            if n.lower().endswith(fasta_ext):
                files.append(os.path.join(root, n))
    logger.info("  [AA_Comparison] Scanning %d FASTA-like files under %s", len(files), aa_folder)

    known_genes = {g.lower(): g for g in list(DEFAULT_GENES) + list(GENE_CONFIG.keys())}
    n_mapped = 0

    sources = {}
    for path in sorted(files):
        parent = os.path.basename(os.path.dirname(path))
        stem = os.path.splitext(os.path.basename(path))[0]
        file_iso_keys = _aa_isolate_keys(stem)
        stem2 = re.sub(r"^\d+[\.\s_]+", "", stem)
        file_iso_keys = list(dict.fromkeys(file_iso_keys + _aa_isolate_keys(stem2)))

        gene_hint = None
        for g_l, g_orig in known_genes.items():
            if _gene_token_in_text(parent, g_orig) or _gene_token_in_text(stem, g_orig):
                gene_hint = g_orig
                break

        try:
            records = list(SeqIO.parse(path, "fasta"))
        except Exception as exc:
            logger.warning("  [AA_Comparison] Could not parse %s (%s)", path, exc)
            continue
        if not records:
            try:
                with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                    text = re.sub(r"\s+", "", fh.read())
                if text and re.fullmatch(r"[A-Za-z*]+", text):
                    records = [type("R", (), {"id": stem, "description": stem, "seq": text})()]
            except OSError:
                continue

        for rec in records:
            rid = getattr(rec, "id", "") or ""
            desc = getattr(rec, "description", "") or rid
            seq = _aa_clean_seq(str(getattr(rec, "seq", "")))
            if not seq:
                continue
            gene = _gene_from_aa_header(desc, rid) or gene_hint
            if not gene:
                continue
            gkey = gene.lower()
            bucket = index.setdefault(gkey, {})
            for k in file_iso_keys:
                pair = (gkey, k)
                if pair in sources:
                    logger.warning("AA collision %s: %s and %s", pair, sources[pair], path)
                    bucket[k] = None
                else:
                    bucket[k] = seq
                    sources[pair] = path
            n_mapped += 1

    total = sum(len(v) for v in index.values())
    logger.info(
        "  [AA_Comparison] Indexed %d gene(s), %d isolate-keyed AA sequences "
        "(%d header→gene mappings via product synonyms)",
        len(index), total, n_mapped,
    )
    return index


def _lookup_original_aa(index, gene, isolate):
    return index.get(gene.lower(), {}).get(canonical(isolate))


def compare_aa_pair(pipeline_aa: str, original_aa: str) -> dict:
    p = _aa_clean_seq(pipeline_aa)
    o = _aa_clean_seq(original_aa)
    result = {
        "pipeline_len": len(p),
        "original_len": len(o),
        "identical": False,
        "status": "MISSING_ORIGINAL",
        "first_mismatch": None,
        "n_mismatches": 0,
        "diff_summary": "",
        "pipeline_aa": p,
        "original_aa": o,
    }
    if not o:
        result["status"] = "MISSING_ORIGINAL"
        result["diff_summary"] = "Original AA sequence not found for this isolate/gene"
        return result
    if not p:
        result["status"] = "MISSING_PIPELINE"
        result["diff_summary"] = "Pipeline did not produce a translated AA sequence"
        return result
    if p == o:
        result["identical"] = True
        result["status"] = "IDENTICAL"
        result["diff_summary"] = "Sequences are identical (same-to-same)"
        return result
    n_mis = 0
    first = None
    limit = min(len(p), len(o))
    for i in range(limit):
        if p[i] != o[i]:
            n_mis += 1
            if first is None:
                first = i + 1
    if len(p) != len(o):
        n_mis += abs(len(p) - len(o))
        if first is None:
            first = limit + 1
    result["n_mismatches"] = n_mis
    result["first_mismatch"] = first
    result["status"] = "LENGTH_DIFF" if len(p) != len(o) else "SEQUENCE_DIFF"
    result["diff_summary"] = (
        f"NOT identical — pipeline_len={len(p)}, original_len={len(o)}, "
        f"mismatches≈{n_mis}, first_mismatch_pos={first}"
    )
    return result


def create_aa_comparison_gene_pdf(
    gene_name: str,
    rows: List[dict],
    output_folder: str,
) -> Optional[str]:
    if not REPORTLAB_OK:
        logger.info("    [AA_Comparison] reportlab not installed — gene PDF skipped for %s", gene_name)
        return None
    os.makedirs(output_folder, exist_ok=True)
    safe_gene = re.sub(r"[^A-Za-z0-9._-]+", "_", gene_name)
    pdf_path = os.path.join(output_folder, f"AA_Comparison_{safe_gene}.pdf")
    try:
        c = rl_canvas.Canvas(pdf_path, pagesize=A4)
        title = f"AA Comparison — {gene_name} (pipeline vs original)"
        page = _PDFPage(c, title=title)

        n_total = len(rows)
        n_id = sum(1 for r in rows if r.get("status") == "IDENTICAL")
        n_diff = sum(1 for r in rows if r.get("status") in ("SEQUENCE_DIFF", "LENGTH_DIFF"))
        n_miss_o = sum(1 for r in rows if r.get("status") == "MISSING_ORIGINAL")
        n_miss_p = sum(1 for r in rows if r.get("status") == "MISSING_PIPELINE")
        n_miss_gene = sum(1 for r in rows if r.get("status") == "GENE_MISSING")

        page.c.setFont("Helvetica-Bold", 13)
        page.c.setFillColor(_PDF_HEADER)
        page.c.drawString(_PDF_MARGIN, page.y - 13,
                          f"Neisseria gonorrhoeae  {gene_name}  — AA Identity QC")
        page.y -= 7 * mm
        page.text("Pipeline-translated AA  vs  Original AA files", size=10, color=HexColor("#1565c0"))
        page.text("Source: supplied amino-acid assembly directory", size=8, color=HexColor("#616161"))
        page.spacer(2)

        page.section("Summary")
        _pdf_draw_kv_table(page, [
            ("Gene", gene_name),
            ("Isolates evaluated", str(n_total)),
            ("IDENTICAL (same-to-same)", str(n_id)),
            ("SEQUENCE / LENGTH difference", str(n_diff)),
            ("Original AA missing", str(n_miss_o)),
            ("Pipeline AA missing", str(n_miss_p)),
            ("Gene not found in assembly", str(n_miss_gene)),
            ("Identity rate (of comparable pairs)",
             f"{(100.0 * n_id / max(1, n_id + n_diff)):.1f}%  ({n_id}/{n_id + n_diff})"),
        ])

        if n_diff:
            page.banner(f"{n_diff} isolate(s) show AA differences — review below",
                        bg=_PDF_DISRUPT, fg=white)
        elif n_id and not (n_miss_o or n_miss_p):
            page.banner("All comparable isolates are IDENTICAL (same-to-same)",
                        bg=_PDF_CURATED, fg=white)

        page.section("Per-isolate comparison")
        for r in rows:
            page.ensure(28)
            st = r.get("status", "")
            color = {
                "IDENTICAL": _PDF_CURATED,
                "SEQUENCE_DIFF": _PDF_DISRUPT,
                "LENGTH_DIFF": _PDF_DISRUPT,
                "MISSING_ORIGINAL": _PDF_NONCUR,
                "MISSING_PIPELINE": _PDF_NONCUR,
                "GENE_MISSING": _PDF_MUTED,
            }.get(st, black)
            page.text(f"Isolate: {r.get('isolate', '')}", size=10, bold=True, color=_PDF_HEADER)
            page.text(f"Status: {st}", size=9, bold=True, color=color)
            page.text(
                f"Pipeline len={r.get('pipeline_len', 0)}  |  "
                f"Original len={r.get('original_len', 0)}  |  "
                f"Mismatches≈{r.get('n_mismatches', 0)}  |  "
                f"First mismatch={r.get('first_mismatch')}",
                size=8, color=HexColor("#424242"),
            )
            page.text(r.get("diff_summary") or "", size=8, color=HexColor("#424242"))
            if st in ("SEQUENCE_DIFF", "LENGTH_DIFF") and r.get("first_mismatch"):
                pos = int(r["first_mismatch"]) - 1
                pseq = r.get("pipeline_aa") or ""
                oseq = r.get("original_aa") or ""
                lo = max(0, pos - 15)
                hi = pos + 16
                page.mono(f"  PIPE  [{lo+1}..]: {pseq[lo:hi]}", size=7, color=HexColor("#006400"))
                page.mono(f"  ORIG  [{lo+1}..]: {oseq[lo:hi]}", size=7, color=HexColor("#c62828"))
            page.hline(HexColor("#e0e0e0"), 0.4)
            page.spacer(1.5)

        diffs = [r for r in rows if r.get("status") in ("SEQUENCE_DIFF", "LENGTH_DIFF")]
        if diffs:
            page.new_page()
            page.section("Appendix — full AA sequences (non-identical only)")
            for r in diffs:
                page.ensure(20)
                page.text(f"{r.get('isolate')} — {r.get('status')}", size=9, bold=True)
                page.text("Pipeline AA:", size=8, color=HexColor("#1565c0"))
                _pdf_write_seq(page, r.get("pipeline_aa") or "", size=7, label="AA")
                page.spacer(1)
                page.text("Original AA:", size=8, color=HexColor("#c62828"))
                _pdf_write_seq(page, r.get("original_aa") or "", size=7, label="AA")
                page.spacer(2)

        c.save()
        logger.info("    [AA_Comparison] %s", os.path.basename(pdf_path))
        return pdf_path
    except Exception as e:
        logger.error("    [AA_Comparison] PDF failed for %s: %s", gene_name, e)
        return None


def create_aa_comparison_master_pdf(
    gene_summaries: List[dict],
    output_folder: str,
) -> Optional[str]:
    if not REPORTLAB_OK:
        return None
    os.makedirs(output_folder, exist_ok=True)
    pdf_path = os.path.join(output_folder, "AA_Comparison_MASTER_Summary.pdf")
    try:
        c = rl_canvas.Canvas(pdf_path, pagesize=A4)
        page = _PDFPage(c, title="AA Comparison — Master Summary")
        page.c.setFont("Helvetica-Bold", 13)
        page.c.setFillColor(_PDF_HEADER)
        page.c.drawString(_PDF_MARGIN, page.y - 13,
                          "NG AMR — Pipeline AA vs Original AA  (Master Identity QC)")
        page.y -= 7 * mm
        page.text("Source: supplied amino-acid assembly directory", size=8, color=HexColor("#616161"))
        page.spacer(2)
        page.section("Per-gene identity summary")
        table_rows = []
        for g in gene_summaries:
            table_rows.append((
                g.get("gene", ""),
                f"id={g.get('n_identical', 0)}  "
                f"diff={g.get('n_diff', 0)}  "
                f"miss_orig={g.get('n_missing_original', 0)}  "
                f"miss_pipe={g.get('n_missing_pipeline', 0)}  "
                f"n={g.get('n_total', 0)}",
            ))
        _pdf_draw_kv_table(page, table_rows or [("—", "no genes")])
        total_id = sum(g.get("n_identical", 0) for g in gene_summaries)
        total_diff = sum(g.get("n_diff", 0) for g in gene_summaries)
        page.spacer(3)
        page.section("Overall")
        _pdf_draw_kv_table(page, [
            ("Total IDENTICAL pairs", str(total_id)),
            ("Total DIFFERENT pairs", str(total_diff)),
            ("Overall identity rate",
             f"{(100.0 * total_id / max(1, total_id + total_diff)):.1f}%"),
        ])
        if total_diff == 0 and total_id > 0:
            page.banner("All comparable pipeline ↔ original AA pairs are IDENTICAL",
                        bg=_PDF_CURATED, fg=white)
        elif total_diff > 0:
            page.banner(f"{total_diff} pair(s) differ — see per-gene PDFs in AA_Comparison/",
                        bg=_PDF_DISRUPT, fg=white)
        c.save()
        logger.info("    [AA_Comparison] %s", os.path.basename(pdf_path))
        return pdf_path
    except Exception as e:
        logger.error("    [AA_Comparison] Master PDF failed: %s", e)
        return None


def generate_aa_comparison_package(
    gene_results: List[dict],
    output_root: str,
    original_aa_folder: Optional[str] = None,
) -> Optional[str]:
    aa_folder = original_aa_folder or ORIGINAL_AA_FOLDER
    out_dir = os.path.join(output_root, AA_COMPARISON_FOLDER_NAME)
    os.makedirs(out_dir, exist_ok=True)
    logger.info("  [AA_Comparison] Writing reports to %s", out_dir)

    index = load_original_aa_index(aa_folder)
    if not index:
        logger.warning("  [AA_Comparison] No original AA sequences indexed — "
                       "PDFs will mark all isolates as MISSING_ORIGINAL")

    gene_summaries = []
    csv_rows = []

    for gr in gene_results or []:
        gene = gr.get("gene") or "gene"
        rows = []
        for d in gr.get("all_data") or []:
            iso = d.get("isolate") or ""
            if not d.get("gene_found"):
                cmp = {
                    "isolate": iso, "status": "GENE_MISSING",
                    "pipeline_len": 0, "original_len": 0,
                    "n_mismatches": 0, "first_mismatch": None,
                    "diff_summary": "Gene not found in assembly — no pipeline AA",
                    "pipeline_aa": "", "original_aa": "",
                }
            else:
                pipe_aa = d.get("query_protein") or ""
                orig_aa = _lookup_original_aa(index, gene, iso) or ""
                cmp = compare_aa_pair(pipe_aa, orig_aa)
                cmp["isolate"] = iso
            rows.append(cmp)
            csv_rows.append({
                "Gene": gene,
                "Isolate": iso,
                "Status": cmp.get("status"),
                "Pipeline_AA_Length": cmp.get("pipeline_len"),
                "Original_AA_Length": cmp.get("original_len"),
                "Identical": cmp.get("identical", False),
                "N_Mismatches": cmp.get("n_mismatches"),
                "First_Mismatch_Pos": cmp.get("first_mismatch"),
                "Diff_Summary": cmp.get("diff_summary"),
            })

        create_aa_comparison_gene_pdf(gene, rows, out_dir)
        gene_summaries.append({
            "gene": gene,
            "n_total": len(rows),
            "n_identical": sum(1 for r in rows if r.get("status") == "IDENTICAL"),
            "n_diff": sum(1 for r in rows if r.get("status") in ("SEQUENCE_DIFF", "LENGTH_DIFF")),
            "n_missing_original": sum(1 for r in rows if r.get("status") == "MISSING_ORIGINAL"),
            "n_missing_pipeline": sum(1 for r in rows if r.get("status") == "MISSING_PIPELINE"),
        })

    create_aa_comparison_master_pdf(gene_summaries, out_dir)

    csv_path = os.path.join(out_dir, "AA_Comparison_All_Isolates.csv")
    try:
        if csv_rows and pd is not None:
            pd.DataFrame(csv_rows).to_csv(csv_path, index=False)
        elif csv_rows:
            with open(csv_path, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(csv_rows[0].keys()))
                w.writeheader()
                w.writerows(csv_rows)
        logger.info("  [AA_Comparison] CSV → %s", csv_path)
    except Exception as e:
        logger.warning("  [AA_Comparison] CSV write failed: %s", e)

    logger.info("  [AA_Comparison] Done — %d gene PDF(s) + master summary under %s",
                len(gene_summaries), out_dir)
    return out_dir
