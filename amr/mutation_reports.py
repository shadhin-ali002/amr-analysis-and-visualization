"""Mutation reports routines."""
from __future__ import annotations

from .inputs import (
    AD_MIC,
    ANTIBIOTIC_ORDER,
    ORIGINAL,
    QC,
    ROUNDED,
    defaultdict,
    openpyxl,
    re,
)

import json

import os

import tempfile as _tempfile

from pathlib import Path

from typing import Any, Dict, List, Optional, Tuple, Set

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Align import PairwiseAligner
except ImportError:
    SeqIO = Seq = PairwiseAligner = None

def safe_name(value):
    return re.sub(r'[^\w. -]', '_', str(value)).strip('. ') or 'unnamed'

try:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import inch, mm
    from reportlab.pdfgen import canvas as rl_canvas
    from reportlab.lib.colors import HexColor, black, white, Color
    REPORTLAB_OK = True
except ImportError:
    REPORTLAB_OK = False
    A4 = (595.2755905511812, 841.8897637795277)
    inch = 72.0
    mm = 2.834645669291339
    def HexColor(_x):
        return None
    black = white = Color = None
    rl_canvas = None

def translate_dna_to_protein(dna_seq: str) -> str:
    return str(Seq(dna_seq).translate(to_stop=False))

def _pairwise_global(a: str, b: str, open_gap: float = -10) -> Tuple[str, str]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(a, b)
    if not alignments:
        return a, b
    aln = alignments[0]
    a_out, b_out = [], []
    a_prev = b_prev = 0
    for (a_s, a_e), (b_s, b_e) in zip(aln.aligned[0], aln.aligned[1]):
        while a_prev < a_s:
            a_out.append(a[a_prev]); b_out.append("-"); a_prev += 1
        while b_prev < b_s:
            a_out.append("-"); b_out.append(b[b_prev]); b_prev += 1
        for x, y in zip(a[a_s:a_e], b[b_s:b_e]):
            a_out.append(x); b_out.append(y)
        a_prev, b_prev = a_e, b_e
    while a_prev < len(a):
        a_out.append(a[a_prev]); b_out.append("-"); a_prev += 1
    while b_prev < len(b):
        a_out.append("-"); b_out.append(b[b_prev]); b_prev += 1
    return "".join(a_out), "".join(b_out)

def align_protein_sequences(ref_prot: str, qry_prot: str,
                            open_gap: float = -5) -> Tuple[str, str]:
    return _pairwise_global(ref_prot, qry_prot, open_gap=open_gap)

def call_protein_mutations(ref_aln: str, qry_aln: str) -> List[dict]:
    mutations = []
    ref_pos = 0
    i = 0
    while i < len(ref_aln):
        r, q = ref_aln[i], qry_aln[i]
        if r != "-" and q != "-":
            ref_pos += 1
            if r != q:
                if q == "*":
                    mutations.append({
                        "type": "stop_gain", "position": ref_pos,
                        "ref_aa": r, "alt_aa": q, "notation": f"{r}{ref_pos}Ter",
                    })
                else:
                    mutations.append({
                        "type": "substitution", "position": ref_pos,
                        "ref_aa": r, "alt_aa": q, "notation": f"{r}{ref_pos}{q}",
                    })
            i += 1
        elif r == "-" and q != "-":
            inserted = []
            start_pos = ref_pos
            while i < len(ref_aln) and ref_aln[i] == "-" and qry_aln[i] != "-":
                inserted.append(qry_aln[i]); i += 1
            aa_seq = "".join(inserted)
            ins_pos = start_pos + 1
            notation = (f"ins{ins_pos}{aa_seq}" if len(aa_seq) == 1
                        else f"ins{ins_pos}_{ins_pos + len(aa_seq) - 1}{aa_seq}")
            mutations.append({
                "type": "insertion", "position": ins_pos, "ref_position": start_pos,
                "aa_seq": aa_seq, "notation": notation, "length": len(aa_seq),
            })
        elif r != "-" and q == "-":
            deleted = []
            start_pos = ref_pos
            while i < len(ref_aln) and ref_aln[i] != "-" and qry_aln[i] == "-":
                deleted.append(ref_aln[i]); ref_pos += 1; i += 1
            end_pos = ref_pos - 1
            aa_seq = "".join(deleted)
            notation = (f"del{start_pos}" if len(aa_seq) == 1
                        else f"del{start_pos}_{end_pos}")
            mutations.append({
                "type": "deletion", "position": start_pos, "end_position": end_pos,
                "aa_seq": aa_seq, "notation": notation, "length": len(aa_seq),
            })
        else:
            i += 1
    return mutations

def detect_frameshifts(protein_mutations: List[dict], nt_mutations: List[dict],
                       window: int = 12) -> List[dict]:
    nt_indels = [m for m in (nt_mutations or [])
                 if m.get("type") in ("insertion", "deletion")]
    if not nt_indels:
        return []
    nt_indels = sorted(nt_indels, key=lambda m: (m.get("position") or 0))

    stop_positions = sorted(
        m["position"] for m in (protein_mutations or [])
        if m.get("type") == "stop_gain" and m.get("position") is not None
    )

    frameshifts: List[dict] = []
    net = 0
    open_fs: Optional[dict] = None

    for nm in nt_indels:
        length = int(nm.get("length") or 0)
        if length <= 0:
            continue
        nt_pos = int(nm.get("position") or 1)
        net += length if nm["type"] == "insertion" else -length

        if net % 3 != 0:
            if open_fs is None:
                aa_pos = max(1, (nt_pos + 2) // 3)
                stop_pos = next((s for s in stop_positions if s > aa_pos), None)
                open_fs = {
                    "type": f"frameshift_{nm['type']}",
                    "position": aa_pos,
                    "nt_position": nt_pos,
                    "length": abs(net),
                    "aa_length": max(1, abs(net) // 3),
                    "net_shift": net,
                    "notation": f"fs{aa_pos}Ter{stop_pos}" if stop_pos else f"fs{aa_pos}",
                }
                if stop_pos:
                    open_fs["stop_position"] = stop_pos
                if nm["type"] == "deletion":
                    end_nt = int(nm.get("end_position") or nt_pos)
                    open_fs["end_position"] = max(open_fs["position"], (end_nt + 2) // 3)
                frameshifts.append(open_fs)
            else:
                open_fs["length"] = abs(net)
                open_fs["aa_length"] = max(1, abs(net) // 3)
                open_fs["net_shift"] = net
        else:
            if open_fs is not None:
                nt_open = int(open_fs.get("nt_position") or 0)
                if abs(nt_pos - nt_open) <= window and open_fs in frameshifts:
                    frameshifts.remove(open_fs)
            open_fs = None

    return frameshifts

def analyze_protein(gene_name: str, ref_dna: str, query_dna: str,
                    nt_mutations: Optional[List[dict]] = None,
                    open_gap: float = -5) -> dict:
    ref_prot = translate_dna_to_protein(ref_dna)
    qry_prot = translate_dna_to_protein(query_dna)
    ref_aln, qry_aln = align_protein_sequences(ref_prot, qry_prot, open_gap=open_gap)
    mutations = call_protein_mutations(ref_aln, qry_aln)
    frameshifts = detect_frameshifts(mutations, nt_mutations or [])
    silent: List[dict] = []
    all_muts = mutations + frameshifts + silent
    seen = set()
    unique = []
    for m in all_muts:
        if m["notation"] not in seen:
            seen.add(m["notation"])
            unique.append(m)
    return {
        "gene": gene_name,
        "ref_protein": ref_prot,
        "query_protein": qry_prot,
        "ref_alignment": ref_aln,
        "query_alignment": qry_aln,
        "mutations": unique,
        "substitutions": [m for m in unique if m["type"] == "substitution"],
        "insertions": [m for m in unique if m["type"] == "insertion"],
        "deletions": [m for m in unique if m["type"] == "deletion"],
        "stop_gains": [m for m in unique if m["type"] == "stop_gain"],
        "frameshifts": frameshifts,
        "silent_mutations": silent,
    }

def _inject_disruption_banner(html_path: str):
    if not os.path.exists(html_path):
        return
    with open(html_path, encoding="utf-8") as fh:
        html = fh.read()
    banner = ('<div style="background:#c62828;color:#fff;padding:10px;font-weight:bold;'
              'margin-bottom:12px;">⚠ DISRUPTED / LARGE STRUCTURAL CHANGE DETECTED</div>')
    if "DISRUPTED / LARGE" in html:
        return
    if "<body>" in html:
        html = html.replace("<body>", "<body>\n" + banner, 1)
    elif "<body " in html:
        import re as _re
        html = _re.sub(r"(<body[^>]*>)", r"\1\n" + banner, html, count=1)
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(html)

def generate_mutation_list_html(mutations, curated_set, silent_mutations, stop_gains=None, silent_notations=None):
    silent_set = set()
    silent_list = []
    if silent_mutations:
        for m in silent_mutations:
            silent_set.add(m['notation'])
            silent_list.append(m)
    if stop_gains:
        for m in stop_gains:
            silent_set.add(m['notation'])
    if silent_notations:
        for n in silent_notations:
            silent_set.add(n)
    curated = []
    non_curated = []
    for m in mutations:
        if m['notation'] in silent_set:
            continue
        if curated_set and m['notation'] in curated_set:
            curated.append(m)
        else:
            non_curated.append(m)
    html = '<div class="mutation-info">\n'
    html += '<h3>🔴 Amino Acid Mutations Found:</h3>\n'
    if curated:
        html += '<div class="curated-section">\n'
        html += '<div class="section-title curated">[CURATED] Mutations:</div>\n'
        for idx, m in enumerate(curated, 1):
            html += f'<div class="mutation-item">{idx}. {m["notation"]} <span class="true-badge">[TRUE]</span></div>\n'
        html += '</div>\n'
    if non_curated:
        html += '<div class="non-curated-section">\n'
        html += '<div class="section-title non-curated">[NON-CURATED] Mutations:</div>\n'
        for idx, m in enumerate(non_curated, 1):
            html += f'<div class="mutation-item">{idx}. {m["notation"]} <span class="false-badge">[FALSE]</span></div>\n'
        html += '</div>\n'
    if silent_list:
        html += '<div class="silent-section">\n'
        html += '<div class="section-title silent">[SILENT] Mutations:</div>\n'
        for idx, m in enumerate(silent_list, 1):
            html += f'<div class="mutation-item">{idx}. {m["notation"]} <span class="silent-badge">[SILENT]</span></div>\n'
        html += '</div>\n'
    html += '</div>\n'
    return html

_PDF_MARGIN = 1 * inch

_PDF_PAGE_W, _PDF_PAGE_H = A4

_PDF_PAGE_W, _PDF_PAGE_H = A4

_PDF_CONTENT_W = _PDF_PAGE_W - 2 * _PDF_MARGIN

_PDF_CURATED = HexColor("#2e7d32") if REPORTLAB_OK else None

_PDF_NONCUR  = HexColor("#e65100") if REPORTLAB_OK else None

_PDF_SILENT  = HexColor("#1565c0") if REPORTLAB_OK else None

_PDF_HEADER  = HexColor("#1565c0") if REPORTLAB_OK else None

_PDF_DISRUPT = HexColor("#c62828") if REPORTLAB_OK else None

_PDF_LIGHT   = HexColor("#f5f5f5") if REPORTLAB_OK else None

_PDF_SNP_REF = HexColor("#006400") if REPORTLAB_OK else None

_PDF_SNP_QRY = HexColor("#FF0000") if REPORTLAB_OK else None

_PDF_AA_SNP_QRY = HexColor("#FFD700") if REPORTLAB_OK else None

_PDF_INDEL   = HexColor("#800080") if REPORTLAB_OK else None

_PDF_POS     = HexColor("#9e9e9e") if REPORTLAB_OK else None

_PDF_MUTED   = HexColor("#757575") if REPORTLAB_OK else None

class _PDFPage:
    def __init__(self, c, title=""):
        self.c = c
        self.y = _PDF_PAGE_H - _PDF_MARGIN
        self.page = 1
        self.title = title

    def new_page(self):
        self.c.showPage()
        self.page += 1
        self.c.setFont("Helvetica", 8)
        self.c.setFillColor(_PDF_MUTED)
        self.c.drawString(_PDF_MARGIN, _PDF_PAGE_H - 0.45 * inch, self.title[:95])
        self.c.drawRightString(_PDF_PAGE_W - _PDF_MARGIN, _PDF_PAGE_H - 0.45 * inch, f"p. {self.page}")
        self.y = _PDF_PAGE_H - _PDF_MARGIN - 0.15 * inch

    def ensure(self, need_mm):
        if self.y - need_mm * mm < _PDF_MARGIN + 0.25 * inch:
            self.new_page()

    def hline(self, color=None, thick=0.6):
        if color is None:
            color = HexColor("#bdbdbd")
        self.c.setStrokeColor(color)
        self.c.setLineWidth(thick)
        self.c.line(_PDF_MARGIN, self.y, _PDF_PAGE_W - _PDF_MARGIN, self.y)
        self.y -= 2.5 * mm

    def text(self, s, size=12, bold=False, color=None, leading=None):
        if color is None:
            color = black
        if leading is None:
            leading = size * 1.35
        self.ensure(leading / mm + 1)
        font = "Helvetica-Bold" if bold else "Helvetica"
        self.c.setFont(font, size)
        self.c.setFillColor(color)
        max_w = _PDF_CONTENT_W
        words = str(s).split()
        line = ""
        for w in words:
            test = (line + " " + w).strip()
            if self.c.stringWidth(test, font, size) <= max_w:
                line = test
            else:
                if line:
                    self.c.drawString(_PDF_MARGIN, self.y - size, line)
                    self.y -= leading
                    self.ensure(leading / mm + 1)
                line = w
        if line:
            self.c.drawString(_PDF_MARGIN, self.y - size, line)
            self.y -= leading

    def mono(self, s, size=10, color=None, leading=None):
        if color is None:
            color = black
        if leading is None:
            leading = size * 1.4
        self.ensure(leading / mm + 1)
        self.c.setFont("Courier", size)
        self.c.setFillColor(color)
        self.c.drawString(_PDF_MARGIN, self.y - size, s)
        self.y -= leading

    def spacer(self, mm_h=3):
        self.y -= mm_h * mm

    def section(self, title, size=12):
        self.ensure(12)
        self.spacer(2)
        self.text(title, size=size, bold=True, color=_PDF_HEADER)
        self.hline(HexColor("#90caf9"), 1.0)

    def banner(self, msg, bg=None, fg=None, size=10):
        if bg is None:
            bg = _PDF_DISRUPT
        if fg is None:
            fg = white
        self.ensure(14)
        h = 8 * mm
        self.c.setFillColor(bg)
        self.c.roundRect(_PDF_MARGIN, self.y - h, _PDF_CONTENT_W, h, 3, fill=1, stroke=0)
        self.c.setFillColor(fg)
        self.c.setFont("Helvetica-Bold", size)
        self.c.drawCentredString(_PDF_PAGE_W / 2, self.y - h + 2.4 * mm, msg[:110])
        self.y -= h + 2.5 * mm

def _pdf_draw_side_by_side_summary(page, nt_rows, aa_rows):
    n = max(len(nt_rows), len(aa_rows))
    while len(nt_rows) < n:
        nt_rows.append(("", ""))
    while len(aa_rows) < n:
        aa_rows.append(("", ""))

    page.ensure(n * 5.5 + 14)
    hdr_h = 6 * mm
    half = _PDF_CONTENT_W / 2
    x0 = _PDF_MARGIN
    page.c.setFillColor(HexColor("#e3f2fd"))
    page.c.rect(x0, page.y - hdr_h, half, hdr_h, fill=1, stroke=0)
    page.c.setFillColor(HexColor("#fff8e1"))
    page.c.rect(x0 + half, page.y - hdr_h, half, hdr_h, fill=1, stroke=0)
    page.c.setStrokeColor(HexColor("#90a4ae"))
    page.c.setLineWidth(0.6)
    page.c.rect(x0, page.y - hdr_h, _PDF_CONTENT_W, hdr_h, fill=0, stroke=1)
    page.c.line(x0 + half, page.y, x0 + half, page.y - hdr_h)
    page.c.setFillColor(_PDF_HEADER)
    page.c.setFont("Helvetica-Bold", 10)
    page.c.drawCentredString(x0 + half / 2, page.y - hdr_h + 1.8 * mm, "Nucleotide Summary")
    page.c.drawCentredString(x0 + half + half / 2, page.y - hdr_h + 1.8 * mm, "Amino Acid Summary")
    page.y -= hdr_h

    row_h = 5.0 * mm
    lab_w = half * 0.58
    for i in range(n):
        page.ensure(6)
        bg = _PDF_LIGHT if i % 2 == 0 else white
        page.c.setFillColor(bg)
        page.c.rect(x0, page.y - row_h, _PDF_CONTENT_W, row_h, fill=1, stroke=0)
        page.c.setStrokeColor(HexColor("#e0e0e0"))
        page.c.setLineWidth(0.3)
        page.c.rect(x0, page.y - row_h, _PDF_CONTENT_W, row_h, fill=0, stroke=1)
        page.c.line(x0 + half, page.y, x0 + half, page.y - row_h)

        nlab, nval = nt_rows[i]
        alab, aval = aa_rows[i]
        page.c.setFillColor(black)
        page.c.setFont("Helvetica-Bold", 8)
        page.c.drawString(x0 + 1.2 * mm, page.y - row_h + 1.3 * mm, str(nlab)[:28])
        page.c.setFont("Helvetica", 8)
        page.c.drawString(x0 + lab_w + 1 * mm, page.y - row_h + 1.3 * mm, str(nval)[:18])
        page.c.setFont("Helvetica-Bold", 8)
        page.c.drawString(x0 + half + 1.2 * mm, page.y - row_h + 1.3 * mm, str(alab)[:28])
        page.c.setFont("Helvetica", 8)
        page.c.drawString(x0 + half + lab_w + 1 * mm, page.y - row_h + 1.3 * mm, str(aval)[:18])
        page.y -= row_h
    page.spacer(3)

def _pdf_mutation_blocks(page, mutations, curated_set, silent_notations=None):
    silent_set = set(silent_notations or [])
    curated, non_curated, silent = [], [], []
    for m in mutations:
        n = m.get("notation", "")
        if n in silent_set:
            silent.append(m)
        elif curated_set and n in curated_set:
            curated.append(m)
        else:
            non_curated.append(m)

    def _block(title, items, color):
        if not items:
            return
        page.ensure(10)
        page.text(title, size=10, bold=True, color=color)
        for idx, m in enumerate(items, 1):
            page.mono(f"  {idx}. {m.get('notation', '')}", size=9, color=black)

    _block("[CURATED] Mutations", curated, _PDF_CURATED)
    _block("[NON-CURATED] Mutations", non_curated, _PDF_NONCUR)
    if silent:
        _block("[SILENT] Mutations", silent, _PDF_SILENT)
    if not (curated or non_curated or silent):
        page.text("No amino-acid mutations recorded.", size=10)

def _pdf_nt_mutation_blocks(page, nt_categorized, is_coding=True):
    syn = nt_categorized.get("syn_snps") or []
    nonsyn = nt_categorized.get("non_syn_snps") or []
    ins = nt_categorized.get("ins_list") or []
    dels = nt_categorized.get("del_list") or []

    def _fmt_snp(m):
        return f"{m['position']}: {m['ref_nt']} → {m['alt_nt']}"

    def _block(title, items, color):
        if not items:
            return
        page.ensure(10)
        page.text(f"{title} ({len(items)})", size=10, bold=True, color=color)
        for idx, m in enumerate(items, 1):
            txt = _fmt_snp(m) if "ref_nt" in m else m.get("notation", "")
            page.mono(f"  {idx}. {txt}", size=9)

    if is_coding:
        _block("Synonymous SNPs", syn, _PDF_CURATED)
        _block("Non-synonymous SNPs", nonsyn, _PDF_NONCUR)
    else:
        _block("SNPs (noncoding)", nonsyn or syn, _PDF_NONCUR)
    _block("Insertions", ins, _PDF_SILENT)
    _block("Deletions", dels, HexColor("#6a1b9a"))
    if not (syn or nonsyn or ins or dels):
        page.text("No nucleotide mutations found.", size=10)

def _pdf_build_aa_mut_index(ref_aln, mutations):
    ref_pos_map = []
    rc = 0
    for aa in ref_aln:
        if aa != "-":
            rc += 1
            ref_pos_map.append(rc)
        else:
            ref_pos_map.append(None)
    idx = {}
    for m in mutations or []:
        t = m.get("type", "")
        if t in ("substitution", "stop_gain", "silent_mutation"):
            pos = m.get("position")
            for i, p in enumerate(ref_pos_map):
                if p == pos:
                    idx[i] = "snp"
                    break
        elif t in ("insertion", "frameshift_insertion"):
            pos = m.get("position")
            for i in range(len(ref_aln)):
                if ref_aln[i] == "-":
                    prev = None
                    for j in range(i - 1, -1, -1):
                        if ref_pos_map[j] is not None:
                            prev = ref_pos_map[j]
                            break
                    if (prev or 0) + 1 == pos:
                        idx[i] = "indel"
        elif t in ("deletion", "frameshift_deletion"):
            pos = m.get("position")
            end = m.get("end_position", pos)
            for i, p in enumerate(ref_pos_map):
                if p is not None and pos <= p <= end:
                    idx[i] = "indel"
    return idx, ref_pos_map

def _pdf_build_nt_mut_index(ref_aln, qry_aln, nt_mutations=None):
    idx = {}
    for i in range(len(ref_aln)):
        r, q = ref_aln[i], qry_aln[i]
        if r == "-" or q == "-":
            if r != q:
                idx[i] = "indel"
        elif r != q:
            idx[i] = "snp"
    if nt_mutations:
        ref_pos_map = []
        rc = 0
        for nt in ref_aln:
            if nt != "-":
                rc += 1
                ref_pos_map.append(rc)
            else:
                ref_pos_map.append(None)
        for m in nt_mutations:
            t = m.get("type", "")
            if t == "SNP":
                pos = m.get("position")
                for i, p in enumerate(ref_pos_map):
                    if p == pos:
                        idx[i] = "snp"
                        break
            elif t in ("insertion", "deletion"):
                pos = m.get("position")
                end = m.get("end_position", pos)
                for i, p in enumerate(ref_pos_map):
                    if p is not None and pos <= p <= (end or pos):
                        idx[i] = "indel"
    return idx

def _pdf_draw_colored_chars(c, x, y, chars, styles, size, char_w):
    for i, ch in enumerate(chars):
        st = styles[i] if i < len(styles) else None
        cx = x + i * char_w
        if st:
            if st in ("snp_ref", "aa_snp_ref"):
                bg, fg = _PDF_SNP_REF, white
            elif st == "snp_qry":
                bg, fg = _PDF_SNP_QRY, white
            elif st == "aa_snp_qry":
                bg, fg = _PDF_AA_SNP_QRY, black
            elif st in ("indel_ref", "indel_qry"):
                bg, fg = _PDF_INDEL, white
            else:
                bg, fg = None, black
            if bg is not None:
                c.setFillColor(bg)
                c.rect(cx - 0.3, y - 1.2, char_w + 0.2, size + 2.2, fill=1, stroke=0)
                c.setFillColor(fg)
            else:
                c.setFillColor(black)
        else:
            c.setFillColor(black)
        c.setFont("Courier-Bold" if st else "Courier", size)
        c.drawString(cx, y, ch)

def _pdf_seq_alignment_separator(page):
    page.spacer(2.0)
    page.ensure(10)
    stars = "*" * 48
    page.c.setFont("Courier-Bold", 10)
    page.c.setFillColor(black)
    page.c.drawCentredString(_PDF_PAGE_W / 2, page.y - 10, stars)
    page.y -= 12
    page.c.setFont("Courier-Bold", 10)
    page.c.setFillColor(black)
    page.c.drawCentredString(_PDF_PAGE_W / 2, page.y - 10, stars)
    page.y -= 14
    page.spacer(1.5)

def _pdf_combined_alignment_chunks(page, ref_prot_aln, qry_prot_aln,
                                   ref_nt_aln, qry_nt_aln,
                                   aa_mutations, nt_mutations, sample_name,
                                   codons_per_line=None):
    aa_mut_idx, aa_pos_map = _pdf_build_aa_mut_index(ref_prot_aln, aa_mutations)
    nt_mut_idx = _pdf_build_nt_mut_index(ref_nt_aln, qry_nt_aln, nt_mutations)

    nt_pos_map = []
    rc = 0
    for nt in ref_nt_aln:
        if nt != "-":
            rc += 1
            nt_pos_map.append(rc)
        else:
            nt_pos_map.append(None)

    mono_size = 10
    pos_size = 9
    char_w = page.c.stringWidth("M", "Courier", mono_size)
    if codons_per_line is None:
        codons_per_line = max(10, int(_PDF_CONTENT_W / (4 * char_w)))

    total_nt = len(ref_nt_aln)
    total_aa = len(ref_prot_aln)
    x0 = _PDF_MARGIN

    def _draw_aa_pos_ruler(n_aa, aa_start):
        page.ensure(8)
        pos_y = page.y - pos_size
        page.c.setFont("Courier", pos_size)
        page.c.setFillColor(_PDF_POS)
        for i in range(n_aa):
            pos = aa_pos_map[aa_start + i] if (aa_start + i) < len(aa_pos_map) else None
            if pos is None:
                continue
            label = str(pos)
            slot_cx = x0 + (i * 4 + 1.5) * char_w
            lw = page.c.stringWidth(label, "Courier", pos_size)
            page.c.drawString(slot_cx - lw / 2.0, pos_y, label)
        page.y -= pos_size * 1.4

    def _draw_nt_pos_ruler(n_nt, nt_start):
        page.ensure(8)
        pos_y = page.y - pos_size
        page.c.setFont("Courier", pos_size)
        page.c.setFillColor(_PDF_POS)
        for ci in range(0, n_nt, 3):
            gi = nt_start + ci
            pos = nt_pos_map[gi] if gi < len(nt_pos_map) else None
            if pos is None:
                for k in range(3):
                    if ci + k < n_nt and (nt_start + ci + k) < len(nt_pos_map):
                        p = nt_pos_map[nt_start + ci + k]
                        if p is not None:
                            pos = p
                            break
            if pos is None:
                continue
            label = str(pos)
            codon_idx = ci // 3
            slot_cx = x0 + (codon_idx * 4 + 1.5) * char_w
            lw = page.c.stringWidth(label, "Courier", pos_size)
            page.c.drawString(slot_cx - lw / 2.0, pos_y, label)
        page.y -= pos_size * 1.35

    for nt_start in range(0, total_nt, codons_per_line * 3):
        nt_end = min(nt_start + codons_per_line * 3, total_nt)
        aa_start = nt_start // 3
        aa_end = min(aa_start + codons_per_line, total_aa)
        n_aa = aa_end - aa_start
        n_nt = nt_end - nt_start

        page.ensure(54)
        page.c.setStrokeColor(HexColor("#cccccc"))
        page.c.setLineWidth(0.7)
        page.c.line(x0, page.y, _PDF_PAGE_W - _PDF_MARGIN, page.y)
        page.y -= 2.5 * mm

        def _aa_line(aln, mut_idx, is_ref):
            chars, styles = [], []
            for i in range(n_aa):
                gi = aa_start + i
                ch = aln[gi] if gi < len(aln) else "-"
                mt = mut_idx.get(gi)
                if mt == "snp" and gi < len(qry_prot_aln) and ref_prot_aln[gi] == qry_prot_aln[gi]:
                    mt = None
                if mt == "snp":
                    st = "aa_snp_ref" if is_ref else "aa_snp_qry"
                elif mt == "indel" or ch == "-":
                    st = "indel_ref" if is_ref else "indel_qry"
                elif (not is_ref) and ch == "*":
                    st = "aa_snp_qry"
                else:
                    st = None
                chars.extend([" ", ch, " ", " "])
                styles.extend([None, st, None, None])
            return chars, styles

        def _nt_line(aln, mut_idx, is_ref):
            chars, styles = [], []
            for i in range(n_nt):
                gi = nt_start + i
                ch = aln[gi] if gi < len(aln) else "-"
                mt = mut_idx.get(gi)
                if mt == "snp":
                    st = "snp_ref" if is_ref else "snp_qry"
                elif mt == "indel" or ch == "-":
                    st = ("indel_ref" if is_ref else "indel_qry") if ch != "-" else None
                else:
                    st = None
                chars.append(ch)
                styles.append(st)
                if (i + 1) % 3 == 0 and i + 1 < n_nt:
                    chars.append(" ")
                    styles.append(None)
            return chars, styles

        _draw_aa_pos_ruler(n_aa, aa_start)

        page.ensure(10)
        ch, st = _aa_line(ref_prot_aln, aa_mut_idx, True)
        _pdf_draw_colored_chars(page.c, x0, page.y - mono_size, ch, st, mono_size, char_w)
        page.y -= mono_size * 1.3

        _draw_nt_pos_ruler(n_nt, nt_start)

        page.ensure(10)
        ch, st = _nt_line(ref_nt_aln, nt_mut_idx, True)
        _pdf_draw_colored_chars(page.c, x0, page.y - mono_size, ch, st, mono_size, char_w)
        page.y -= mono_size * 1.25

        page.ensure(14)
        label_w = max(
            page.c.stringWidth("REF:", "Helvetica", 9),
            page.c.stringWidth(f"QRY: {sample_name}", "Helvetica", 9),
        ) + 8
        page.c.setFillColor(HexColor("#ee3f2fd"))
        page.c.roundRect(x0 - 2, page.y - 22, label_w, 20, 2, fill=1, stroke=0)
        page.c.setFont("Helvetica", 9)
        page.c.setFillColor(HexColor("#1565c0"))
        page.c.drawString(x0, page.y - 9, "REF:")
        page.y -= 11
        page.c.setFont("Helvetica", 9)
        page.c.setFillColor(HexColor("#1565c0"))
        page.c.drawString(x0, page.y - 9, f"QRY: {sample_name}")
        page.y -= 13

        page.ensure(10)
        ch, st = _nt_line(qry_nt_aln, nt_mut_idx, False)
        _pdf_draw_colored_chars(page.c, x0, page.y - mono_size, ch, st, mono_size, char_w)
        page.y -= mono_size * 1.3

        _draw_nt_pos_ruler(n_nt, nt_start)

        page.ensure(10)
        ch, st = _aa_line(qry_prot_aln, aa_mut_idx, False)
        _pdf_draw_colored_chars(page.c, x0, page.y - mono_size, ch, st, mono_size, char_w)
        page.y -= mono_size * 1.3

        _draw_aa_pos_ruler(n_aa, aa_start)

        page.y -= 1.5 * mm
        page.c.setStrokeColor(HexColor("#cccccc"))
        page.c.setLineWidth(0.7)
        page.c.line(x0, page.y, _PDF_PAGE_W - _PDF_MARGIN, page.y)
        page.y -= 3.0 * mm

def _pdf_legend(page, kind="aa"):
    page.ensure(12)
    page.text("Legend:", size=9, bold=True)
    items = []
    if kind in ("aa", "both"):
        items += [
            (_PDF_SNP_REF, white, "AA Ref change"),
            (_PDF_AA_SNP_QRY, black, "AA Qry change"),
            (_PDF_INDEL, white, "AA Indel"),
        ]
    if kind in ("nt", "both"):
        items += [
            (_PDF_SNP_REF, white, "NT SNP (Ref)"),
            (_PDF_SNP_QRY, white, "NT SNP (Qry)"),
            (_PDF_INDEL, white, "NT Indel"),
        ]
    x = _PDF_MARGIN
    y = page.y - 8
    for bg, fg, lab in items:
        page.c.setFillColor(bg)
        page.c.roundRect(x, y - 1, 11, 9, 1.5, fill=1, stroke=0)
        page.c.setFillColor(fg)
        page.c.setFont("Helvetica-Bold", 7)
        page.c.drawCentredString(x + 5.5, y + 1, "A")
        page.c.setFillColor(black)
        page.c.setFont("Helvetica", 8)
        page.c.drawString(x + 13, y, lab)
        x += page.c.stringWidth(lab, "Helvetica", 8) + 22
        if x > _PDF_PAGE_W - _PDF_MARGIN - 60:
            x = _PDF_MARGIN
            y -= 11
            page.y -= 11
    page.y = y - 6

def _pdf_write_seq(page, seq: str, size: int = 8, color=None, chunk: int = None,
                  label: str = ""):
    if color is None:
        color = HexColor("#212121")
    s = seq or ""
    if not s:
        s = "(empty)"

    page.c.setFont("Courier", size)
    char_w = page.c.stringWidth("M", "Courier", size)
    if char_w <= 0:
        char_w = size * 0.6
    if chunk is None or chunk <= 0:
        chunk = max(10, int(_PDF_CONTENT_W / char_w))

    label_w = 0.0
    if label:
        page.ensure(size / mm + 3)
        lab = label if label.endswith(":") else (label + ":")
        page.c.setFont("Helvetica-Bold", size)
        label_w = page.c.stringWidth(lab + " ", "Helvetica-Bold", size)
        page.c.setFillColor(HexColor("#e3f2fd"))
        page.c.roundRect(
            _PDF_MARGIN - 1, page.y - size - 1.5,
            label_w + 2, size + 3, 1.5, fill=1, stroke=0,
        )
        page.c.setFillColor(HexColor("#0d47a1"))
        page.c.setFont("Helvetica-Bold", size)
        page.c.drawString(_PDF_MARGIN, page.y - size, lab)
        remaining_w = max(char_w * 10, _PDF_CONTENT_W - label_w)
        first_budget = max(10, int(remaining_w / char_w))
        first = s[:first_budget]
        rest = s[first_budget:]
        page.c.setFont("Courier", size)
        page.c.setFillColor(color)
        page.c.drawString(_PDF_MARGIN + label_w, page.y - size, first)
        page.y -= size * 1.4
        s = rest
        if not s:
            return

    for i in range(0, len(s), chunk):
        page.mono(s[i:i + chunk], size=size, color=color)

def _draw_combined_pdf(c, sample_name, ref_dna, query_dna, protein_result, nt_result, curated_set=None,
                       is_disrupted=False, gene_name="gene", ref_name="", conceptual=False):
    title = f"{gene_name} Combined Alignment — {sample_name}"
    page = _PDFPage(c, title=title)

    try:
        nt_cat = categorize_nt_mutations(nt_result["mutations"], ref_dna, query_dna)
    except Exception:
        nt_cat = {
            "syn_count": 0, "non_syn_count": 0,
            "snp_count": nt_result.get("snp_count", 0),
            "ins_count": nt_result.get("insertion_count", 0),
            "del_count": nt_result.get("deletion_count", 0),
            "syn_snps": [], "non_syn_snps": [], "ins_list": [], "del_list": [],
        }

    muts = protein_result["mutations"]
    silent = protein_result.get("silent_mutations") or []
    silent_set = {m["notation"] for m in silent}
    ref_prot = protein_result.get("ref_protein", "")
    qry_prot = protein_result.get("query_protein", "")

    page.c.setFont("Helvetica-Bold", 12)
    page.c.setFillColor(_PDF_HEADER)
    page.c.drawString(_PDF_MARGIN, page.y - 12,
                      f"Neisseria gonorrhoeae  {gene_name}  (NT + AA) Alignment Summary")
    page.y -= 6.5 * mm
    page.text(f"Reference: {ref_name or gene_name}", size=10, color=HexColor("#1565c0"))
    page.text(f"Isolate: {sample_name}", size=10, color=HexColor("#1565c0"))
    page.spacer(1.5)

    if is_disrupted:
        page.banner("GENE DISRUPTED – large deletion / frameshift / stop-gain detected")
    if conceptual:
        page.banner("Conceptual translation (non-coding region) – visualisation only",
                    bg=HexColor("#ff8f00"), fg=black)

    page.section("Alignment Summary")
    nt_rows = [
        ("Reference length (nt)", str(len(ref_dna))),
        ("Query length (nt)", str(len(query_dna))),
        ("Total SNPs", str(nt_cat.get("snp_count", nt_result.get("snp_count", 0)))),
        ("Synonymous (silent)", str(nt_cat.get("syn_count", 0))),
        ("Non-synonymous", str(nt_cat.get("non_syn_count", 0))),
        ("Indels", str(nt_cat.get("ins_count", 0) + nt_cat.get("del_count", 0))),
        ("Alignment length (gaps)", str(len(nt_result.get("ref_alignment") or ""))),
    ]
    aa_rows = [
        ("Reference AA length", str(len(ref_prot))),
        ("Query AA length", str(len(qry_prot))),
        ("Substitutions", str(len(protein_result["substitutions"]))),
        ("Insertions", str(len(protein_result["insertions"]))),
        ("Deletions", str(len(protein_result["deletions"]))),
        ("Frameshifts", str(len(protein_result["frameshifts"]))),
        ("Total non-silent muts",
         str(len([m for m in muts if m["notation"] not in silent_set]))),
    ]
    _pdf_draw_side_by_side_summary(page, nt_rows, aa_rows)

    page.section("Amino-Acid Mutations")
    _pdf_mutation_blocks(page, muts, curated_set or set(), silent_set)

    page.section("Nucleotide Mutations (brief)")
    _pdf_nt_mutation_blocks(page, nt_cat, is_coding=True)

    page.new_page()
    page.section("Sequences")
    page.text(f"Reference: {ref_name or gene_name}", size=10, color=HexColor("#1565c0"))
    page.text("Reference Sequence:", size=9, color=HexColor("#424242"))
    _pdf_write_seq(page, ref_prot, size=8, label="AA")
    page.spacer(1.0)
    _pdf_write_seq(page, ref_dna, size=8, label="NT")
    page.spacer(2.0)
    page.text(f"Isolate: {sample_name}", size=10, color=HexColor("#1565c0"))
    page.text("Query Sequence:", size=9, color=HexColor("#424242"))
    _pdf_write_seq(page, qry_prot, size=8, label="AA")
    page.spacer(1.0)
    _pdf_write_seq(page, query_dna, size=8, label="NT")

    _pdf_seq_alignment_separator(page)
    page.section("Combined Alignment (AA + NT)")
    _pdf_legend(page, "both")
    page.spacer(1)
    _pdf_combined_alignment_chunks(
        page,
        protein_result["ref_alignment"],
        protein_result["query_alignment"],
        nt_result["ref_alignment"],
        nt_result["query_alignment"],
        muts,
        nt_result.get("mutations"),
        sample_name,
        codons_per_line=19,
    )

def create_combined_html_alignment(sample_name, ref_dna, query_dna, protein_result, nt_result, output_folder, prefix, curated_set=None, is_disrupted=False, gene_name="gene", ref_name="", conceptual=False, output_path=None):
    ref_prot = protein_result['ref_protein']
    qry_prot = protein_result['query_protein']
    ref_prot_aln = protein_result['ref_alignment']
    qry_prot_aln = protein_result['query_alignment']
    mutations = protein_result['mutations']
    substitutions = protein_result['substitutions']
    insertions = protein_result['insertions']
    deletions = protein_result['deletions']
    stop_gains = protein_result['stop_gains']
    frameshifts = protein_result['frameshifts']
    silent_mutations = protein_result['silent_mutations']
    ref_nt_aln = nt_result['ref_alignment']
    qry_nt_aln = nt_result['query_alignment']
    ref_aa_pos_map = []
    ref_counter = 0
    for aa in ref_prot_aln:
        if aa != '-':
            ref_counter += 1
            ref_aa_pos_map.append(ref_counter)
        else:
            ref_aa_pos_map.append(None)
    ref_nt_pos_map = []
    ref_counter = 0
    for nt in ref_nt_aln:
        if nt != '-':
            ref_counter += 1
            ref_nt_pos_map.append(ref_counter)
        else:
            ref_nt_pos_map.append(None)
    mutation_positions_set = set()
    for m in mutations:
        if m['type'] in ['substitution', 'stop_gain', 'silent_mutation']:
            mutation_positions_set.add(m['position'])
        elif m['type'] in ['insertion', 'frameshift_insertion']:
            mutation_positions_set.add(m['position'] - 1)
        elif m['type'] in ['deletion', 'frameshift_deletion']:
            for pos in range(m['position'], m.get('end_position', m['position']) + 1):
                mutation_positions_set.add(pos)
    aa_mutation_indices = {}
    for m in mutations:
        if m['type'] in ['substitution', 'stop_gain']:
            pos = m['position']
            for idx, aa_pos in enumerate(ref_aa_pos_map):
                if aa_pos == pos:
                    aa_mutation_indices[idx] = 'snp'
                    break
    for m in mutations:
        if m['type'] in ['insertion', 'frameshift_insertion']:
            pos = m['position']
            for idx in range(len(ref_prot_aln)):
                if ref_prot_aln[idx] == '-' and qry_prot_aln[idx] != '-':
                    prev_pos = None
                    for j in range(idx-1, -1, -1):
                        if ref_aa_pos_map[j] is not None:
                            prev_pos = ref_aa_pos_map[j]
                            break
                    if prev_pos is None:
                        prev_pos = 0
                    insert_pos = prev_pos + 1
                    if insert_pos == pos:
                        aa_mutation_indices[idx] = 'indel'
                        break
    for m in mutations:
        if m['type'] in ['deletion', 'frameshift_deletion']:
            pos = m['position']
            for idx, aa_pos in enumerate(ref_aa_pos_map):
                if aa_pos == pos and qry_prot_aln[idx] == '-':
                    aa_mutation_indices[idx] = 'indel'
                    break
    nt_mutation_indices = {}
    for i in range(len(ref_nt_aln)):
        r = ref_nt_aln[i]
        q = qry_nt_aln[i]
        if r == '-' and q != '-':
            nt_mutation_indices[i] = 'indel'
        elif r != '-' and q == '-':
            nt_mutation_indices[i] = 'indel'
        elif r != q:
            nt_mutation_indices[i] = 'snp'
    for m in nt_result['mutations']:
        if m['type'] == 'SNP':
            for idx, pos in enumerate(ref_nt_pos_map):
                if pos == m['position']:
                    nt_mutation_indices[idx] = 'snp'
        elif m['type'] == 'insertion':
            ins_pos = m['position']
            for idx, pos in enumerate(ref_nt_pos_map):
                if pos == ins_pos - 1:
                    for j in range(idx + 1, min(idx + 1 + m['length'], len(ref_nt_aln))):
                        if ref_nt_aln[j] == '-':
                            nt_mutation_indices[j] = 'indel'
                    break
            for idx, nt in enumerate(qry_nt_aln):
                if nt != '-' and ref_nt_aln[idx] == '-':
                    if idx > 0 and ref_nt_pos_map[idx-1] == ins_pos - 1:
                        nt_mutation_indices[idx] = 'indel'
        elif m['type'] == 'deletion':
            for idx, pos in enumerate(ref_nt_pos_map):
                if pos and m['position'] <= pos <= m.get('end_position', m['position']):
                    nt_mutation_indices[idx] = 'indel'
    silent_mutation_set = set()
    for m in silent_mutations:
        silent_mutation_set.add(m['notation'])
    for m in stop_gains:
        silent_mutation_set.add(m['notation'])
    nt_categorized = categorize_nt_mutations(nt_result['mutations'], ref_dna, query_dna)
    nt_syn_count = nt_categorized['syn_count']
    nt_non_syn_count = nt_categorized['non_syn_count']
    nt_snp_count = nt_result['snp_count']
    nt_indel_count = nt_result['insertion_count'] + nt_result['deletion_count']
    nt_ref_len = len(ref_dna)
    nt_qry_len = len(query_dna)
    nt_aln_len = len(ref_nt_aln)
    aa_ref_len = len(ref_prot)
    aa_qry_len = len(qry_prot)
    aa_sub_count = len(substitutions)
    aa_ins_count = len(insertions)
    aa_del_count = len(deletions)
    aa_fs_count = len(frameshifts)
    aa_silent_count = len(silent_mutations)
    aa_total_muts = len([m for m in mutations if m['notation'] not in silent_mutation_set])
    disruption_html = ""
    if is_disrupted:
        disruption_html = '<div style="background:#ff4444; color:white; padding:15px; border-radius:8px; text-align:center; font-size:1.5em; margin-bottom:20px;">🛑 GENE DISRUPTED – Large deletion (≥15 aa), frequent deletions (≥3 events), frameshift, or stop-gain detected</div>'
    conceptual_html = ""
    if conceptual:
        conceptual_html = '<div style="background:#ffcc00; color:#333; padding:10px; border-radius:5px; margin-bottom:15px; border:1px solid #cc9900;">⚠️ <b>Conceptual translation</b> – this region is annotated as non‑coding; the amino acid alignment is generated for visualization purposes only and should not be interpreted as a real protein product.</div>'
    ref_info = f" (Reference: {ref_name})" if ref_name else ""
    title_suffix = " (conceptual)" if conceptual else ""
    ref_label = ref_name if ref_name else gene_name
    html = """<html>
<head>
<meta charset="UTF-8">
<title>{gene_name} Combined Alignment - {sample}</title>
<style>
body{{font-family:Arial; margin:20px; background:#f5f5f5;}}
.container{{background:white; padding:20px; border-radius:10px; box-shadow:0 2px 10px rgba(0,0,0,0.1);}}
h2{{color:#2c3e50; border-bottom:2px solid #3498db; padding-bottom:10px;}}
h3{{color:#2c3e50; margin-top:20px;}}
.summary-side-by-side{{display:flex; gap:0; margin:15px 0;}}
.summary-box{{flex:1; background:#f8f9fa; border:1px solid #dee2e6; padding:15px;}}
.summary-box:first-child{{border-right:4px solid #333; border-radius:5px 0 0 5px;}}
.summary-box:last-child{{border-left:none; border-radius:0 5px 5px 0;}}
.summary-box table{{border-collapse:collapse; width:100%; table-layout:fixed;}}
.summary-box td{{padding:8px 12px; border:1px solid #dee2e6;}}
.summary-box td.col-label-left{{width:22%; font-weight:bold; background:#e9ecef;}}
.summary-box td.col-label-right{{width:17%; font-weight:bold; background:#e9ecef;}}
.summary-box td.col-value{{width:3%; text-align:center;}}
.summary-box td.col-spacer{{width:1%; background:transparent; border:none;}}
.summary-box .snp-color{{color:#CC0000; font-weight:bold;}}
.summary-box .indel-color{{color:#800080; font-weight:bold;}}
.summary-box .aa-snp-color{{color:#B8860B; font-weight:bold;}}
.summary-box .aa-indel-color{{color:#800080; font-weight:bold;}}
.summary-box .silent-color{{color:#006400; font-weight:bold;}}
.pre-wrap{{font-size:12px; line-height:1.45; overflow-x:auto; font-family:'Courier New',Courier,monospace; background:#fafafa; padding:15px; border-radius:5px; border:1px solid #e0e0e0; white-space:pre; letter-spacing:0; word-spacing:0;}}
.snp-ref{{background:#006400; font-weight:bold; color:#FFFFFF; padding:0px; display:inline-block; min-width:1px;}}
.snp-qry{{background:#FF0000; font-weight:bold; color:#FFFFFF; padding:0px; display:inline-block; min-width:1px;}}
.indel-ref{{background:#800080; font-weight:bold; color:#FFFFFF; padding:0px; display:inline-block; min-width:1px;}}
.indel-qry{{background:#800080; font-weight:bold; color:#FFFFFF; padding:0px; display:inline-block; min-width:1px;}}
.aa-snp-ref{{background:#006400; font-weight:bold; color:#FFFFFF; padding:0px 2px; display:inline-block; min-width:1px; border-radius:3px;}}
.aa-snp-qry{{background:#FFD700; font-weight:bold; color:#000000; padding:0px 2px; display:inline-block; min-width:1px; border-radius:3px;}}
.aa-indel-ref{{background:#800080; font-weight:bold; color:#FFFFFF; padding:0px 2px; display:inline-block; min-width:1px; border-radius:3px;}}
.aa-indel-qry{{background:#800080; font-weight:bold; color:#FFFFFF; padding:0px 2px; display:inline-block; min-width:1px; border-radius:3px;}}
.ref{{background:#eeeeee; font-weight:bold; color:#2c3e50;}}
.position{{color:#444444; font-size:12px; font-family:'Courier New', monospace; white-space:pre; letter-spacing:0; font-variant-numeric:tabular-nums;}}
.position-mutation{{color:#CC0000; font-weight:bold; font-size:12px; font-family:'Courier New', monospace; white-space:pre; letter-spacing:0; font-variant-numeric:tabular-nums;}}
.query{{background:#f0f8ff; font-weight:bold;}}
.aa-ref{{background:#eeeeee; font-weight:bold; color:#2c3e50;}}
.aa-query{{background:#f0f8ff; font-weight:bold;}}
.legend{{display:flex; flex-wrap:wrap; gap:10px; align-items:center; font-size:14px; margin-bottom:20px; padding:10px; background:#f9f9f9; border:1px solid #ddd; border-radius:5px;}}
.legend span{{padding:2px 8px; border-radius:3px; white-space:nowrap;}}
.legend .legend-label{{font-weight:bold; margin-right:5px;}}
.mutation-info{{margin:20px 0; padding:15px; background:#fff3cd; border-left:4px solid #ffc107; border-radius:5px;}}
.mutation-info h3{{margin-top:0; color:#856404;}}
.mutation-item{{margin:5px 0; font-family:'Courier New', monospace;}}
.mutation-item .true-badge{{color:#2e7d32; font-weight:bold;}}
.mutation-item .false-badge{{color:#c62828; font-weight:bold;}}
.mutation-item .silent-badge{{color:#006400; font-weight:bold;}}
.curated-section{{background-color:#e8f5e9; padding:10px; margin:10px 0; border-radius:5px; border-left:4px solid #2e7d32;}}
.non-curated-section{{background-color:#fff3e0; padding:10px; margin:10px 0; border-radius:5px; border-left:4px solid #e65100;}}
.silent-section{{background-color:#e3f2fd; padding:10px; margin:10px 0; border-radius:5px; border-left:4px solid #1565c0;}}
.section-title{{font-weight:bold; font-size:16px; margin-bottom:10px;}}
.section-title.curated{{color:#2e7d32;}}
.section-title.non-curated{{color:#e65100;}}
.section-title.silent{{color:#1565c0;}}
.alignment-separator{{
    border-top: 2px solid #cccccc;
    margin: 5px 0;
    padding: 0;
}}
.side-by-side {{
    display: flex;
    justify-content: center;
    align-items: stretch;
    padding: 0;
    gap: 0;
    margin: 20px 0;
}}
.side-by-side .aa-col,
.side-by-side .nt-col {{
    flex: 0 0 35%;
    min-width: 280px;
    background: #fafafa;
    padding: 15px;
    border: 1px solid #e0e0e0;
    box-sizing: border-box;
}}
.side-by-side .aa-col {{
    border-right: 4px solid #333;
    border-radius: 8px 0 0 8px;
}}
.side-by-side .nt-col {{
    border-left: none;
    border-radius: 0 8px 8px 0;
}}
.nt-syn-snp {{ background-color: #e8f5e9; padding: 10px; margin: 10px 0; border-radius: 5px; border-left: 4px solid #2e7d32; }}
.nt-nonsyn-snp {{ background-color: #fff3e0; padding: 10px; margin: 10px 0; border-radius: 5px; border-left: 4px solid #e65100; }}
.nt-ins-section {{ background-color: #e3f2fd; padding: 10px; margin: 10px 0; border-radius: 5px; border-left: 4px solid #1565c0; }}
.nt-del-section {{ background-color: #f3e5f5; padding: 10px; margin: 10px 0; border-radius: 5px; border-left: 4px solid #6a1b9a; }}
</style>
</head>
<body>
<div class="container">
{disruption}
{conceptual}
<h2>🔬🦠 <i>Neisseria gonorrhoeae</i> {gene_name} (NT+AA) Alignment_Summary{ref_info}{title_suffix}</h2>
<p><b>Reference:</b> {ref_label}</p>
<p style="color:#424242; margin:2px 0 2px 0;"><b>Reference Sequence:</b></p>
<p style="color:#616161; margin:0; font-size:9px;"><b>AA:</b></p>
<pre style="font-size:12px; line-height:1.35; white-space:pre-wrap; word-break:break-all; background:#fafafa; padding:8px; border-radius:4px; border:1px solid #e0e0e0;">{ref_aa_display}</pre>
<p style="color:#616161; margin:0; font-size:9px;"><b>NT:</b></p>
<pre style="font-size:12px; line-height:1.35; white-space:pre-wrap; word-break:break-all; background:#fafafa; padding:8px; border-radius:4px; border:1px solid #e0e0e0;">{ref_nt_display}</pre>
<p><b>Isolate:</b> {sample}</p>
<p style="color:#424242; margin:2px 0 2px 0;"><b>Query Sequence:</b></p>
<p style="color:#616161; margin:0; font-size:9px;"><b>AA:</b></p>
<pre style="font-size:12px; line-height:1.35; white-space:pre-wrap; word-break:break-all; background:#fafafa; padding:8px; border-radius:4px; border:1px solid #e0e0e0;">{qry_aa_display}</pre>
<p style="color:#616161; margin:0; font-size:9px;"><b>NT:</b></p>
<pre style="font-size:12px; line-height:1.35; white-space:pre-wrap; word-break:break-all; background:#fafafa; padding:8px; border-radius:4px; border:1px solid #e0e0e0;">{qry_nt_display}</pre>
<div class="summary-side-by-side">
    <div class="summary-box">
        <h3 style="margin-top:0; color:#2c3e50;">Nucleotide Alignment Summary</h3>
        <table>
            <tr>
                <td class="col-label-left">Reference Length (nt)</td>
                <td class="col-value">{nt_ref_len}</td>
                <td class="col-spacer"></td>
                <td class="col-label-right">Query Length (nt)</td>
                <td class="col-value">{nt_qry_len}</td>
            </tr>
            <tr>
                <td class="col-label-left">Total SNPs</td>
                <td class="col-value" style="color:#CC0000; font-weight:bold;">{nt_snp_count}</td>
                <td class="col-spacer"></td>
                <td class="col-label-right">Synonymous (silent)</td>
                <td class="col-value" style="color:#2e7d32; font-weight:bold;">{nt_syn_count}</td>
            </tr>
            <tr>
                <td class="col-label-left">Non-synonymous</td>
                <td class="col-value" style="color:#e65100; font-weight:bold;">{nt_non_syn_count}</td>
                <td class="col-spacer"></td>
                <td class="col-label-right">Indels</td>
                <td class="col-value" style="color:#800080; font-weight:bold;">{nt_indel_count}</td>
            </tr>
            <tr>
                <td class="col-label-left">Alignment Length (with gaps)</td>
                <td class="col-value">{nt_aln_len}</td>
                <td class="col-spacer"></td>
                <td class="col-label-right"></td>
                <td class="col-value"></td>
            </tr>
        </table>
    </div>
    <div class="summary-box">
        <h3 style="margin-top:0; color:#2c3e50;">Amino Acid Alignment Summary</h3>
        <table>
            <tr>
                <td class="col-label-left">Reference AA Length</td>
                <td class="col-value">{aa_ref_len}</td>
                <td class="col-spacer"></td>
                <td class="col-label-right">Query AA Length</td>
                <td class="col-value">{aa_qry_len}</td>
            </tr>
            <tr>
                <td class="col-label-left"><span class="aa-snp-color">Substitutions</span></td>
                <td class="col-value aa-snp-color">{aa_sub_count}</td>
                <td class="col-spacer"></td>
                <td class="col-label-right"><span class="aa-indel-color">Insertions</span></td>
                <td class="col-value aa-indel-color">{aa_ins_count}</td>
            </tr>
            <tr>
                <td class="col-label-left"><span class="aa-indel-color">Deletions</span></td>
                <td class="col-value aa-indel-color">{aa_del_count}</td>
                <td class="col-spacer"></td>
                <td class="col-label-right"><span class="snp-color">Frameshifts</span></td>
                <td class="col-value snp-color">{aa_fs_count}</td>
            </tr>
            <tr>
                <td class="col-label-left"><span class="silent-color">Silent Mutations</span></td>
                <td class="col-value silent-color">{aa_silent_count}</td>
                <td class="col-spacer"></td>
                <td class="col-label-right">Total Mutations</td>
                <td class="col-value" style="font-weight:bold;">{aa_total_muts}</td>
            </tr>
        </table>
    </div>
</div>
<div class="legend">
    <span class="legend-label">Legend:</span>
    <span style="background:#006400; color:#FFFFFF;">AA Ref Change</span>
    <span style="background:#FFD700; color:#000000;">AA Qry Change</span>
    <span style="background:#800080; color:#FFFFFF;">AA Indel</span>
    <span style="background:#006400; color:#FFFFFF;">NT SNP (Ref)</span>
    <span style="background:#FF0000; color:#FFFFFF;">NT SNP (Qry)</span>
    <span style="background:#800080; color:#FFFFFF;">NT Indel</span>
</div>
"""
    html = html.format(
        disruption=disruption_html,
        conceptual=conceptual_html,
        gene_name=gene_name,
        ref_info=ref_info,
        title_suffix=title_suffix,
        ref_label=ref_label,
        sample=sample_name,
        ref_aa_display=ref_prot,
        ref_nt_display=ref_dna,
        qry_aa_display=qry_prot,
        qry_nt_display=query_dna,
        nt_ref_len=nt_ref_len,
        nt_qry_len=nt_qry_len,
        nt_snp_count=nt_snp_count,
        nt_syn_count=nt_syn_count,
        nt_non_syn_count=nt_non_syn_count,
        nt_indel_count=nt_indel_count,
        nt_aln_len=nt_aln_len,
        aa_ref_len=aa_ref_len,
        aa_qry_len=aa_qry_len,
        aa_sub_count=aa_sub_count,
        aa_ins_count=aa_ins_count,
        aa_del_count=aa_del_count,
        aa_fs_count=aa_fs_count,
        aa_silent_count=aa_silent_count,
        aa_total_muts=aa_total_muts
    )
    html += generate_side_by_side_mutation_summary(
        mutations,
        curated_set,
        silent_mutations,
        stop_gains,
        nt_categorized
    )
    html += f'<h3>🦠 <i>Neisseria gonorrhoeae</i> {gene_name} Alignment (AA + NT)</h3>\n'
    html += '<div class="pre-wrap">\n'
    CODONS_PER_LINE = 46
    total_nt_positions = len(ref_nt_aln)
    total_aa_positions = len(ref_prot_aln)
    for nt_chunk_start in range(0, total_nt_positions, CODONS_PER_LINE * 3):
        nt_chunk_end = min(nt_chunk_start + CODONS_PER_LINE * 3, total_nt_positions)
        chunk_ref_nt = ref_nt_aln[nt_chunk_start:nt_chunk_end]
        chunk_qry_nt = qry_nt_aln[nt_chunk_start:nt_chunk_end]
        aa_start = nt_chunk_start // 3
        aa_end = min(aa_start + CODONS_PER_LINE, total_aa_positions)
        chunk_ref_aa = ref_prot_aln[aa_start:aa_end]
        chunk_qry_aa = qry_prot_aln[aa_start:aa_end]
        html += '<div class="alignment-separator"></div>\n'
        pos_ruler_display = ""
        for i in range(aa_start, aa_end):
            pos_val = ref_aa_pos_map[i]
            if pos_val is None:
                pos_ruler_display += "    "
            else:
                lab = str(pos_val)
                if len(lab) == 1:
                    pos_ruler_display += f" {lab}  "
                elif len(lab) == 2:
                    pos_ruler_display += f" {lab} "
                elif len(lab) == 3:
                    pos_ruler_display += f"{lab} "
                else:
                    pos_ruler_display += lab[:4]
        html += f'<span class="position">{pos_ruler_display}</span>\n'
        ref_aa_display = ""
        for i, aa in enumerate(chunk_ref_aa):
            idx = aa_start + i
            if idx in aa_mutation_indices:
                if aa_mutation_indices[idx] == 'snp':
                    cell = f'<span style="background:#006400;font-weight:bold;color:#FFFFFF;">{aa}</span>'
                else:
                    cell = f'<span style="background:#800080;font-weight:bold;color:#FFFFFF;">{aa}</span>'
            else:
                cell = f'<span style="font-weight:bold;">{aa}</span>'
            ref_aa_display += f" {cell}  "
        html += ref_aa_display + '\n'
        ref_nt_display = ""
        for i in range(0, len(chunk_ref_nt), 3):
            codon = chunk_ref_nt[i:i+3]
            codon_start = nt_chunk_start + i
            for j, char in enumerate(codon):
                pos = codon_start + j
                if pos in nt_mutation_indices and char != '-':
                    if nt_mutation_indices[pos] == 'snp':
                        ref_nt_display += f'<span style="background:#006400;font-weight:bold;color:#FFFFFF;">{char}</span>'
                    else:
                        ref_nt_display += f'<span style="background:#800080;font-weight:bold;color:#FFFFFF;">{char}</span>'
                else:
                    ref_nt_display += f'<span style="font-weight:bold;">{char}</span>'
            ref_nt_display += " "
        html += ref_nt_display + '\n'
        html += '<span class="ref">REF:</span>\n'
        html += '<span class="query">QRY: ' + sample_name + '</span>\n'
        qry_nt_display = ""
        for i in range(0, len(chunk_qry_nt), 3):
            codon = chunk_qry_nt[i:i+3]
            codon_start = nt_chunk_start + i
            for j, char in enumerate(codon):
                pos = codon_start + j
                if pos in nt_mutation_indices and char != '-':
                    if nt_mutation_indices[pos] == 'snp':
                        qry_nt_display += f'<span style="background:#FF0000;font-weight:bold;color:#FFFFFF;">{char}</span>'
                    else:
                        qry_nt_display += f'<span style="background:#800080;font-weight:bold;color:#FFFFFF;">{char}</span>'
                else:
                    qry_nt_display += f'<span style="font-weight:bold;">{char}</span>'
            qry_nt_display += " "
        html += qry_nt_display + '\n'
        qry_aa_display = ""
        for i, aa in enumerate(chunk_qry_aa):
            idx = aa_start + i
            if idx in aa_mutation_indices:
                if aa_mutation_indices[idx] == 'snp':
                    cell = f'<span style="background:#FFD700;font-weight:bold;color:#000000;">{aa}</span>'
                else:
                    cell = f'<span style="background:#800080;font-weight:bold;color:#FFFFFF;">{aa}</span>'
            else:
                cell = f'<span style="font-weight:bold;">{aa}</span>'
            qry_aa_display += f" {cell}  "
        html += qry_aa_display + '\n\n'
        html += '<div class="alignment-separator"></div>\n\n'
    html += '</div>\n'
    html += f'<p><i>Generated by {gene_name} mutation analysis pipeline</i></p>\n'
    html += '</div>\n'
    html += '</body>\n</html>'
    html_filename = f"{prefix}_{sample_name}_combined.html"
    html_path = output_path or os.path.join(output_folder, html_filename)
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return html_path

def analyze_nucleotide(ref_dna, qry_dna):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(ref_dna, qry_dna)
    if not alignments:
        ref_aln = ref_dna
        qry_aln = qry_dna
    else:
        alignment = alignments[0]
        ref_aln_list = []
        qry_aln_list = []
        r_blocks = alignment.aligned[0]
        q_blocks = alignment.aligned[1]
        r_prev = 0
        q_prev = 0
        for (r_start, r_end), (q_start, q_end) in zip(r_blocks, q_blocks):
            while r_prev < r_start:
                ref_aln_list.append(ref_dna[r_prev])
                qry_aln_list.append("-")
                r_prev += 1
            while q_prev < q_start:
                ref_aln_list.append("-")
                qry_aln_list.append(qry_dna[q_prev])
                q_prev += 1
            for r, q in zip(ref_dna[r_start:r_end], qry_dna[q_start:q_end]):
                ref_aln_list.append(r)
                qry_aln_list.append(q)
            r_prev = r_end
            q_prev = q_end
        while r_prev < len(ref_dna):
            ref_aln_list.append(ref_dna[r_prev])
            qry_aln_list.append("-")
            r_prev += 1
        while q_prev < len(qry_dna):
            ref_aln_list.append("-")
            qry_aln_list.append(qry_dna[q_prev])
            q_prev += 1
        ref_aln = "".join(ref_aln_list)
        qry_aln = "".join(qry_aln_list)
    ref_pos_map = []
    ref_counter = 0
    for nt in ref_aln:
        if nt != '-':
            ref_counter += 1
            ref_pos_map.append(ref_counter)
        else:
            ref_pos_map.append(None)
    mutations = []
    i = 0
    while i < len(ref_aln):
        r = ref_aln[i]
        q = qry_aln[i]
        if r != '-' and q != '-':
            if r != q:
                pos = ref_pos_map[i]
                mutations.append({
                    'type': 'SNP',
                    'position': pos,
                    'ref_nt': r,
                    'alt_nt': q,
                    'notation': f"{r}{pos}{q}"
                })
            i += 1
        elif r == '-' and q != '-':
            inserted = []
            pos_before = None
            for j in range(i-1, -1, -1):
                if ref_pos_map[j] is not None:
                    pos_before = ref_pos_map[j]
                    break
            if pos_before is None:
                pos_before = 0
            while i < len(ref_aln) and ref_aln[i] == '-' and qry_aln[i] != '-':
                inserted.append(qry_aln[i])
                i += 1
            alt_seq = ''.join(inserted)
            ins_pos = pos_before + 1
            if len(alt_seq) == 1:
                notation = f"ins{ins_pos}{alt_seq}"
            else:
                notation = f"ins{ins_pos}_{ins_pos+len(alt_seq)-1}{alt_seq}"
            mutations.append({
                'type': 'insertion',
                'position': ins_pos,
                'ref_nt': '-',
                'alt_nt': alt_seq,
                'notation': notation,
                'length': len(alt_seq)
            })
        elif r != '-' and q == '-':
            deleted = []
            start_pos = ref_pos_map[i]
            while i < len(ref_aln) and ref_aln[i] != '-' and qry_aln[i] == '-':
                deleted.append(ref_aln[i])
                i += 1
            end_pos = ref_pos_map[i-1] if i-1 >= 0 else start_pos
            ref_seq = ''.join(deleted)
            if len(deleted) == 1:
                notation = f"del{start_pos}"
            else:
                notation = f"del{start_pos}_{end_pos}"
            mutations.append({
                'type': 'deletion',
                'position': start_pos,
                'end_position': end_pos,
                'ref_nt': ref_seq,
                'alt_nt': '-',
                'notation': notation,
                'length': len(deleted)
            })
        else:
            i += 1
    return {
        'ref_sequence': ref_dna,
        'query_sequence': qry_dna,
        'ref_alignment': ref_aln,
        'query_alignment': qry_aln,
        'mutations': mutations,
        'snp_count': sum(1 for m in mutations if m['type'] == 'SNP'),
        'insertion_count': sum(1 for m in mutations if m['type'] == 'insertion'),
        'deletion_count': sum(1 for m in mutations if m['type'] == 'deletion')
    }

def categorize_nt_mutations(nt_mutations, ref_dna, query_dna, is_coding=True):
    snp_list = []
    ins_list = []
    del_list = []
    ref_dna = (ref_dna or "").upper()
    query_dna = (query_dna or "").upper()

    subst_at: Dict[int, str] = {}
    deleted_ref_positions = set()
    shift_events: List[Tuple[int, int]] = []

    for m in nt_mutations or []:
        mtype = m.get("type")
        if mtype == "SNP":
            if m.get("position") is not None:
                subst_at[int(m["position"])] = str(m.get("alt_nt") or "")
        elif mtype == "insertion":
            shift_events.append((int(m.get("position") or 1),
                                 int(m.get("length") or 0)))
        elif mtype == "deletion":
            start = int(m.get("position") or 1)
            end = int(m.get("end_position") or start)
            for p in range(start, end + 1):
                deleted_ref_positions.add(p)
            shift_events.append((start, -int(m.get("length") or (end - start + 1))))

    shift_events.sort(key=lambda t: t[0])

    def _net_shift_before(ref_position: int) -> int:
        total = 0
        for pos, delta in shift_events:
            if pos <= ref_position:
                total += delta
            else:
                break
        return total

    def _query_base_for(ref_position: int) -> Optional[str]:
        if ref_position in deleted_ref_positions:
            return None
        if ref_position in subst_at:
            return subst_at[ref_position]
        if 1 <= ref_position <= len(ref_dna):
            return ref_dna[ref_position - 1]
        return None

    missense_count = 0
    nonsense_count = 0
    syn_count = 0
    undetermined_count = 0

    for m in nt_mutations or []:
        if m["type"] == "SNP":
            pos = m.get("position")
            if not is_coding:
                m["effect"] = "SNP (noncoding)"
                m["effect_class"] = "noncoding"
                undetermined_count += 1
                snp_list.append(m)
                continue
            if pos is None:
                m["effect"] = "SNP"
                m["effect_class"] = "undetermined"
                undetermined_count += 1
                snp_list.append(m)
                continue

            pos = int(pos)
            codon_index = (pos - 1) // 3
            codon_start = codon_index * 3 + 1
            codon_positions = [codon_start, codon_start + 1, codon_start + 2]

            if codon_positions[-1] > len(ref_dna):
                m["effect"] = "SNP (incomplete reference codon)"
                m["effect_class"] = "undetermined"
                undetermined_count += 1
                snp_list.append(m)
                continue

            if _net_shift_before(codon_start) % 3 != 0:
                m["effect"] = "SNP (downstream of frameshift)"
                m["effect_class"] = "undetermined"
                undetermined_count += 1
                snp_list.append(m)
                continue

            qry_bases = [_query_base_for(p) for p in codon_positions]
            if any(b is None or len(b) != 1 for b in qry_bases):
                m["effect"] = "SNP (indel-affected codon)"
                m["effect_class"] = "undetermined"
                undetermined_count += 1
                snp_list.append(m)
                continue

            ref_codon = ref_dna[codon_start - 1: codon_start + 2]
            qry_codon = "".join(qry_bases)
            ref_aa = translate_dna_to_protein(ref_codon)
            qry_aa = translate_dna_to_protein(qry_codon)
            m["ref_codon"] = ref_codon
            m["query_codon"] = qry_codon
            m["codon_number"] = codon_index + 1
            m["ref_aa"] = ref_aa
            m["query_aa"] = qry_aa
            if ref_aa == qry_aa:
                m["effect"] = "Synonymous (silent)"
                m["effect_class"] = "synonymous"
                syn_count += 1
            elif qry_aa == "*":
                m["effect"] = "Nonsense"
                m["effect_class"] = "nonsense"
                nonsense_count += 1
            else:
                m["effect"] = "Non-synonymous (missense)"
                m["effect_class"] = "missense"
                missense_count += 1
            snp_list.append(m)
        elif m["type"] == "insertion":
            ins_list.append(m)
        elif m["type"] == "deletion":
            del_list.append(m)

    syn_snps = [m for m in snp_list if m.get("effect_class") == "synonymous"]
    non_syn_snps = [m for m in snp_list if m.get("effect_class") != "synonymous"]
    return {
        "syn_snps": syn_snps,
        "non_syn_snps": non_syn_snps,
        "ins_list": ins_list,
        "del_list": del_list,
        "syn_count": syn_count,
        "non_syn_count": len(non_syn_snps),
        "missense_count": missense_count,
        "nonsense_count": nonsense_count,
        "undetermined_count": undetermined_count,
        "snp_count": len(snp_list),
        "ins_count": len(ins_list),
        "del_count": len(del_list),
    }

def generate_side_by_side_mutation_summary(protein_mutations, curated_set,
                                           silent_mutations, stop_gains,
                                           nt_categorized):
    aa_html = generate_mutation_list_html(
        protein_mutations,
        curated_set,
        silent_mutations,
        stop_gains=stop_gains
    )
    syn_snps = nt_categorized['syn_snps']
    non_syn_snps = nt_categorized['non_syn_snps']
    ins_list = nt_categorized['ins_list']
    del_list = nt_categorized['del_list']
    def format_snp(m):
        return f"{m['position']}: {m['ref_nt']} → {m['alt_nt']}"
    def format_indel(m):
        return m['notation']
    nt_html = '<div class="nt-mutation-summary">\n'
    nt_html += '<h3>🔴 Nucleotide Mutations Found:</h3>\n'
    if syn_snps:
        nt_html += '<div class="nt-syn-snp">\n'
        nt_html += f'<div class="section-title">Synonymous SNPs ({len(syn_snps)})</div>\n'
        for idx, m in enumerate(syn_snps, 1):
            nt_html += f'<div class="mutation-item">{idx}. {format_snp(m)}</div>\n'
        nt_html += '</div>\n'
    if non_syn_snps:
        nt_html += '<div class="nt-nonsyn-snp">\n'
        nt_html += f'<div class="section-title">Non-synonymous SNPs ({len(non_syn_snps)})</div>\n'
        for idx, m in enumerate(non_syn_snps, 1):
            nt_html += f'<div class="mutation-item">{idx}. {format_snp(m)}</div>\n'
        nt_html += '</div>\n'
    if ins_list:
        nt_html += '<div class="nt-ins-section">\n'
        nt_html += f'<div class="section-title">Insertions ({len(ins_list)})</div>\n'
        for idx, m in enumerate(ins_list, 1):
            nt_html += f'<div class="mutation-item">{idx}. {format_indel(m)}</div>\n'
        nt_html += '</div>\n'
    if del_list:
        nt_html += '<div class="nt-del-section">\n'
        nt_html += f'<div class="section-title">Deletions ({len(del_list)})</div>\n'
        for idx, m in enumerate(del_list, 1):
            nt_html += f'<div class="mutation-item">{idx}. {format_indel(m)}</div>\n'
        nt_html += '</div>\n'
    if not (syn_snps or non_syn_snps or ins_list or del_list):
        nt_html += '<div class="mutation-item">No nucleotide mutations found.</div>\n'
    nt_html += '</div>\n'
    combined = """
    <div class="side-by-side">
        <div class="aa-col">
            {aa_html}
        </div>
        <div class="nt-col">
            {nt_html}
        </div>
    </div>
    """.format(aa_html=aa_html, nt_html=nt_html)
    return combined

CODING = {"penA", "ponA1", "porB1b", "gyrA", "parC", "parE", "mtrR", "mtrC", "folP", "rpoB", "rpoD",
              "rplD", "rplV", "rpsE", "rpsJ", "blaTEM", "tetM", "ermA", "ermB", "ermC", "ermF", "ereA", "ereB", "mefA"}

_CACHE = {"hits": {}, "contigs": {}, "refs": None}

def _inject_calls_banner(html_path, disp, gene, r):
    if not html_path or not os.path.exists(html_path):
        return
    with open(html_path, encoding="utf-8") as fh:
        html = fh.read()
    if "VALIDATED CALLS" in html:
        return
    calls = r.get("mutations", [])
    items = "".join(
        f'<li><b>{m["mutation"]}</b> &mdash; {m["type"]} &mdash; '
        f'<span style="color:{"#b45309" if m["curated"] else "#374151"}">{"Curated" if m["curated"] else "Non-Curated"}</span></li>'
        for m in calls) or "<li>No mutations relative to the reference</li>"
    banner = ('<div style="border:2px solid #1f3864;background:#eef2fb;padding:10px;margin-bottom:12px;'
              'font-family:Arial,sans-serif;font-size:13px;">'
              f'<div style="font-weight:bold;color:#1f3864;">VALIDATED CALLS &mdash; {disp} &mdash; {gene} '
              f'(reference: {r.get("ref_name")}; status: {r.get("status")})</div>'
              f'<ul style="margin:6px 0 6px 18px;">{items}</ul>'
              '<div style="color:#555;">The alignment below is drawn by the alignment renderer from the extracted '
              'sequence. Labels inside the alignment use the renderer\'s own notation (raw reference index for promoters/rRNA, '
              'insertions numbered +1, residues after a frameshift not meaningful). Where they differ, the calls above are '
              'authoritative.</div></div>')
    if "<body>" in html:
        html = html.replace("<body>", "<body>\n" + banner, 1)
    else:
        html = re.sub(r"(<body[^>]*>)", lambda mm: mm.group(1) + "\n" + banner, html, count=1)
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(html)

def _prepare_report_section(gene, row, curated, path, include_html):
    """Compute the display alignment once and reuse it in both report groupings."""
    from html import escape
    if not row["gene_found"] or not row.get("query_dna"):
        return None
    ref = _CACHE["refs"][row["ref_name"]]
    query = row["query_dna"]
    nt = analyze_nucleotide(ref, query)
    aa = analyze_protein(gene, ref, query, nt_mutations=nt["mutations"],
                         open_gap=-10.0 if gene == "porB1b" else -5.0)
    conceptual = gene not in CODING and gene != "mtrD"
    flags = dict(is_disrupted=bool(row.get("is_disrupted")), gene_name=gene,
                 ref_name=row["ref_name"], conceptual=conceptual)
    if include_html:
        sample = escape(row["isolate"])
        target = path.with_suffix(".html")
        create_combined_html_alignment(sample, ref, query, aa, nt, str(path.parent), gene, curated,
                                       output_path=target, **flags)
        _inject_calls_banner(target, sample, gene, row)
        if flags["is_disrupted"]:
            _inject_disruption_banner(target)
    with path.open("w", encoding="utf-8") as handle:
        json.dump([ref, query, aa, nt, flags], handle)
    return path

def _group_report_stems(labels):
    stems, used = {}, set()
    for label in labels:
        stem = safe_name(label)[:100].rstrip(" .") or "unnamed"
        if re.fullmatch(r"(?i:CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])", stem.split(".")[0]):
            stem = "_" + stem
        base = stem
        index = 1
        while stem.casefold() in used:
            index += 1
            stem = f"{base}_{index}"
        used.add(stem.casefold())
        stems[label] = stem
    return stems

def _write_report_group(folder, stem, title, records, book, include_html):
    from contextlib import ExitStack
    from html import escape
    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter
    from textwrap import wrap
    if not REPORTLAB_OK:
        raise ImportError("Grouped PDF reports require reportlab")
    workbook = openpyxl.Workbook()
    summary = workbook.active
    summary.title = "Summary"
    summary.append(["Gene", "Isolate", "Cohort", "Call status", "Reference", "Identity (%)",
                    "Coverage (%)", "Called variants", "Alignment"])
    calls = workbook.create_sheet("Mutation calls")
    calls.append(["Gene", "Isolate", "Position", "Mutation", "Type", "Curated", "Detail"])
    align = workbook.create_sheet("NT and AA alignments")
    align.append(["Gene", "Isolate", "Molecule", "Alignment column (1-based)", "Reference", "Observed"])
    mic = workbook.create_sheet("MICs")
    mic.append(["Isolate", "Antibiotic", "Original EST MIC", "Rounded EST MIC", "AD MIC"])
    for iso in dict.fromkeys(row["isolate"] for _, row, _, _ in records):
        for drug in ANTIBIOTIC_ORDER:
            values = book.values.get((iso, drug), {})
            mic.append([iso, drug, *(values.get(field) for field in (ORIGINAL, ROUNDED, AD_MIC))])

    page_path = folder / (stem + ".pdf")
    canvas = rl_canvas.Canvas(str(page_path), pagesize=A4)
    canvas.setTitle(title)
    with ExitStack() as stack:
        html_out = None
        if include_html:
            html_out = stack.enter_context((folder / (stem + ".html")).open("w", encoding="utf-8"))
            html_out.write('<!doctype html><html lang="en"><meta charset="utf-8"><title>' + escape(title)
                           + '</title><style>body{font:16px Arial,sans-serif;margin:2rem;color:#24334a}'
                             'section{margin:2rem 0}iframe{display:block;width:100%;height:80vh;border:1px solid #ced5df}'
                             'nav a{display:inline-block;margin:.25em .65em .25em 0}</style><h1>'
                           + escape(title) + '</h1><nav aria-label="Report sections">')
            for i, (gene, row, _, _) in enumerate(records):
                html_out.write(f'<a href="#s{i}">{escape(gene)} / {escape(row["isolate"])}</a>')
            html_out.write("</nav>")

        for i, (gene, row, curated, data_path) in enumerate(records):
            iso = row["isolate"]
            status = str(row.get("status") or "Unreported")
            count = len(row["mutations"]) if row["gene_found"] else None
            align_status = "Available" if data_path else "Unavailable"
            summary.append([gene, iso, "QC" if iso == QC else "Clinical", status, row.get("ref_name"),
                            row.get("pident"), row.get("qcov"), count, align_status])
            heading = f"{gene} / {iso}"
            canvas.bookmarkPage(f"section{i}")
            canvas.addOutlineEntry(heading, f"section{i}", level=0)
            page = _PDFPage(canvas, title=heading)
            page.section(title)
            page.text(heading, size=14, bold=True)
            page.text(f"Call status: {status}", size=11)
            page.text(f"Reference: {row.get('ref_name') or 'Unavailable'}", size=11)
            page.section("Mutation caller results")
            for item in row["mutations"]:
                calls.append([gene, iso, item.get("position"), item["mutation"], item["type"],
                              "Yes" if item["curated"] else "No", item.get("detail")])
                description = f"{item['mutation']} | {item['type']} | {item.get('detail', '')}"
                for line in wrap(description, 88):
                    page.text(line, size=9)
            if not row["mutations"]:
                page.text("No variants called." if row["gene_found"] else "Locus unassessed.", size=11)
            if not data_path:
                page.text("Alignment unavailable.", size=11)
            canvas.showPage()
            if html_out:
                html_out.write(f'<section id="s{i}"><h2>{escape(heading)}</h2><p>{escape(status)}</p>')
            if data_path:
                with data_path.open(encoding="utf-8") as handle:
                    ref, query, aa, nt, flags = json.load(handle)
                _draw_combined_pdf(canvas, iso, ref, query, aa, nt, curated, **flags)
                canvas.showPage()
                for mol, result in (("NT", nt), ("AA (conceptual)" if flags["conceptual"] else "AA", aa)):
                    ref_aln, obs_aln = result["ref_alignment"], result["query_alignment"]
                    for start in range(0, len(ref_aln), 120):
                        align.append([gene, iso, mol, start + 1, ref_aln[start:start+120], obs_aln[start:start+120]])
                if html_out:
                    fragment = data_path.with_suffix(".html").read_text(encoding="utf-8")
                    html_out.write(f'<iframe loading="lazy" sandbox title="{escape(heading)}" '
                                   f'srcdoc="{escape(fragment, quote=True)}"></iframe>')
            elif html_out:
                html_out.write("<p>Alignment unavailable.</p>")
                for item in row["mutations"]:
                    html_out.write("<p>" + escape(str(item["mutation"])) + "</p>")
            if html_out:
                html_out.write("</section>")
        if html_out:
            html_out.write("</html>")
    canvas.save()
    for sheet, widths in zip(workbook, ([22, 22, 14, 32, 25, 16, 16, 18, 18], [22, 22, 14, 30, 22, 13, 68],
                                        [22, 22, 20, 26, 72, 72], [22, 25, 22, 22, 22])):
        sheet.freeze_panes = "C2" if sheet != mic else "B2"
        sheet.auto_filter.ref = sheet.dimensions
        for col, width in enumerate(widths, 1):
            sheet.column_dimensions[get_column_letter(col)].width = width
        for cell in sheet[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="1F4E79")
        for cells in sheet:
            for cell in cells:
                if cell.data_type == "f":
                    cell.data_type = "s"
    # Preserve the original presence/absence matrix alongside the two alignment views.
    variants = list(dict.fromkeys((gene, item["mutation"]) for gene, row, _, _ in records
                                  for item in row["mutations"]))
    for offset in range(0, max(len(variants), 1), 16000):
        name = "Mutation Matrix" if offset == 0 else f"Mutation Matrix {offset // 16000 + 1}"
        table = workbook.create_sheet(name)
        panel = variants[offset:offset + 16000]
        table.append(["Gene", "Isolate", "Call status", *(f"{g} | {m}" for g, m in panel)])
        for gene, row, _, _ in records:
            observed = {item["mutation"] for item in row["mutations"]}
            table.append([gene, row["isolate"], row.get("status"), *(
                ("YES" if marker in observed else "NO") if g == gene and row["gene_found"] else None
                for g, marker in panel)])
        table.freeze_panes = "D2"
        table.auto_filter.ref = table.dimensions
        for cell in table[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="1F4E79")
            table.column_dimensions[cell.column_letter].width = min(45, max(16, len(str(cell.value)) + 3))
        for cells in table.iter_rows(min_row=2, max_col=3):
            for cell in cells:
                if cell.data_type == "f":
                    cell.data_type = "s"
    workbook.save(folder / (stem + ".xlsx"))
    workbook.close()

def grouped_mutation_reports(results, book, root, include_html=True):
    """Publish By_Gene and By_Isolate together once all report formats have succeeded."""
    target = Path(root) / "Mutation_Reports"
    if target.exists():
        raise FileExistsError(target)
    with _tempfile.TemporaryDirectory(dir=root, prefix=".mutation-reports-") as work:
        cache = Path(work)
        ready = cache / "Mutation_Reports"
        by_gene, by_isolate = defaultdict(list), defaultdict(list)
        index = 0
        for result in results:
            gene, curated = result["gene"], result["curated_set"]
            for row in result["all_data"]:
                section = _prepare_report_section(gene, row, curated, cache / f"{index}.json", include_html)
                index += 1
                record = (gene, row, curated, section)
                by_gene[gene].append(record)
                by_isolate[row["isolate"]].append(record)
        for name, groups, label in (("By_Gene", by_gene, "Gene"), ("By_Isolate", by_isolate, "Isolate")):
            folder = ready / name
            folder.mkdir(parents=True)
            stems = _group_report_stems(groups)
            for key, records in groups.items():
                _write_report_group(folder, stems[key], f"{label} — {key}", records, book, include_html)
        ready.rename(target)
    return target
