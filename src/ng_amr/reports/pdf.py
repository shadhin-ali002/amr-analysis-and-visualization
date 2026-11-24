"""Low-level PDF drawing helpers for mutation reports."""
from __future__ import annotations

from ng_amr.genomics.alignment import (
    A4,
    REPORTLAB_OK,
    HexColor,
    black,
    categorize_nt_mutations,
    inch,
    mm,
    white,
)

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


def _pdf_draw_kv_table(page, rows, col1_w=0.48):
    page.ensure(len(rows) * 5.5 + 3)
    row_h = 5.0 * mm
    x0 = _PDF_MARGIN
    w1 = _PDF_CONTENT_W * col1_w
    for i, (lab, val) in enumerate(rows):
        page.ensure(6)
        bg = _PDF_LIGHT if i % 2 == 0 else white
        page.c.setFillColor(bg)
        page.c.rect(x0, page.y - row_h, _PDF_CONTENT_W, row_h, fill=1, stroke=0)
        page.c.setStrokeColor(HexColor("#e0e0e0"))
        page.c.setLineWidth(0.3)
        page.c.rect(x0, page.y - row_h, _PDF_CONTENT_W, row_h, fill=0, stroke=1)
        page.c.setFillColor(black)
        page.c.setFont("Helvetica-Bold", 9)
        page.c.drawString(x0 + 1.5 * mm, page.y - row_h + 1.4 * mm, str(lab)[:50])
        page.c.setFont("Helvetica", 9)
        page.c.drawString(x0 + w1 + 1.5 * mm, page.y - row_h + 1.4 * mm, str(val)[:60])
        page.y -= row_h
    page.spacer(2.5)
