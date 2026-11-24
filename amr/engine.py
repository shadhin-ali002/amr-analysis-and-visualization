"""Engine routines."""
from __future__ import annotations

import matplotlib

from .inputs import (
    AD_CAT,
    AD_MCF,
    AD_MIC,
    ALIASES,
    ANTIBIOTIC_ORDER,
    ANTIBIOTIC_SHORT,
    Audit,
    CATEGORIES,
    Counter,
    DATE_FIELDS,
    EST_CAT,
    EST_MCF,
    MIC_FIELDS,
    Phenotypes,
    QC,
    ROUNDED,
    WGS_CAT,
    assembly_id,
    canonical,
    csv,
    csv_write,
    defaultdict,
    discover_assemblies,
    load_aliases,
    logger,
    math,
    measurement,
    openpyxl,
    re,
    read_master,
)
from .mutation_reports import (
    A4,
    CODING,
    Dict,
    HexColor,
    List,
    Optional,
    PairwiseAligner,
    Path,
    REPORTLAB_OK,
    Seq,
    SeqIO,
    Tuple,
    _CACHE,
    _PDFPage,
    _PDF_CONTENT_W,
    _PDF_CURATED,
    _PDF_DISRUPT,
    _PDF_HEADER,
    _PDF_LIGHT,
    _PDF_MARGIN,
    _PDF_MUTED,
    _PDF_NONCUR,
    _pdf_write_seq,
    _tempfile,
    analyze_nucleotide,
    analyze_protein,
    black,
    call_protein_mutations,
    categorize_nt_mutations,
    grouped_mutation_reports,
    json,
    mm,
    os,
    rl_canvas,
    safe_name,
    white,
)
from .visualization import (
    MATPLOTLIB_OK,
    _fig_short_label,
    _mutation_sort_key,
    create_master_analysis_figures,
    plt,
)

import atexit

import hashlib

import shutil

import subprocess

import sys

import statistics as _stats

import numpy as np

import pandas as pd

try:
    import parasail as _parasail
except ImportError:
    _parasail = None

ORIGINAL_AA_FOLDER = "Assembly/AA"

AA_COMPARISON_FOLDER_NAME = "AA_Comparison"

EXTRACTED_SEQUENCES_FOLDER_NAME = "Extracted_Sequences"

GENE_AA_PRODUCT_SYNONYMS: Dict[str, List[str]] = {
    "penA": [
        "penicillin-binding protein 2",
        "penicillin binding protein 2",
        "peptidoglycan d,d-transpeptidase pbp2",
        "pbp2",
        "pbp-2",
    ],
    "ponA1": [
        "penicillin-binding protein 1a",
        "penicillin binding protein 1a",
        "penicillin-binding protein 1",
        "pbp1a",
        "pbp-1a",
        "pona",
    ],
    "porB1b": [
        "trimeric porin porb.ib",
        "trimeric porin porb",
        "porb.ib",
        "porb1b",
        "porin porb",
        "outer membrane protein porb",
    ],
    "blaTEM": [
        "tem family class a beta-lactamase",
        "tem family",
        "class a beta-lactamase",
        "beta-lactamase tem",
        "blatem",
        "tem-1",
        "tem1",
    ],
    "parC": [
        "dna topoisomerase iv subunit a",
        "topoisomerase iv subunit a",
        "parc",
    ],
    "parE": [
        "dna topoisomerase iv subunit b",
        "topoisomerase iv subunit b",
        "pare",
    ],
    "gyrA": [
        "dna gyrase subunit a",
        "gyrase subunit a",
        "gyra",
    ],
    "mtrR": [
        "multidrug efflux system transcriptional repressor mtrr",
        "transcriptional repressor mtrr",
        "mtrr",
    ],
    "mtrC": [
        "multidrug efflux rnd transporter periplasmic adaptor subunit mtrc",
        "periplasmic adaptor subunit mtrc",
        "mtrc",
    ],
    "mtrD": [
        "multidrug efflux rnd transporter permease subunit mtrd",
        "permease subunit mtrd",
        "mtrd",
    ],
    "rpsJ": [
        "30s ribosomal protein s10",
        "ribosomal protein s10",
        "rpsj",
    ],
    "folP": [
        "dihydropteroate synthase",
        "dihydropteroate",
        "folp",
    ],
    "rpsE": [
        "30s ribosomal protein s5",
        "ribosomal protein s5",
        "rpse",
    ],
    "rpoB": [
        "dna-directed rna polymerase subunit beta",
        "rna polymerase subunit beta",
        "rpob",
    ],
    "rpoD": [
        "rna polymerase sigma factor rpod",
        "sigma factor rpod",
        "rpod",
    ],
    "rplD": [
        "50s ribosomal protein l4",
        "ribosomal protein l4",
        "rpld",
    ],
    "rplV": [
        "50s ribosomal protein l22",
        "ribosomal protein l22",
        "rplv",
    ],
    "macAB_promoter": [],
    "mtrR_promoter": [],
    "16S_rDNA": [],
    "23S_rDNA": [],
    "ermA": ["ermA", "rRNA methyltransferase ermA", "23s rrna methyltransferase erma"],
    "ermB": ["ermB", "rRNA methyltransferase ermB", "23s rrna methyltransferase ermb"],
    "ermC": ["ermC", "rRNA methyltransferase ermC", "23s rrna methyltransferase ermc"],
    "ermF": ["ermF", "rRNA methyltransferase ermF", "23s rrna methyltransferase ermf"],
    "ereA": ["ereA", "erythromycin esterase", "erythromycin esterase erea"],
    "ereB": ["ereB", "erythromycin esterase ereb"],
    "mefA": ["mefA", "macrolide efflux", "mef(a)", "macrolide-efflux"],
}

DEFAULT_GENES = [
    "penA", "ponA1", "porB1b", "blaTEM", "parC", "parE", "gyrA", "mtrR", "mtrR_promoter", "mtrC", "mtrD",
    "macAB_promoter", "rpsJ", "folP", "rpsE", "16S_rDNA", "rpoB", "rpoD", "rplD", "rplV", "23S_rDNA",
    "ermA", "ermB", "ermC", "ermF", "ereA", "ereB", "mefA",
    "norM_promoter", "tetM",
]

GENE_FOLDER_NAMES = {
    "penA":           "1. penA",
    "ponA1":          "2. ponA1",
    "porB1b":         "3. porB1b",
    "blaTEM":         "4. blaTEM",
    "parC":           "5. parC",
    "parE":           "6. parE",
    "gyrA":           "7. gyrA",
    "mtrR":           "8. mtrR",
    "mtrR_promoter":  "9. mtrR_promoter",
    "mtrC":           "10. mtrC",
    "mtrD":           "11. mtrD",
    "macAB_promoter": "12. macAB_promoter",
    "rpsJ":           "13. rpsJ",
    "folP":           "14. folP",
    "rpsE":           "15. rpsE",
    "16S_rDNA":       "16. 16S_rDNA",
    "rpoB":           "17. rpoB",
    "rpoD":           "18. rpoD",
    "rplD":           "19. rplD",
    "rplV":           "20. rplV",
    "23S_rDNA":       "21. 23S_rDNA",
    "ermA":           "22. ermA",
    "ermB":           "23. ermB",
    "ermC":           "24. ermC",
    "ermF":           "25. ermF",
    "ereA":           "26. ereA",
    "ereB":           "27. ereB",
    "mefA":           "28. mefA",
    "norM_promoter":  "29. norM_promoter",
    "tetM":           "30. tetM",
}

GENE_CONFIG = {
    "penA":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ponA1":           {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "porB1b":          {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "blaTEM":          {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "gyrA":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "parC":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "parE":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "mtrR":            {"gene_type": "CDS",       "multi_ref": False, "disruption": True},
    "mtrR_promoter":   {"gene_type": "noncoding", "multi_ref": True,  "disruption": True},
    "mtrC":            {"gene_type": "CDS",       "multi_ref": False, "disruption": True},
    "mtrD":            {"gene_type": "CDS",       "multi_ref": True,  "disruption": True},
    "macAB_promoter":  {"gene_type": "noncoding", "multi_ref": True,  "disruption": True},
    "norM_promoter":   {"gene_type": "noncoding", "multi_ref": False, "disruption": False},
    "23S_rDNA":        {"gene_type": "rRNA",      "multi_ref": False, "disruption": False},
    "16S_rDNA":        {"gene_type": "rRNA",      "multi_ref": False, "disruption": False},
    "tetM":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rpsJ":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rpoB":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "folP":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rpsE":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rpoD":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rplD":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rplV":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ermA":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ermB":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ermC":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ermF":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ereA":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ereB":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "mefA":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
}

DISRUPTION_MIN_DELETION_AA = 15

DISRUPTION_MIN_DELETION_EVENTS = 3

DISRUPTION_MAX_LENGTH_CHANGE = 0.20

def short_isolate_label(iso: str, max_len: int = 18) -> str:
    m = re.search(r"IsolateID[_-]?(.+?)(?:_GC[AF]\d|$)", iso or "", re.I)
    if m:
        return m.group(1).strip("_")[:max_len]
    m = re.match(r"^(\d+)", iso or "")
    return m.group(1) if m else (iso or "")[:max_len]

def get_gene_type(gene_name: str) -> str:
    return GENE_CONFIG.get(gene_name, {}).get("gene_type", "CDS")

def _gene_token_in_text(text: str, gene_name: str) -> bool:
    if not text or not gene_name:
        return False
    return re.search(
        r"(?<![A-Za-z0-9_])" + re.escape(gene_name) + r"(?![A-Za-z0-9_])",
        text,
        flags=re.I,
    ) is not None

def assess_disruption(protein_result: dict) -> bool:
    for m in protein_result["deletions"]:
        if m.get("length", 0) >= DISRUPTION_MIN_DELETION_AA:
            return True
    if len(protein_result["deletions"]) >= DISRUPTION_MIN_DELETION_EVENTS:
        return True
    if protein_result["stop_gains"] or protein_result["frameshifts"]:
        return True

    def _orf_len(prot: str) -> int:
        prot = prot or ""
        cut = prot.find("*")
        return len(prot) if cut < 0 else cut

    ref_len = _orf_len(protein_result["ref_protein"])
    qry_len = _orf_len(protein_result["query_protein"])
    if ref_len > 0 and abs(ref_len - qry_len) / ref_len > DISRUPTION_MAX_LENGTH_CHANGE:
        return True
    return False

def normalize_sample_id(value):
    return assembly_id(value).lower()

def load_curated(path: str) -> dict:
    with open(path) as f:
        return json.load(f)

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

def create_overall_mutation_charts(mutation_counts, total_isolates, curated_set,
                                   output_dir, prefix, gene_name="penA"):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.info("  Skipping overall mutation charts: matplotlib not installed.")
        return None

    if total_isolates <= 0 or not mutation_counts:
        return None

    items_sorted = sorted(mutation_counts.items(), key=lambda x: -x[1])
    all_labels = [mut for mut, cnt in items_sorted]
    all_counts = [cnt for mut, cnt in items_sorted]
    all_pcts = [100.0 * cnt / total_isolates for cnt in all_counts]
    curated_flags = [mut in curated_set for mut in all_labels]
    curated_color = "#2e7d32"
    noncurated_color = "#e65100"

    os.makedirs(output_dir, exist_ok=True)

    def save_fig(fig, filename):
        fig.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, max(6, len(all_labels)*0.3)))
    y_pos = range(len(all_labels))
    width = 0.35
    curated_vals = [pct if flag else 0 for pct, flag in zip(all_pcts, curated_flags)]
    noncurated_vals = [pct if not flag else 0 for pct, flag in zip(all_pcts, curated_flags)]
    ax.barh([i - width/2 for i in y_pos], curated_vals, height=width, color=curated_color, label='Curated')
    ax.barh([i + width/2 for i in y_pos], noncurated_vals, height=width, color=noncurated_color, label='Non-curated')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(all_labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("% of isolates")
    ax.set_title(f"{gene_name}: ALL Mutations (n={total_isolates} isolates)", fontweight='bold')
    if ax.get_legend_handles_labels()[0]:
        ax.legend()
    ax.grid(axis='x', linestyle='--', alpha=0.3)
    save_fig(fig, f"{prefix}_overall_1_grouped_bar.png")

    fig, ax = plt.subplots(figsize=(12, max(6, len(all_labels)*0.3)))
    for _i, _pct in enumerate(all_pcts):
        ax.plot([0, _pct], [_i, _i], color="gray", linewidth=1, alpha=0.6)
    for i, (pct, flag) in enumerate(zip(all_pcts, curated_flags)):
        color = curated_color if flag else noncurated_color
        ax.plot(pct, i, 'o', markersize=10, color=color)
        ax.text(pct + 1, i, f"{pct:.0f}% ({all_counts[i]}/{total_isolates})", va='center', fontsize=8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(all_labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("% of isolates")
    ax.set_title(f"{gene_name}: ALL Mutations Lollipop Plot (n={total_isolates} isolates)", fontweight='bold')
    if ax.get_legend_handles_labels()[0]:
        ax.legend()
    save_fig(fig, f"{prefix}_overall_2_lollipop.png")

    fig, ax = plt.subplots(figsize=(12, max(6, len(all_labels)*0.3)))
    for i, (pct, flag) in enumerate(zip(all_pcts, curated_flags)):
        color = curated_color if flag else noncurated_color
        ax.barh(i, pct, color=color, alpha=0.7, height=0.6)
        ax.text(0.5, i, f"{all_labels[i]}  ({pct:.0f}%)", va='center', ha='left', fontsize=9,
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))
    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.set_xlabel("% of isolates")
    ax.set_title(f"{gene_name}: ALL Mutations Strip (color intensity = prevalence, n={total_isolates} isolates)", fontweight='bold')
    save_fig(fig, f"{prefix}_overall_3_heatmap_strip.png")

    fig, ax = plt.subplots(figsize=(12, max(6, len(all_labels)*0.3)))
    for i, (pct, cnt, flag) in enumerate(zip(all_pcts, all_counts, curated_flags)):
        color = curated_color if flag else noncurated_color
        ax.plot([0, pct], [i, i], color='gray', linewidth=1, linestyle='-')
        ax.plot(pct, i, 'o', markersize=10, color=color)
        ax.text(pct + 1, i, f"{pct:.0f}% ({cnt}/{total_isolates})", va='center', fontsize=8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(all_labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("% of isolates")
    ax.set_title(f"{gene_name}: ALL Mutations Dumbbell Plot (n={total_isolates} isolates)", fontweight='bold')
    ax.grid(axis='x', linestyle='--', alpha=0.3)
    save_fig(fig, f"{prefix}_overall_4_dumbbell.png")

    curated_items = [(m, c, p) for m, c, p, f in zip(all_labels, all_counts, all_pcts, curated_flags) if f]
    noncurated_items = [(m, c, p) for m, c, p, f in zip(all_labels, all_counts, all_pcts, curated_flags) if not f]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(6, len(all_labels)*0.3)))
    if curated_items:
        curated_labels = [m for m, c, p in curated_items]
        curated_pcts = [p for m, c, p in curated_items]
        y1 = range(len(curated_items))
        ax1.barh(y1, curated_pcts, color=curated_color, height=0.6)
        ax1.set_yticks(y1)
        ax1.set_yticklabels(curated_labels, fontsize=9)
        ax1.invert_yaxis()
        ax1.set_title(f"Curated ({len(curated_items)})", fontweight='bold')
        ax1.set_xlabel("%")
        for i, pct in enumerate(curated_pcts):
            ax1.text(pct+1, i, f"{pct:.0f}%", va='center', fontsize=8)
    else:
        ax1.axis('off')
        ax1.text(0.5, 0.5, "No curated mutations", ha='center', va='center', fontsize=11, color='gray')

    if noncurated_items:
        non_labels = [m for m, c, p in noncurated_items]
        non_pcts = [p for m, c, p in noncurated_items]
        y2 = range(len(noncurated_items))
        ax2.barh(y2, non_pcts, color=noncurated_color, height=0.6)
        ax2.set_yticks(y2)
        ax2.set_yticklabels(non_labels, fontsize=9)
        ax2.invert_yaxis()
        ax2.set_title(f"Non-curated ({len(noncurated_items)})", fontweight='bold')
        ax2.set_xlabel("%")
        for i, pct in enumerate(non_pcts):
            ax2.text(pct+1, i, f"{pct:.0f}%", va='center', fontsize=8)
    else:
        ax2.axis('off')
        ax2.text(0.5, 0.5, "No non-curated mutations", ha='center', va='center', fontsize=11, color='gray')

    fig.suptitle(f"{gene_name}: ALL Mutations Card (n={total_isolates} isolates)", fontweight='bold')
    save_fig(fig, f"{prefix}_overall_5_mutation_card.png")

    return output_dir

def create_gene_mutation_cooccurrence_figure(gene_results: List[dict], output_root: str) -> Optional[str]:
    if not MATPLOTLIB_OK:
        logger.info("  matplotlib missing – skip gene mutation co-occurrence figure")
        return None
    if not gene_results:
        return None

    import numpy as np
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch

    isolate_set = set()
    for gr in gene_results:
        for d in gr.get("all_data") or []:
            isolate_set.add(d["isolate"])
        for d in gr.get("all_mutation_data") or []:
            isolate_set.add(d["isolate"])

    def _iso_order(name: str):
        m = re.match(r"^(\d+)", (name or "").strip())
        idx = int(m.group(1)) if m else 10**9
        return (idx, name)

    isolates = sorted(isolate_set, key=_iso_order)
    genes = [gr["gene"] for gr in gene_results]
    if not isolates or not genes:
        return None

    n_iso = len(isolates)
    n_genes = len(genes)
    iso_idx = {iso: i for i, iso in enumerate(isolates)}

    matrix = np.zeros((n_iso, n_genes), dtype=int)

    for gi, gr in enumerate(gene_results):
        data_by_iso = {d["isolate"]: d for d in (gr.get("all_data") or [])}
        mutdata_by_iso = {
            d["isolate"]: set(d.get("mutations") or [])
            for d in (gr.get("all_mutation_data") or [])
        }
        for iso in isolates:
            ad = data_by_iso.get(iso)
            if ad is None or not bool(ad.get("gene_found", ad.get("status") == "OK")):
                matrix[iso_idx[iso], gi] = 2
            elif mutdata_by_iso.get(iso):
                matrix[iso_idx[iso], gi] = 1

    mutation_gene_burden = (matrix == 1).sum(axis=1)
    missing_burden = (matrix == 2).sum(axis=1)
    iso_order = np.lexsort((np.arange(n_iso), -missing_burden, -mutation_gene_burden))
    matrix = matrix[iso_order]
    ordered_isolates = [isolates[i] for i in iso_order]

    def _short(iso: str) -> str:
        m = re.search(r"IsolateID[_-]?(.+?)(?:_GC[AF]\d|$)", iso, re.I)
        if m:
            return m.group(1).strip("_")[:18]
        m = re.match(r"^(\d+)", iso)
        return m.group(1) if m else iso[:18]

    short_isolates = [_short(i) for i in ordered_isolates]
    gene_mutation_frequency = (matrix == 1).sum(axis=0)
    isolate_mutated_gene_count = (matrix == 1).sum(axis=1)

    mutation_binary = (matrix == 1).astype(int)
    pairwise_cooccurrence = mutation_binary.T @ mutation_binary
    np.fill_diagonal(pairwise_cooccurrence, gene_mutation_frequency)

    out_dir = os.path.join(output_root, "master_analysis_figures")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(
        out_dir, "14_gene_mutation_cooccurrence_per_isolate.png"
    )

    fig_h = max(9.5, 5.5 + n_iso * 0.30 + n_genes * 0.34)
    fig_w = max(11.5, n_genes * 0.78 + 4.0)

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(
        3, 2,
        width_ratios=[n_genes, 2.2],
        height_ratios=[
            1.7,
            max(4.5, n_iso * 0.30),
            max(3.6, n_genes * 0.52),
        ],
        wspace=0.08,
        hspace=0.28,
    )

    ax_top = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[1, 0])
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_heat)
    ax_pair = fig.add_subplot(gs[2, 0])

    ax_top.bar(
        range(n_genes),
        gene_mutation_frequency,
        width=0.8,
        color="#1565c0",
    )
    ax_top.set_xlim(-0.5, n_genes - 0.5)
    ax_top.set_ylabel("No. isolates", fontsize=9)
    ax_top.set_xticks(range(n_genes))
    ax_top.set_xticklabels(genes, rotation=45, ha="left", fontsize=8)
    ax_top.set_title(
        "Gene Mutation Co-occurrence per Isolate",
        fontweight="bold",
        pad=12,
    )
    ax_top.grid(axis="y", linestyle="--", alpha=0.25)
    for x, value in enumerate(gene_mutation_frequency):
        if value > 0:
            ax_top.text(
                x,
                value + max(0.15, n_iso * 0.015),
                str(int(value)),
                ha="center",
                va="bottom",
                fontsize=7,
            )

    cmap = ListedColormap(["#ffffff", "#1565c0", "#bdbdbd"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    ax_heat.imshow(
        matrix,
        aspect="auto",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
    )
    ax_heat.set_yticks(range(n_iso))
    ax_heat.set_yticklabels(short_isolates, fontsize=7)
    ax_heat.set_xticks(range(n_genes))
    ax_heat.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax_heat.set_ylabel("Isolates", fontsize=10)
    ax_heat.set_xlabel("Genes", fontsize=10)
    ax_heat.set_title(
        "Mutation status: isolate × gene",
        fontsize=10,
        pad=8,
    )

    ax_right.barh(
        range(n_iso),
        isolate_mutated_gene_count,
        height=0.8,
        color="#1565c0",
    )
    ax_right.set_xlabel("Mutated genes", fontsize=9)
    ax_right.grid(axis="x", linestyle="--", alpha=0.25)
    ax_right.tick_params(axis="y", left=False, labelleft=False)
    ax_right.invert_yaxis()
    max_count = max(
        1,
        int(isolate_mutated_gene_count.max()) if n_iso else 1,
    )
    ax_right.set_xlim(0, max_count + max(1, max_count * 0.15))
    for y, value in enumerate(isolate_mutated_gene_count):
        ax_right.text(
            value + 0.08,
            y,
            str(int(value)),
            va="center",
            fontsize=7,
        )

    pair_im = ax_pair.imshow(
        pairwise_cooccurrence,
        aspect="auto",
        interpolation="nearest",
        cmap="Blues",
        vmin=0,
        vmax=max(1, int(pairwise_cooccurrence.max())),
    )
    ax_pair.set_xticks(range(n_genes))
    ax_pair.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax_pair.set_yticks(range(n_genes))
    ax_pair.set_yticklabels(genes, fontsize=8)
    ax_pair.set_xlabel("Gene", fontsize=10)
    ax_pair.set_ylabel("Gene", fontsize=10)
    ax_pair.set_title(
        "Pairwise mutation co-occurrence: number of isolates carrying both genes",
        fontsize=10,
        fontweight="bold",
        pad=10,
    )

    if n_genes <= 20:
        vmax_pair = max(1, int(pairwise_cooccurrence.max()))
        for i in range(n_genes):
            for j in range(n_genes):
                value = int(pairwise_cooccurrence[i, j])
                text_color = "white" if value > vmax_pair * 0.55 else "black"
                ax_pair.text(
                    j,
                    i,
                    str(value),
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=text_color,
                )

    cbar = fig.colorbar(
        pair_im,
        ax=ax_pair,
        fraction=0.046,
        pad=0.04,
    )
    cbar.set_label("No. isolates with mutations in both genes", fontsize=8)

    legend_handles = [
        Patch(
            facecolor="#1565c0",
            edgecolor="none",
            label="≥1 mutation",
        ),
        Patch(
            facecolor="#ffffff",
            edgecolor="#777777",
            label="No mutation detected",
        ),
        Patch(
            facecolor="#bdbdbd",
            edgecolor="none",
            label="Gene not found / not called",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, 0.005),
        frameon=False,
        fontsize=8,
    )

    fig.savefig(
        out_path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)

    logger.info(f"  Gene mutation co-occurrence figure → {out_path}")
    return out_path

def _parse_lab_id(value):
    return canonical(value)

def _parse_accession(full_name: str) -> str:
    m = re.search(r"(GC[AF]\d[\d.]*)", full_name, re.IGNORECASE)
    return m.group(1) if m else ""

def _parse_assembly_name(full_name: str) -> str:
    m = re.search(r"(ASM\S+)", full_name, re.IGNORECASE)
    return m.group(1) if m else ""

def write_master_matrix(gene_results: List[dict], output_root: str, thresholds: List[int], book):
    if pd is None:
        logger.info("  pandas missing – skip master matrix")
        return None

    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
        from openpyxl.utils import get_column_letter
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

    thin = Border(
        left=Side(style="thin", color="CCCCCC"),
        right=Side(style="thin", color="CCCCCC"),
        top=Side(style="thin", color="CCCCCC"),
        bottom=Side(style="thin", color="CCCCCC"),
    )
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
        n_called = sum(1 for iso in isolates if gene_iso[g][iso]["found"])
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
        n_called = sum(1 for iso in isolates if gene_iso[g][iso]["found"])
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

def create_multigene_mutation_cooccurrence_count_heatmap(
    gene_results: List[dict], output_root: str
) -> Optional[str]:
    if not MATPLOTLIB_OK:
        logger.info("  matplotlib missing – skip multi-gene co-occurrence count heatmap")
        return None
    if not gene_results:
        return None

    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap
    from collections import Counter

    isolate_set = set()
    for gr in gene_results:
        for d in gr.get("all_data") or []:
            isolate_set.add(d["isolate"])
        for d in gr.get("all_mutation_data") or []:
            isolate_set.add(d["isolate"])

    def _iso_order(name: str):
        m = re.match(r"^(\d+)", (name or "").strip())
        idx = int(m.group(1)) if m else 10**9
        return (idx, name)

    isolates = sorted(isolate_set, key=_iso_order)
    genes = [gr["gene"] for gr in gene_results]
    if not isolates or not genes:
        return None

    n_iso = len(isolates)
    n_genes = len(genes)
    iso_idx = {iso: i for i, iso in enumerate(isolates)}

    matrix = np.zeros((n_iso, n_genes), dtype=int)
    for gi, gr in enumerate(gene_results):
        data_by_iso = {d["isolate"]: d for d in (gr.get("all_data") or [])}
        mutdata_by_iso = {
            d["isolate"]: set(d.get("mutations") or [])
            for d in (gr.get("all_mutation_data") or [])
        }
        for iso in isolates:
            ad = data_by_iso.get(iso)
            found = bool(ad and ad.get("gene_found", ad.get("status") == "OK"))
            if found and mutdata_by_iso.get(iso):
                matrix[iso_idx[iso], gi] = 1

    cooc = matrix.T @ matrix

    pattern_counter: Counter = Counter()
    pattern_isolates: Dict[tuple, List[str]] = {}
    for i, iso in enumerate(isolates):
        key = tuple(int(x) for x in matrix[i])
        if sum(key) == 0:
            continue
        pattern_counter[key] += 1
        pattern_isolates.setdefault(key, []).append(iso)

    ranked_patterns = sorted(
        pattern_counter.items(), key=lambda x: (-x[1], -sum(x[0]), x[0])
    )

    out_dir = os.path.join(output_root, "master_analysis_figures")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(
        out_dir, "19_multigene_mutation_cooccurrence_count_heatmap.png"
    )

    n_top = min(12, len(ranked_patterns))
    fig_h = max(7.5, n_genes * 0.55 + 2.5)
    fig_w = max(11.0, n_genes * 0.75 + 5.5)

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(
        1, 2,
        width_ratios=[max(5.5, n_genes * 0.72), 4.2],
        wspace=0.28,
    )
    ax_hm = fig.add_subplot(gs[0, 0])
    ax_tbl = fig.add_subplot(gs[0, 1])
    ax_tbl.axis("off")

    vmax = max(1, int(cooc.max()))
    cmap = LinearSegmentedColormap.from_list(
        "cooc_blues",
        ["#ffffff", "#bbdefb", "#64b5f6", "#1976d2", "#0d47a1"],
        N=256,
    )
    im = ax_hm.imshow(
        cooc,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=vmax,
    )
    ax_hm.set_xticks(range(n_genes))
    ax_hm.set_xticklabels(genes, rotation=45, ha="right", fontsize=9)
    ax_hm.set_yticks(range(n_genes))
    ax_hm.set_yticklabels(genes, fontsize=9)
    ax_hm.set_xlabel("Gene", fontsize=10)
    ax_hm.set_ylabel("Gene", fontsize=10)
    ax_hm.set_title(
        "Multi-gene mutation co-occurrence\n"
        "(cell = number of isolates with ≥1 mutation in BOTH genes)",
        fontsize=11,
        fontweight="bold",
        pad=10,
    )

    for i in range(n_genes):
        for j in range(n_genes):
            val = int(cooc[i, j])
            text_color = "white" if val > vmax * 0.55 else "black"
            weight = "bold" if i == j else "normal"
            ax_hm.text(
                j, i, str(val),
                ha="center", va="center",
                fontsize=8, color=text_color, fontweight=weight,
            )

    cbar = fig.colorbar(im, ax=ax_hm, fraction=0.046, pad=0.04)
    cbar.set_label("No. isolates with co-occurring mutations", fontsize=8)

    lines = [
        "Most frequent multi-gene",
        "mutation combinations",
        "(exact gene profile → count)",
        "",
    ]
    for rank, (pattern, count) in enumerate(ranked_patterns[:n_top], 1):
        gene_names = [genes[k] for k, bit in enumerate(pattern) if bit]
        if not gene_names:
            continue
        if len(gene_names) == 1:
            label = gene_names[0] + " only"
        else:
            label = " + ".join(gene_names)
        if len(label) > 42:
            label = label[:39] + "..."
        lines.append(f"{rank:2d}. {label}")
        lines.append(f"     → {count} isolate{'s' if count != 1 else ''}")

    if n_top == 0:
        lines.append("(no mutations detected)")

    y0 = 0.98
    for i, line in enumerate(lines):
        is_header = i < 3
        ax_tbl.text(
            0.02, y0 - i * 0.055,
            line,
            transform=ax_tbl.transAxes,
            fontsize=8 if not is_header else 9,
            fontweight="bold" if is_header or line.startswith("     →") else "normal",
            family="monospace" if not is_header else "sans-serif",
            va="top",
            color="#1565c0" if line.startswith("     →") else "#212121",
        )

    fig.suptitle(
        f"Gene mutation co-occurrence counts across {n_iso} isolates",
        fontsize=12,
        fontweight="bold",
        y=1.02,
    )

    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    logger.info(f"  Multi-gene mutation co-occurrence count heatmap → {out_path}")
    return out_path

_DOMAIN_CONFIG = None

_DEFAULT_NG_AMR_LITERATURE_HOTSPOTS: Dict[str, set] = {
    "penA":        {"A501V", "A501T", "A502T", "A502V", "G542S", "G545S",
                     "P551S", "P551L", "A311V"},
    "ponA1":       {"L421P"},
    "porB1b":      {"G120K", "G120D", "A121D", "A121N", "A121G",
                     "G101K", "G101T", "A102D", "A102N"},
    "gyrA":        {"S91F", "D95A", "D95G", "D95N", "D95Y"},
    "parC":        {"S87R", "S87W", "S87N", "S87I", "D86N", "E91K", "S88P"},
    "mtrR":        {"A39T", "G45D", "H105Y", "A40D"},
    "mtrR_promoter": set(),
    "rpsJ":        {"V57M"},
    "23S_rDNA":    {"A2059G", "C2611T", "C2599T"},
    "rpsE":        {"T24P"},
    "16S_rDNA":    {"C1192T"},
    "blaTEM":      {"M182T"},
    "mtrC":        set(),
    "mtrD":        set(),
    "macAB_promoter": set(),
    "rplD":        set(),
    "rplV":        set(),
    "rpoB":        set(),
    "rpoD":        set(),
    "folP":        set(),
    "ermA":        set(),
    "ermB":        set(),
    "ermC":        set(),
    "ermF":        set(),
    "ereA":        set(),
    "ereB":        set(),
    "mefA":        set(),
    "tetM":        set(),
}

_DEFAULT_AA_PROPERTY_GROUPS: Dict[str, str] = {
    "G": "Nonpolar", "A": "Nonpolar", "V": "Nonpolar", "L": "Nonpolar",
    "I": "Nonpolar", "M": "Nonpolar", "P": "Nonpolar",
    "F": "Aromatic", "W": "Aromatic", "Y": "Aromatic",
    "S": "Polar", "T": "Polar", "C": "Polar", "N": "Polar", "Q": "Polar",
    "D": "Acidic", "E": "Acidic",
    "K": "Basic", "R": "Basic", "H": "Basic",
}

_DEFAULT_PUB_PALETTE = {
    "primary": "#1565c0", "curated": "#2e7d32", "noncurated": "#e65100",
    "accent": "#6a1b9a", "warn": "#c62828", "neutral": "#607d8b",
    "grid": "#e0e0e0",
}

NG_AMR_LITERATURE_HOTSPOTS: Dict[str, set] = (
    _DOMAIN_CONFIG["literature_hotspots"] if _DOMAIN_CONFIG else _DEFAULT_NG_AMR_LITERATURE_HOTSPOTS
)

AA_PROPERTY_GROUPS: Dict[str, str] = (
    _DOMAIN_CONFIG["aa_property_groups"] if _DOMAIN_CONFIG else _DEFAULT_AA_PROPERTY_GROUPS
)

PUB_PALETTE = (
    _DOMAIN_CONFIG["pub_palette"] if _DOMAIN_CONFIG else _DEFAULT_PUB_PALETTE
)

def _pub_rc():
    if not MATPLOTLIB_OK:
        return
    matplotlib.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "font.family": "DejaVu Sans",
        "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold",
        "axes.labelsize": 10, "axes.edgecolor": "#424242", "axes.linewidth": 0.8,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": "#e0e0e0", "grid.linewidth": 0.5,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })

def _safe_savefig(fig, path: str):
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    finally:
        plt.close(fig)

def _no_data_figure(path: str, title: str, note: str = "No qualifying mutations for this view"):
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axis("off")
    ax.set_title(title, pad=14)
    ax.text(0.5, 0.5, note, ha="center", va="center", fontsize=11,
            color="#757575", transform=ax.transAxes)
    _safe_savefig(fig, path)

def _wilson_ci(count: int, n: int, z: float = 1.96):
    if n <= 0:
        return 0.0, 0.0, 0.0
    p = count / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z * math.sqrt((p * (1 - p) / n) + (z * z / (4 * n * n)))) / denom
    return p, max(0.0, centre - half), min(1.0, centre + half)

def _get_literature_hotspots(gene_name: str) -> set:
    return set(NG_AMR_LITERATURE_HOTSPOTS.get(gene_name, set()))

def _split_relevant(gene_name: str, mutation_counts: Dict[str, int], curated_set) -> Tuple[Dict[str, int], Dict[str, int]]:
    relevant_names = set(curated_set or set()) | _get_literature_hotspots(gene_name)
    all_counts = dict(mutation_counts or {})
    relevant_counts = {k: v for k, v in all_counts.items() if k in relevant_names}
    return all_counts, relevant_counts

def _notation_metadata(all_data: List[dict]):
    position_map, type_map = {}, {}
    for entry in all_data or []:
        for m in entry.get("mutations", []) or []:
            notation = m.get("mutation")
            if not notation:
                continue
            if notation not in position_map:
                position_map[notation] = m.get("position", 0) or 0
            if notation not in type_map:
                type_map[notation] = m.get("type", "substitution")
    return position_map, type_map

def _classify_substitution(notation: str):
    m = re.match(r"^([A-Z])(\d+)([A-Z])$", notation or "")
    if not m:
        return None
    ref_aa, pos, alt_aa = m.group(1), int(m.group(2)), m.group(3)
    g1, g2 = AA_PROPERTY_GROUPS.get(ref_aa), AA_PROPERTY_GROUPS.get(alt_aa)
    if g1 is None or g2 is None:
        change = "Unclassified"
    else:
        change = "Conservative (same class)" if g1 == g2 else "Radical (class change)"
    return ref_aa, pos, alt_aa, change

def _chart_bar(counts: Dict[str, int], total: int, path: str, title: str,
                mode: str = "count", orientation: str = "h", top_n: Optional[int] = None,
                color=None, xlabel=None):
    color = color or PUB_PALETTE["primary"]
    if not counts:
        return _no_data_figure(path, title)
    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    if top_n:
        items = items[:top_n]
    labels = [k for k, _ in items]
    if mode == "percent":
        values = [100.0 * v / total if total else 0 for _, v in items]
        xlabel = xlabel or "Isolates carrying mutation (%)"
    else:
        values = [v for _, v in items]
        xlabel = xlabel or "Number of isolates"
    n = len(labels)
    fig, ax = plt.subplots(figsize=(8, max(3, 0.32 * n + 1)) if orientation == "h" else (max(6, 0.4 * n), 5))
    if orientation == "h":
        y = range(n)
        ax.barh(y, values, color=color, edgecolor="white", height=0.7)
        ax.set_yticks(list(y)); ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel(xlabel)
    else:
        x = range(n)
        ax.bar(x, values, color=color, edgecolor="white", width=0.7)
        ax.set_xticks(list(x)); ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
        ax.set_ylabel(xlabel)
    ax.set_title(title)
    fig.text(0.99, 0.01, f"n = {total} isolates", ha="right", fontsize=7, color="#9e9e9e")
    _safe_savefig(fig, path)

def _chart_donut_categories(counts: Dict[str, int], total: int, path: str, title: str):
    if not counts or not total:
        return _no_data_figure(path, title)
    bins = {"Rare (<10%)": 0, "Low (10–25%)": 0, "Moderate (25–50%)": 0, "High (50–75%)": 0, "Majority (>75%)": 0}
    for _, c in counts.items():
        pct = 100.0 * c / total
        if pct < 10: bins["Rare (<10%)"] += 1
        elif pct < 25: bins["Low (10–25%)"] += 1
        elif pct < 50: bins["Moderate (25–50%)"] += 1
        elif pct < 75: bins["High (50–75%)"] += 1
        else: bins["Majority (>75%)"] += 1
    bins = {k: v for k, v in bins.items() if v > 0}
    if not bins:
        return _no_data_figure(path, title)
    colors = ["#90a4ae", "#4fc3f7", "#ffb74d", "#ef5350", "#8e24aa"][:len(bins)]
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.pie(bins.values(), labels=bins.keys(), autopct="%1.0f%%", colors=colors,
           startangle=90, pctdistance=0.78, wedgeprops=dict(width=0.42, edgecolor="white"))
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_isolate_mutation_heatmap(all_mutation_data: List[dict], counts: Dict[str, int],
                                     path: str, title: str, top_n: int = 25):
    if not counts or not all_mutation_data:
        return _no_data_figure(path, title)
    mut_order = [m for m, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]]
    isolates = [d["isolate"] for d in all_mutation_data]
    if not mut_order or not isolates:
        return _no_data_figure(path, title)
    import numpy as np
    mat = np.zeros((len(isolates), len(mut_order)))
    for i, d in enumerate(all_mutation_data):
        present = set(d.get("mutations", []))
        for j, m in enumerate(mut_order):
            mat[i, j] = 1.0 if m in present else 0.0
    fig_h = max(4, 0.18 * len(isolates) + 1.5)
    fig, ax = plt.subplots(figsize=(max(8, 0.4 * len(mut_order)), min(fig_h, 22)))
    ax.imshow(mat, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(mut_order))); ax.set_xticklabels(mut_order, rotation=90, fontsize=7)
    if len(isolates) <= 60:
        ax.set_yticks(range(len(isolates)))
        ax.set_yticklabels([short_isolate_label(x, 22) for x in isolates], fontsize=6)
    else:
        ax.set_yticks([])
        ax.set_ylabel(f"{len(isolates)} isolates")
    ax.set_title(title)
    ax.grid(False)
    _safe_savefig(fig, path)

def _chart_type_distribution(counts: Dict[str, int], type_map: Dict[str, str], path: str,
                              title: str, kind: str = "pie"):
    if not counts:
        return _no_data_figure(path, title)
    tally: Dict[str, int] = {}
    for m, c in counts.items():
        t = type_map.get(m, "substitution").replace("_", " ").title()
        tally[t] = tally.get(t, 0) + c
    colors = {"Substitution": "#1565c0", "Stop Gain": "#c62828", "Insertion": "#2e7d32",
              "Deletion": "#ef6c00", "Frameshift Insertion": "#6a1b9a", "Frameshift Deletion": "#4527a0"}
    cs = [colors.get(k, "#607d8b") for k in tally]
    fig, ax = plt.subplots(figsize=(6.5, 6.5) if kind == "pie" else (6.5, 4.5))
    if kind == "pie":
        ax.pie(tally.values(), labels=tally.keys(), autopct="%1.0f%%", colors=cs,
               startangle=90, wedgeprops=dict(edgecolor="white"))
    else:
        ax.bar(list(tally.keys()), list(tally.values()), color=cs, edgecolor="white")
        ax.set_ylabel("Isolate-mutation events"); plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_position_lollipop(counts: Dict[str, int], position_map: Dict[str, int], path: str, title: str):
    pts = [(position_map.get(m, 0), c, m) for m, c in counts.items() if position_map.get(m, 0)]
    if not pts:
        return _no_data_figure(path, title, "No positional (residue-level) mutations to plot")
    pts.sort()
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]; labels = [p[2] for p in pts]
    fig, ax = plt.subplots(figsize=(max(8, 0.15 * len(pts) + 4), 5))
    ax.vlines(xs, 0, ys, color=PUB_PALETTE["neutral"], linewidth=1.2, alpha=0.8)
    ax.plot(xs, ys, "o", color=PUB_PALETTE["accent"], markersize=6, zorder=3)
    for x, y, lab in zip(xs, ys, labels):
        if y == max(ys) or len(pts) <= 20:
            ax.annotate(lab, (x, y), textcoords="offset points", xytext=(0, 5), fontsize=6, ha="center")
    ax.set_xlabel("Residue / nucleotide position (gene-relative)")
    ax.set_ylabel("Number of isolates")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_position_histogram(counts: Dict[str, int], position_map: Dict[str, int], path: str, title: str):
    positions = []
    for m, c in counts.items():
        p = position_map.get(m, 0)
        if p:
            positions.extend([p] * c)
    if not positions:
        return _no_data_figure(path, title, "No positional data available")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bins = min(30, max(5, len(set(positions))))
    ax.hist(positions, bins=bins, color=PUB_PALETTE["primary"], edgecolor="white")
    ax.set_xlabel("Residue / nucleotide position (gene-relative)")
    ax.set_ylabel("Isolate-mutation events")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_cumulative_curve(counts: Dict[str, int], total: int, path: str, title: str):
    if not counts or not total:
        return _no_data_figure(path, title)
    items = sorted(counts.items(), key=lambda kv: -kv[1])
    cum = []
    running = 0
    for _, c in items:
        running = min(total, running + c)
        cum.append(100.0 * running / total)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(range(1, len(cum) + 1), cum, marker="o", markersize=3, color=PUB_PALETTE["primary"])
    ax.fill_between(range(1, len(cum) + 1), cum, alpha=0.15, color=PUB_PALETTE["primary"])
    ax.set_xlabel("Number of distinct mutations considered (rank order, most prevalent first)")
    ax.set_ylabel("Cumulative isolate coverage (%)")
    ax.set_ylim(0, 105)
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_rank_frequency(counts: Dict[str, int], path: str, title: str):
    if not counts:
        return _no_data_figure(path, title)
    values = sorted(counts.values(), reverse=True)
    ranks = range(1, len(values) + 1)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.plot(ranks, values, "o-", color=PUB_PALETTE["accent"], markersize=4, linewidth=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Mutation rank (log scale)")
    ax.set_ylabel("Isolate count (log scale)")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_radar_top(counts: Dict[str, int], total: int, path: str, title: str, top_n: int = 8):
    items = sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]
    if len(items) < 3:
        return _no_data_figure(path, title, "Fewer than 3 mutations — radar chart not informative")
    labels = [k for k, _ in items]
    values = [100.0 * v / total if total else 0 for _, v in items]
    values += values[:1]
    angles = [n / float(len(labels)) * 2 * math.pi for n in range(len(labels))]
    angles += angles[:1]
    fig, ax = plt.subplots(figsize=(6.5, 6.5), subplot_kw=dict(polar=True))
    ax.plot(angles, values, color=PUB_PALETTE["accent"], linewidth=1.5)
    ax.fill(angles, values, color=PUB_PALETTE["accent"], alpha=0.25)
    ax.set_xticks(angles[:-1]); ax.set_xticklabels(labels, fontsize=8)
    ax.set_title(title, pad=20)
    _safe_savefig(fig, path)

def _chart_isolate_metric_bar(all_data: List[dict], metric: str, path: str, title: str):
    rows = [(d["isolate"], d.get(metric)) for d in all_data
            if d.get("gene_found") and isinstance(d.get(metric), (int, float))]
    if not rows:
        return _no_data_figure(path, title)
    rows.sort(key=lambda r: -r[1])
    labels = [short_isolate_label(r[0], 18) for r in rows]
    values = [r[1] for r in rows]
    fig, ax = plt.subplots(figsize=(max(8, 0.18 * len(rows)), 4.5))
    ax.bar(range(len(rows)), values, color=PUB_PALETTE["primary"], width=0.8)
    ax.set_xticks([]) if len(rows) > 40 else ax.set_xticks(range(len(rows)))
    if len(rows) <= 40:
        ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel("Percent identity (%)" if metric == "pident" else "Query coverage (%)")
    ax.set_xlabel(f"{len(rows)} isolates, sorted")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_aa_substitution_heatmap(counts: Dict[str, int], path: str, title: str):
    import numpy as np
    aas = list("ACDEFGHIKLMNPQRSTVWY")
    idx = {a: i for i, a in enumerate(aas)}
    mat = None
    found_any = False
    mat = [[0] * len(aas) for _ in aas]
    for notation, c in counts.items():
        parsed = _classify_substitution(notation)
        if not parsed:
            continue
        ref_aa, _, alt_aa, _ = parsed
        if ref_aa in idx and alt_aa in idx:
            mat[idx[ref_aa]][idx[alt_aa]] += c
            found_any = True
    if not found_any:
        return _no_data_figure(path, title, "No single-residue missense substitutions to map")
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(mat, cmap="OrRd")
    ax.set_xticks(range(len(aas))); ax.set_xticklabels(aas, fontsize=8)
    ax.set_yticks(range(len(aas))); ax.set_yticklabels(aas, fontsize=8)
    ax.set_xlabel("Substituted (variant) residue"); ax.set_ylabel("Reference residue")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Isolate-mutation events")
    ax.grid(False)
    _safe_savefig(fig, path)

def _chart_physicochemical_bar(counts: Dict[str, int], path: str, title: str):
    tally: Dict[str, int] = {}
    any_found = False
    for notation, c in counts.items():
        parsed = _classify_substitution(notation)
        if not parsed:
            continue
        any_found = True
        change = parsed[3]
        tally[change] = tally.get(change, 0) + c
    if not any_found:
        return _no_data_figure(path, title, "No single-residue missense substitutions to classify")
    colors = {"Conservative (same class)": "#2e7d32", "Radical (class change)": "#c62828",
              "Unclassified": "#90a4ae"}
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ks = list(tally.keys())
    ax.bar(ks, [tally[k] for k in ks], color=[colors.get(k, "#607d8b") for k in ks], edgecolor="white")
    ax.set_ylabel("Isolate-mutation events")
    ax.set_title(title)
    plt.setp(ax.get_xticklabels(), rotation=15, ha="right")
    _safe_savefig(fig, path)

def _chart_cooccurrence_heatmap(all_mutation_data: List[dict], counts: Dict[str, int], path: str,
                                 title: str, top_n: int = 15):
    import numpy as np
    mut_order = [m for m, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]]
    if len(mut_order) < 2 or not all_mutation_data:
        return _no_data_figure(path, title, "Fewer than 2 mutations — no co-occurrence to show")
    n = len(mut_order)
    mat = np.zeros((n, n))
    for d in all_mutation_data:
        present = [m for m in mut_order if m in set(d.get("mutations", []))]
        for a in present:
            for b in present:
                mat[mut_order.index(a)][mut_order.index(b)] += 1
    fig, ax = plt.subplots(figsize=(max(6, 0.42 * n), max(6, 0.42 * n)))
    im = ax.imshow(mat, cmap="Purples")
    ax.set_xticks(range(n)); ax.set_xticklabels(mut_order, rotation=90, fontsize=7)
    ax.set_yticks(range(n)); ax.set_yticklabels(mut_order, fontsize=7)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Co-occurring isolates")
    ax.grid(False)
    _safe_savefig(fig, path)

def _chart_cluster_or_pca(all_mutation_data: List[dict], counts: Dict[str, int], path: str,
                           title: str, kind: str = "pca", top_n: int = 25):
    import numpy as np
    mut_order = [m for m, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]]
    if len(mut_order) < 2 or len(all_mutation_data) < 3:
        return _no_data_figure(path, title, "Not enough isolates/mutations for clustering")
    X = np.array([[1.0 if m in set(d.get("mutations", [])) else 0.0 for m in mut_order]
                   for d in all_mutation_data])
    if kind == "pca":
        Xc = X - X.mean(axis=0, keepdims=True)
        try:
            U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
            coords = U[:, :2] * S[:2]
        except Exception:
            return _no_data_figure(path, title, "PCA could not be computed for this dataset")
        fig, ax = plt.subplots(figsize=(6.5, 6))
        ax.scatter(coords[:, 0], coords[:, 1], s=45, color=PUB_PALETTE["primary"], alpha=0.75, edgecolor="white")
        ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
        ax.set_title(title)
        _safe_savefig(fig, path)
    else:
        try:
            from scipy.cluster.hierarchy import dendrogram, linkage
            Z = linkage(X, method="average", metric="hamming")
            fig, ax = plt.subplots(figsize=(max(8, 0.25 * len(all_mutation_data)), 5))
            dendrogram(Z, labels=[short_isolate_label(d["isolate"], 16) for d in all_mutation_data],
                       leaf_rotation=90, leaf_font_size=6, ax=ax, color_threshold=0)
            ax.set_ylabel("Hamming distance (mutation profile)")
            ax.set_title(title)
            _safe_savefig(fig, path)
        except Exception:
            return _no_data_figure(path, title, "SciPy not available — dendrogram skipped")

def _chart_manifest(gene_name, counts, total, all_data, all_mutation_data, type_map,
                     position_map, thresholds, label):
    thr_list = thresholds or [25, 30, 35, 40, 50]
    m = []
    m.append(("01_frequency_bar_count.png", lambda p: _chart_bar(counts, total, p,
              f"{gene_name} — mutation frequency (isolate count) [{label}]", mode="count", orientation="h")))
    m.append(("02_frequency_bar_percent.png", lambda p: _chart_bar(counts, total, p,
              f"{gene_name} — mutation prevalence (%) [{label}]", mode="percent", orientation="h")))
    m.append(("03_top10_bar_vertical.png", lambda p: _chart_bar(counts, total, p,
              f"{gene_name} — top 10 most prevalent mutations [{label}]", mode="percent", orientation="v", top_n=10)))
    m.append(("04_top20_bar_vertical.png", lambda p: _chart_bar(counts, total, p,
              f"{gene_name} — top 20 most prevalent mutations [{label}]", mode="percent", orientation="v", top_n=20)))
    m.append(("05_top30_bar_vertical.png", lambda p: _chart_bar(counts, total, p,
              f"{gene_name} — top 30 most prevalent mutations [{label}]", mode="percent", orientation="v", top_n=30)))
    m.append(("06_lollipop_count.png", lambda p: _chart_lollipop(counts, total, p,
              f"{gene_name} — mutation frequency lollipop plot (count) [{label}]", mode="count")))
    m.append(("07_lollipop_percent.png", lambda p: _chart_lollipop(counts, total, p,
              f"{gene_name} — mutation prevalence lollipop plot (%) [{label}]", mode="percent")))
    m.append(("08_pie_top8_mutations.png", lambda p: _chart_pie_top(counts, p,
              f"{gene_name} — share of isolate-mutation events, top 8 [{label}]")))
    m.append(("09_donut_prevalence_categories.png", lambda p: _chart_donut_categories(counts, total, p,
              f"{gene_name} — mutations by prevalence category [{label}]")))
    m.append(("10_prevalence_category_bar.png", lambda p: _chart_category_bar(counts, total, p,
              f"{gene_name} — count of mutations per prevalence category [{label}]")))
    m.append(("11_isolate_mutation_heatmap_top20.png", lambda p: _chart_isolate_mutation_heatmap(
              all_mutation_data, counts, p, f"{gene_name} — isolate × mutation presence heatmap, top 20 [{label}]", top_n=20)))
    m.append(("12_isolate_mutation_heatmap_top40.png", lambda p: _chart_isolate_mutation_heatmap(
              all_mutation_data, counts, p, f"{gene_name} — isolate × mutation presence heatmap, top 40 [{label}]", top_n=40)))
    m.append(("13_mutation_type_pie.png", lambda p: _chart_type_distribution(counts, type_map, p,
              f"{gene_name} — mutation type distribution [{label}]", kind="pie")))
    m.append(("14_mutation_type_bar.png", lambda p: _chart_type_distribution(counts, type_map, p,
              f"{gene_name} — mutation type distribution [{label}]", kind="bar")))
    m.append(("15_protein_position_lollipop.png", lambda p: _chart_position_lollipop(counts, position_map, p,
              f"{gene_name} — mutation needle/lollipop plot by gene position [{label}]")))
    m.append(("16_position_density_histogram.png", lambda p: _chart_position_histogram(counts, position_map, p,
              f"{gene_name} — mutation density along gene position [{label}]")))
    m.append(("17_cumulative_isolate_coverage_curve.png", lambda p: _chart_cumulative_curve(counts, total, p,
              f"{gene_name} — cumulative isolate coverage by ranked mutations [{label}]")))
    m.append(("18_rank_frequency_loglog.png", lambda p: _chart_rank_frequency(counts, p,
              f"{gene_name} — mutation rank–frequency distribution (log–log) [{label}]")))
    m.append(("19_prevalence_wilson_ci.png", lambda p: _chart_ci_bar(counts, total, p,
              f"{gene_name} — mutation prevalence with 95% Wilson CI, top 15 [{label}]")))
    m.append(("20_isolates_with_vs_without_mutation.png", lambda p: _chart_status_pie(all_data, p,
              f"{gene_name} — isolates with vs without a detected mutation [{label}]", mode="mutation_presence")))
    m.append(("21_mutation_burden_per_isolate_hist.png", lambda p: _chart_burden_hist(all_mutation_data, p,
              f"{gene_name} — mutation burden per isolate (histogram) [{label}]", kind="hist")))
    m.append(("22_mutation_burden_per_isolate_box.png", lambda p: _chart_burden_hist(all_mutation_data, p,
              f"{gene_name} — mutation burden per isolate (boxplot) [{label}]", kind="box")))
    m.append(("23_top_mutations_radar.png", lambda p: _chart_radar_top(counts, total, p,
              f"{gene_name} — top mutation prevalence radar [{label}]")))
    m.append(("24_cumulative_waterfall_top12.png", lambda p: _chart_waterfall(counts, total, p,
              f"{gene_name} — cumulative isolate-coverage waterfall, top 12 [{label}]")))
    m.append(("25_dotmatrix_top15_isolates.png", lambda p: _chart_dotmatrix(all_mutation_data, counts, p,
              f"{gene_name} — presence dot-matrix, top 15 mutations [{label}]")))
    m.append(("26_qc_pident_vs_qcov_scatter.png", lambda p: _chart_qc_scatter(all_data, p,
              f"{gene_name} — reference identity vs. query coverage by isolate [{label}]")))
    m.append(("27_isolate_pident_bar_sorted.png", lambda p: _chart_isolate_metric_bar(all_data, "pident", p,
              f"{gene_name} — percent identity to reference, isolates sorted [{label}]")))
    m.append(("28_isolate_qcov_bar_sorted.png", lambda p: _chart_isolate_metric_bar(all_data, "qcov", p,
              f"{gene_name} — reference query coverage, isolates sorted [{label}]")))
    m.append(("29_gene_status_pie.png", lambda p: _chart_status_pie(all_data, p,
              f"{gene_name} — gene detection status across isolates [{label}]", mode="gene_status")))
    for thr in thr_list[:5]:
        m.append((f"30_common_mutations_ge{thr}pct.png", (lambda p, t=thr: _chart_threshold_bar(
                  counts, total, t, p, f"{gene_name} — mutations at ≥{t}% prevalence [{label}]"))))
    m.append(("35_threshold_trend_line.png", lambda p: _chart_threshold_trend(counts, total, thr_list, p,
              f"{gene_name} — number of common mutations vs. prevalence threshold [{label}]")))
    m.append(("36_aa_substitution_matrix_heatmap.png", lambda p: _chart_aa_substitution_heatmap(counts, p,
              f"{gene_name} — amino-acid substitution matrix (reference→variant) [{label}]")))
    m.append(("37_physicochemical_change_bar.png", lambda p: _chart_physicochemical_bar(counts, p,
              f"{gene_name} — conservative vs. radical substitutions [{label}]")))
    m.append(("38_within_gene_cooccurrence_heatmap.png", lambda p: _chart_cooccurrence_heatmap(
              all_mutation_data, counts, p, f"{gene_name} — within-gene mutation co-occurrence, top 15 [{label}]")))
    m.append(("39_isolate_similarity_dendrogram.png", lambda p: _chart_cluster_or_pca(
              all_mutation_data, counts, p, f"{gene_name} — isolate clustering by mutation profile [{label}]", kind="dendro")))
    m.append(("40_isolate_pca_scatter.png", lambda p: _chart_cluster_or_pca(
              all_mutation_data, counts, p, f"{gene_name} — isolate ordination (PCA) by mutation profile [{label}]", kind="pca")))
    m.append(("41_summary_panel_figure.png", lambda p: _chart_summary_panel(
              gene_name, counts, total, all_data, all_mutation_data, p, f"{gene_name} — summary panel [{label}]")))
    m.append(("42_summary_statistics_table.png", lambda p: _chart_summary_table(
              gene_name, counts, total, all_data, p, f"{gene_name} — summary statistics [{label}]", label)))
    m.append(("43_top_mutations_labeled_bar.png", lambda p: _chart_bar(counts, total, p,
              f"{gene_name} — top 15 mutations, isolate counts labelled [{label}]", mode="count", orientation="h", top_n=15,
              color=PUB_PALETTE["curated"] if "AMR" in label else PUB_PALETTE["primary"])))
    m.append(("44_pareto_chart_top15.png", lambda p: _chart_pareto(counts, total, p,
              f"{gene_name} — Pareto chart of mutation events, top 15 [{label}]")))
    m.append(("45_prevalence_category_stacked.png", lambda p: _chart_category_bar(counts, total, p,
              f"{gene_name} — prevalence category breakdown (alt. view) [{label}]")))
    return m

def _run_battery(gene_name, mutation_counts, total_isolates, all_data, all_mutation_data,
                  curated_set, out_dir, thresholds, relevant_only: bool):
    all_counts, relevant_counts = _split_relevant(gene_name, mutation_counts, curated_set)
    counts = relevant_counts if relevant_only else all_counts
    label = "NG AMR-relevant only" if relevant_only else "all detected mutations"

    if relevant_only:
        filtered_amd = []
        for d in all_mutation_data or []:
            keep = [x for x in d.get("mutations", []) if x in counts]
            filtered_amd.append({"isolate": d["isolate"], "mutations": keep})
    else:
        filtered_amd = all_mutation_data or []

    position_map, type_map = _notation_metadata(all_data)
    os.makedirs(out_dir, exist_ok=True)
    manifest = _chart_manifest(gene_name, counts, total_isolates, all_data, filtered_amd,
                                type_map, position_map, thresholds, label)
    ok, failed = 0, 0
    for fname, fn in manifest:
        out_path = os.path.join(out_dir, fname)
        try:
            fn(out_path)
            ok += 1
        except Exception as e:
            failed += 1
            logger.warning(f"    chart failed [{gene_name}] {fname}: {e}")
            try:
                _no_data_figure(out_path, f"{gene_name} — {fname}", f"Chart generation error: {e}")
            except Exception as _exc:
                logger.debug("Ignored non-critical exception: %s", _exc)
                pass
    return ok, failed, len(manifest)

def generate_gene_wise_images(gene_results: List[dict], output_root: str,
                               curated_master: Optional[dict] = None,
                               thresholds: Optional[List[int]] = None) -> str:
    if not MATPLOTLIB_OK:
        logger.info("  [genes_wise_images] matplotlib not available — skipped")
        return ""
    _pub_rc()
    root = os.path.join(output_root, "genes_wise_images")
    os.makedirs(root, exist_ok=True)
    logger.info("\n" + "=" * 70)
    logger.info("GENE-WISE PUBLICATION IMAGE GALLERIES")
    logger.info("=" * 70)
    for res in gene_results:
        if not res:
            continue
        gene_name = res.get("gene")
        counts = res.get("mutation_counts") or {}
        total = res.get("total_isolates") or 0
        all_data = res.get("all_data") or []
        all_mutation_data = res.get("all_mutation_data") or []
        curated_set = res.get("curated_set") or set()
        if curated_master:
            curated_set = set(curated_set) | set(curated_master.get(gene_name, []))
        folder_name = GENE_FOLDER_NAMES.get(gene_name, gene_name)
        gene_root = os.path.join(root, folder_name)
        all_dir = os.path.join(gene_root, "All_Mutations")
        rel_dir = os.path.join(gene_root, "NG_AMR_Relevant")
        ok1, fail1, n1 = _run_battery(gene_name, counts, total, all_data, all_mutation_data,
                                       curated_set, all_dir, thresholds, relevant_only=False)
        ok2, fail2, n2 = _run_battery(gene_name, counts, total, all_data, all_mutation_data,
                                       curated_set, rel_dir, thresholds, relevant_only=True)
        logger.info(f"  {gene_name:16s}  All_Mutations: {ok1}/{n1} images"
              f"   NG_AMR_Relevant: {ok2}/{n2} images"
              + (f"   (failures: {fail1 + fail2})" if (fail1 + fail2) else ""))
    logger.info(f"\nGene-wise galleries written under: {root}")
    return root

def generate_master_publication_figures(gene_results: List[dict], output_root: str,
                                         curated_master: Optional[dict] = None) -> str:
    if not MATPLOTLIB_OK:
        logger.info("  [Publication_Figures] matplotlib not available — skipped")
        return ""
    _pub_rc()
    root = os.path.join(output_root, "Publication_Figures")
    all_dir = os.path.join(root, "All_Mutations")
    rel_dir = os.path.join(root, "NG_AMR_Relevant")
    os.makedirs(all_dir, exist_ok=True)
    os.makedirs(rel_dir, exist_ok=True)

    valid = [r for r in gene_results if r and r.get("mutation_counts")]

    def _gene_totals(relevant_only):
        rows = []
        for r in valid:
            gene = r["gene"]
            cset = set(r.get("curated_set") or set())
            if curated_master:
                cset |= set(curated_master.get(gene, []))
            all_c, rel_c = _split_relevant(gene, r["mutation_counts"], cset)
            c = rel_c if relevant_only else all_c
            rows.append((gene, c, r.get("total_isolates") or 0))
        return rows

    def _build(relevant_only, out_dir, label):
        rows = _gene_totals(relevant_only)
        gene_names = [g for g, _, _ in rows]
        n_distinct = [len(c) for _, c, _ in rows]
        n_events = [sum(c.values()) for _, c, _ in rows]
        mean_burden = [ (sum(c.values()) / t) if t else 0 for _, c, t in rows]

        fig, ax = plt.subplots(figsize=(max(9, 0.4 * len(gene_names)), 5))
        ax.bar(gene_names, n_distinct, color=PUB_PALETTE["primary"], edgecolor="white")
        ax.set_ylabel("Distinct mutations detected")
        ax.set_title(f"Distinct mutations per gene — {label}")
        plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
        _safe_savefig(fig, os.path.join(out_dir, "01_distinct_mutations_per_gene.png"))

        fig, ax = plt.subplots(figsize=(max(9, 0.4 * len(gene_names)), 5))
        ax.bar(gene_names, n_events, color=PUB_PALETTE["accent"], edgecolor="white")
        ax.set_ylabel("Total isolate-mutation events")
        ax.set_title(f"Mutation event burden per gene — {label}")
        plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
        _safe_savefig(fig, os.path.join(out_dir, "02_mutation_events_per_gene.png"))

        fig, ax = plt.subplots(figsize=(max(9, 0.4 * len(gene_names)), 5))
        ax.bar(gene_names, mean_burden, color=PUB_PALETTE["curated"], edgecolor="white")
        ax.set_ylabel("Mean mutations per isolate")
        ax.set_title(f"Mean mutation burden per isolate, by gene — {label}")
        plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
        _safe_savefig(fig, os.path.join(out_dir, "03_mean_burden_per_isolate_per_gene.png"))

        pooled = {}
        gene_of = {}
        for gene, c, _ in rows:
            for m, cnt in c.items():
                key = f"{gene}:{m}"
                pooled[key] = cnt
                gene_of[key] = gene
        top = sorted(pooled.items(), key=lambda kv: -kv[1])[:25]
        if top:
            cmap = plt.get_cmap("tab20")
            genes_uniq = sorted(set(gene_of[k] for k, _ in top))
            gcolor = {g: cmap(i % 20) for i, g in enumerate(genes_uniq)}
            fig, ax = plt.subplots(figsize=(9, max(4, 0.32 * len(top) + 1)))
            labels = [k for k, _ in top]
            vals = [v for _, v in top]
            colors = [gcolor[gene_of[k]] for k, _ in top]
            ax.barh(range(len(labels)), vals, color=colors)
            ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7)
            ax.invert_yaxis()
            ax.set_xlabel("Number of isolates")
            ax.set_title(f"Top 25 mutations across all genes — {label}")
            from matplotlib.patches import Patch
            handles = [Patch(color=gcolor[g], label=g) for g in genes_uniq]
            ax.legend(handles=handles, loc="lower right", fontsize=7, ncol=2, frameon=False)
            _safe_savefig(fig, os.path.join(out_dir, "04_top25_mutations_all_genes.png"))
        else:
            _no_data_figure(os.path.join(out_dir, "04_top25_mutations_all_genes.png"),
                             f"Top 25 mutations across all genes — {label}")

        pct_with = []
        for gene, c, t in rows:
            n_isolates = t
            if n_isolates == 0:
                pct_with.append(0); continue
            src = next((r for r in valid if r["gene"] == gene), None)
            amd = src.get("all_mutation_data") or [] if src else []
            keep_set = set(c.keys())
            n_with = sum(1 for d in amd if set(d.get("mutations", [])) & keep_set)
            pct_with.append(100.0 * n_with / n_isolates if n_isolates else 0)
        fig, ax = plt.subplots(figsize=(max(9, 0.4 * len(gene_names)), 5))
        ax.bar(gene_names, pct_with, color=PUB_PALETTE["warn"], edgecolor="white")
        ax.set_ylabel("Isolates with ≥1 qualifying mutation (%)")
        ax.set_ylim(0, 105)
        ax.set_title(f"Isolate-level mutation prevalence per gene — {label}")
        plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
        _safe_savefig(fig, os.path.join(out_dir, "05_pct_isolates_with_mutation_per_gene.png"))

        found, disrupted, missing = [], [], []
        for r in valid:
            ad = r.get("all_data") or []
            f = sum(1 for d in ad if d.get("gene_found"))
            dch = sum(1 for d in ad if d.get("gene_found") and d.get("is_disrupted"))
            mi = sum(1 for d in ad if not d.get("gene_found"))
            found.append(max(0, f - dch)); disrupted.append(dch); missing.append(mi)
        gnames_all = [r["gene"] for r in valid]
        fig, ax = plt.subplots(figsize=(max(9, 0.4 * len(gnames_all)), 5))
        ax.bar(gnames_all, found, label="Intact", color="#2e7d32")
        ax.bar(gnames_all, disrupted, bottom=found, label="Disrupted", color="#ef6c00")
        bottom2 = [f + d for f, d in zip(found, disrupted)]
        ax.bar(gnames_all, missing, bottom=bottom2, label="Missing", color="#b0bec5")
        ax.set_ylabel("Isolates")
        ax.set_title("Gene detection status across all genes")
        ax.legend(frameon=False, fontsize=8)
        plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
        _safe_savefig(fig, os.path.join(out_dir, "06_gene_status_stacked_bar.png"))

        import numpy as np
        max_per_gene = 6
        col_labels, row_labels = [], gene_names
        cell_data = []
        for gene, c, t in rows:
            top6 = sorted(c.items(), key=lambda kv: -kv[1])[:max_per_gene]
            cell_data.append(top6)
        max_cols = max((len(x) for x in cell_data), default=0)
        if max_cols:
            mat = np.zeros((len(rows), max_cols))
            annot = [["" for _ in range(max_cols)] for _ in range(len(rows))]
            for i, ((gene, c, t), top6) in enumerate(zip(rows, cell_data)):
                for j, (m, cnt) in enumerate(top6):
                    mat[i, j] = 100.0 * cnt / t if t else 0
                    annot[i][j] = m
            fig, ax = plt.subplots(figsize=(max(8, 0.9 * max_cols), max(6, 0.35 * len(rows) + 1)))
            im = ax.imshow(mat, cmap="YlOrRd", aspect="auto")
            ax.set_yticks(range(len(row_labels))); ax.set_yticklabels(row_labels, fontsize=8)
            ax.set_xticks(range(max_cols)); ax.set_xticklabels([f"Rank {i+1}" for i in range(max_cols)], fontsize=8)
            for i in range(len(rows)):
                for j in range(max_cols):
                    if annot[i][j]:
                        ax.text(j, i, _fig_short_label(annot[i][j], 14), ha="center", va="center", fontsize=7)
            fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Prevalence (%)")
            ax.set_title(f"Top-ranked mutation per gene — {label}")
            ax.grid(False)
            _safe_savefig(fig, os.path.join(out_dir, "07_gene_by_top_mutation_heatmap.png"))
        else:
            _no_data_figure(os.path.join(out_dir, "07_gene_by_top_mutation_heatmap.png"),
                             f"Top-ranked mutation per gene — {label}")

        gene_share = {g: sum(c.values()) for g, c, _ in rows}
        gene_share = {g: v for g, v in gene_share.items() if v > 0}
        if gene_share:
            fig, ax = plt.subplots(figsize=(7.5, 7.5))
            cmap = plt.get_cmap("tab20")
            colors = [cmap(i % 20) for i in range(len(gene_share))]
            ax.pie(gene_share.values(), labels=gene_share.keys(), autopct="%1.0f%%",
                   colors=colors, startangle=90, textprops={"fontsize": 7})
            ax.set_title(f"Share of mutation events by gene — {label}")
            _safe_savefig(fig, os.path.join(out_dir, "08_gene_share_pie.png"))
        else:
            _no_data_figure(os.path.join(out_dir, "08_gene_share_pie.png"), f"Share of mutation events by gene — {label}")

        rows_tbl = [["Gene", "Isolates", "Distinct mut.", "Events", "% isolates w/ mut."]]
        for (gene, c, t), pw in zip(rows, pct_with):
            rows_tbl.append([gene, str(t), str(len(c)), str(sum(c.values())), f"{pw:.1f}%"])
        fig, ax = plt.subplots(figsize=(9, 0.4 * len(rows_tbl) + 1))
        ax.axis("off")
        tbl = ax.table(cellText=rows_tbl[1:], colLabels=rows_tbl[0], cellLoc="center", loc="center")
        tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.4)
        for (r, c_), cell in tbl.get_celld().items():
            if r == 0:
                cell.set_facecolor(PUB_PALETTE["primary"]); cell.set_text_props(color="white", fontweight="bold")
        ax.set_title(f"Cross-gene summary table — {label}", pad=16)
        _safe_savefig(fig, os.path.join(out_dir, "09_cross_gene_summary_table.png"))

    _build(False, all_dir, "all detected mutations, whole dataset")
    _build(True, rel_dir, "NG AMR-relevant mutations only, whole dataset")

    logger.info(f"\nMaster publication figures written under: {root}")
    return root

def generate_thesis_publication_package(gene_results: List[dict], output_root: str,
                                         curated_master: Optional[dict] = None,
                                         thresholds: Optional[List[int]] = None) -> None:
    logger.info("\n" + "#" * 70)
    logger.info("# THESIS PUBLICATION IMAGE PACKAGE")
    logger.info("#" * 70)
    generate_master_publication_figures(gene_results, output_root, curated_master)
    generate_gene_wise_images(gene_results, output_root, curated_master, thresholds)
    logger.info("\nDone. New folders created (existing pipeline outputs untouched):")
    logger.info(f"  {os.path.join(output_root, 'Publication_Figures', 'All_Mutations')}")
    logger.info(f"  {os.path.join(output_root, 'Publication_Figures', 'NG_AMR_Relevant')}")
    logger.info(f"  {os.path.join(output_root, 'genes_wise_images', '<gene>', 'All_Mutations')}")
    logger.info(f"  {os.path.join(output_root, 'genes_wise_images', '<gene>', 'NG_AMR_Relevant')}")

def _ext_no_data(path: str, title: str, note: str = "No qualifying data for this chart"):
    if not MATPLOTLIB_OK:
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axis("off")
    ax.set_title(title, pad=12)
    ax.text(0.5, 0.5, note, ha="center", va="center", fontsize=11,
            color="#757575", transform=ax.transAxes)
    _safe_savefig(fig, path)

def _ext_items(counts, top_n=None, sort_asc=False):
    items = sorted(counts.items(), key=lambda kv: (kv[1], kv[0]) if sort_asc else (-kv[1], kv[0]))
    if top_n:
        items = items[:top_n]
    return items

def _render_ranked_bar_chart(counts, total, path, title, top_n=25, orientation="v",
                              mode="count", color="#1565c0", ylabel_count="Number of isolates",
                              ylabel_pct="Isolates carrying mutation (%)"):
    items = _ext_items(counts, top_n)
    if not items or (mode == "pct" and not total):
        return _ext_no_data(path, title)
    labels = [k for k, _ in items]
    if mode == "pct":
        vals = [100.0 * v / total for _, v in items]
        axis_label = ylabel_pct
    else:
        vals = [v for _, v in items]
        axis_label = ylabel_count

    if orientation == "h":
        fig, ax = plt.subplots(figsize=(8, max(3.2, 0.28 * len(labels) + 1)))
        ax.barh(range(len(labels)), vals, color=color, edgecolor="white")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=7.5)
        ax.invert_yaxis()
        ax.set_xlabel(axis_label)
    else:
        fig, ax = plt.subplots(figsize=(max(7, 0.35 * len(labels)), 5))
        ax.bar(range(len(labels)), vals, color=color, edgecolor="white")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
        ax.set_ylabel(axis_label)
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_mut_freq_bar(counts, total, path, title, top_n=25):
    return _render_ranked_bar_chart(counts, total, path, title, top_n=top_n, orientation="v", mode="count", color="#1565c0")

def _chart_mut_freq_hbar(counts, total, path, title, top_n=25):
    return _render_ranked_bar_chart(counts, total, path, title, top_n=top_n, orientation="h", mode="count", color="#2e7d32")

def _chart_mut_pct_bar(counts, total, path, title, top_n=25):
    return _render_ranked_bar_chart(counts, total, path, title, top_n=top_n, orientation="v", mode="pct", color="#e65100")

def _chart_mut_pct_hbar(counts, total, path, title, top_n=25):
    return _render_ranked_bar_chart(counts, total, path, title, top_n=top_n, orientation="h", mode="pct", color="#6a1b9a")

def _chart_sorted_freq(counts, total, path, title, top_n=30):
    items = _ext_items(counts, top_n, sort_asc=True)
    if not items:
        return _ext_no_data(path, title)
    labels, vals = zip(*items)
    fig, ax = plt.subplots(figsize=(8, max(3.2, 0.28 * len(labels) + 1)))
    ax.barh(range(len(labels)), vals, color="#00897b", edgecolor="white")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("Number of isolates")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_dot_freq(counts, total, path, title, top_n=25):
    items = _ext_items(counts, top_n)
    if not items:
        return _ext_no_data(path, title)
    labels, vals = zip(*items)
    fig, ax = plt.subplots(figsize=(8, max(3.2, 0.28 * len(labels) + 1)))
    ax.plot(vals, range(len(labels)), "o", color="#1565c0", markersize=7)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Number of isolates")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_lollipop(counts, total, path, title, top_n=25, mode="count"):
    items = _ext_items(counts, top_n, sort_asc=True)
    if not items:
        return _ext_no_data(path, title)
    labels = [k for k, _ in items]
    vals = [100.0 * v / total if (mode == "percent" and total) else v for _, v in items]
    fig, ax = plt.subplots(figsize=(8, max(3.2, 0.28 * len(labels) + 1)))
    ax.hlines(range(len(labels)), 0, vals, color="#455a64", linewidth=1.3)
    ax.plot(vals, range(len(labels)), "o", color="#c62828", markersize=6)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("Isolates (%)" if mode == "percent" else "Number of isolates")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_stem(counts, total, path, title, top_n=20, mode="count"):
    items = _ext_items(counts, top_n)
    if not items:
        return _ext_no_data(path, title)
    labels = [k for k, _ in items]
    vals = [100.0 * v / total if (mode == "percent" and total) else v for _, v in items]
    fig, ax = plt.subplots(figsize=(max(7, 0.4 * len(labels)), 4.8))
    markerline, stemlines, baseline = ax.stem(range(len(labels)), vals, linefmt="#6a1b9a",
                                               markerfmt="o", basefmt=" ")
    plt.setp(markerline, color="#6a1b9a", markersize=5)
    plt.setp(stemlines, color="#6a1b9a", linewidth=1.2)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Isolates (%)" if mode == "percent" else "Count")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_needle(counts, total, path, title, top_n=20):
    items = _ext_items(counts, top_n)
    if not items:
        return _ext_no_data(path, title)
    labels, vals = zip(*items)
    fig, ax = plt.subplots(figsize=(max(7, 0.4 * len(labels)), 4.8))
    ax.vlines(range(len(labels)), 0, vals, color="#1565c0", linewidth=1.5)
    ax.scatter(range(len(labels)), vals, s=28, color="#c62828", zorder=3)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Number of isolates")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_hist_counts(counts, path, title):
    if not counts:
        return _ext_no_data(path, title)
    vals = list(counts.values())
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(vals, bins=min(20, max(5, len(set(vals)))), color="#1565c0", edgecolor="white")
    ax.set_xlabel("Isolate count per mutation")
    ax.set_ylabel("Number of distinct mutations")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_density_counts(counts, path, title):
    if not counts:
        return _ext_no_data(path, title)
    vals = list(counts.values())
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(vals, bins=min(20, max(5, len(set(vals)))), density=True,
            color="#bbdefb", edgecolor="#1565c0", histtype="stepfilled", alpha=0.85)
    ax.set_xlabel("Isolate count per mutation")
    ax.set_ylabel("Density")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_rug(counts, path, title):
    if not counts:
        return _ext_no_data(path, title)
    vals = sorted(counts.values())
    fig, ax = plt.subplots(figsize=(8, 2.8))
    ax.eventplot(vals, orientation="horizontal", colors="#c62828", lineoffsets=0, linelengths=0.8)
    ax.set_yticks([])
    ax.set_xlabel("Isolate count")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_strip(counts, path, title, top_n=40):
    items = _ext_items(counts, top_n)
    if not items:
        return _ext_no_data(path, title)
    labels, vals = zip(*items)
    fig, ax = plt.subplots(figsize=(8, max(3, 0.22 * len(labels) + 1)))
    ax.scatter(vals, range(len(labels)), s=18, color="#1565c0", alpha=0.8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("Number of isolates")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_bubble(counts, total, path, title, top_n=25):
    items = _ext_items(counts, top_n)
    if not items:
        return _ext_no_data(path, title)
    labels, vals = zip(*items)
    sizes = [max(30, 12 * v) for v in vals]
    fig, ax = plt.subplots(figsize=(max(7, 0.4 * len(labels)), 5))
    ax.scatter(range(len(labels)), vals, s=sizes, c=vals, cmap="YlOrRd",
               alpha=0.75, edgecolors="white")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Number of isolates")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_sized_point(counts, total, path, title, top_n=25):
    items = _ext_items(counts, top_n)
    if not items or not total:
        return _ext_no_data(path, title)
    labels = [k for k, _ in items]
    vals = [100.0 * v / total for _, v in items]
    sizes = [max(40, 4 * v) for v in vals]
    fig, ax = plt.subplots(figsize=(8, max(3.2, 0.28 * len(labels) + 1)))
    ax.scatter(vals, range(len(labels)), s=sizes, c="#0277bd", alpha=0.8, edgecolors="white")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Prevalence (%)")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_presence_heatmap(all_mutation_data, counts, path, title, top_n=20, cmap="Blues"):
    if not counts or not all_mutation_data:
        return _ext_no_data(path, title)
    import numpy as np
    mut_order = [m for m, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]]
    isolates = [d["isolate"] for d in all_mutation_data]
    if not mut_order or not isolates:
        return _ext_no_data(path, title)
    mat = np.zeros((len(isolates), len(mut_order)))
    for i, d in enumerate(all_mutation_data):
        present = set(d.get("mutations", []))
        for j, m in enumerate(mut_order):
            mat[i, j] = 1.0 if m in present else 0.0
    fig_h = max(4, min(22, 0.16 * len(isolates) + 1.2))
    fig, ax = plt.subplots(figsize=(max(7.5, 0.35 * len(mut_order)), fig_h))
    ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(mut_order)))
    ax.set_xticklabels(mut_order, rotation=90, fontsize=6.5)
    if len(isolates) <= 50:
        ax.set_yticks(range(len(isolates)))
        ax.set_yticklabels([short_isolate_label(x, 20) for x in isolates], fontsize=5.5)
    else:
        ax.set_yticks([])
        ax.set_ylabel(f"{len(isolates)} isolates")
    ax.set_title(title)
    ax.grid(False)
    _safe_savefig(fig, path)

def _chart_burden_heatmap(all_mutation_data, counts, path, title, top_n=20):
    if not counts or not all_mutation_data:
        return _ext_no_data(path, title)
    import numpy as np
    mut_order = [m for m, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]]
    isolates = [d["isolate"] for d in all_mutation_data]
    mat = np.zeros((len(isolates), len(mut_order)))
    for i, d in enumerate(all_mutation_data):
        present = set(d.get("mutations", []))
        for j, m in enumerate(mut_order):
            mat[i, j] = counts.get(m, 0) if m in present else 0
    fig, ax = plt.subplots(figsize=(max(7.5, 0.35 * len(mut_order)), max(4, min(18, 0.14 * len(isolates) + 1))))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(mut_order)))
    ax.set_xticklabels(mut_order, rotation=90, fontsize=6.5)
    ax.set_yticks([])
    ax.set_ylabel(f"{len(isolates)} isolates")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Global count")
    ax.grid(False)
    _safe_savefig(fig, path)

def _chart_intensity_heatmap(all_mutation_data, counts, path, title, top_n=18):
    return _chart_presence_heatmap(all_mutation_data, counts, path, title, top_n=top_n, cmap="OrRd")

def _chart_binary_matrix(all_mutation_data, counts, path, title, top_n=25):
    return _chart_presence_heatmap(all_mutation_data, counts, path, title, top_n=top_n, cmap="Greys")

def _chart_dotmatrix(all_mutation_data, counts, path, title, top_n=15):
    if not counts or not all_mutation_data:
        return _ext_no_data(path, title)
    mut_order = [m for m, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]]
    isolates = [d["isolate"] for d in all_mutation_data][:55]
    fig, ax = plt.subplots(figsize=(max(7.5, 0.32 * len(mut_order)), max(3.8, 0.2 * len(isolates) + 1)))
    for i, d in enumerate(all_mutation_data[:55]):
        present = set(d.get("mutations", []))
        for j, m in enumerate(mut_order):
            if m in present:
                ax.scatter(j, i, s=40, color="#1565c0")
    ax.set_xticks(range(len(mut_order)))
    ax.set_xticklabels(mut_order, rotation=90, fontsize=6.5)
    ax.set_yticks(range(len(isolates)))
    ax.set_yticklabels([short_isolate_label(x, 18) for x in isolates], fontsize=5.5)
    ax.invert_yaxis()
    ax.set_title(title + (f" (first 55 of {len(all_mutation_data)} isolates)" if len(all_mutation_data)>55 else ""))
    ax.grid(True, linestyle=":", alpha=0.35)
    _safe_savefig(fig, path)

def _chart_type_pie(counts, type_map, path, title):
    if not counts:
        return _ext_no_data(path, title)
    tally = {}
    for m, c in counts.items():
        t = type_map.get(m, "substitution").replace("_", " ").title()
        tally[t] = tally.get(t, 0) + c
    colors_map = {"Substitution": "#1565c0", "Stop Gain": "#c62828", "Insertion": "#2e7d32",
                  "Deletion": "#ef6c00", "Frameshift Insertion": "#6a1b9a",
                  "Frameshift Deletion": "#4527a0"}
    cs = [colors_map.get(k, "#607d8b") for k in tally]
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.pie(tally.values(), labels=tally.keys(), autopct="%1.0f%%", colors=cs,
           startangle=90, wedgeprops=dict(edgecolor="white"))
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_type_bar(counts, type_map, path, title):
    if not counts:
        return _ext_no_data(path, title)
    tally = {}
    for m, c in counts.items():
        t = type_map.get(m, "substitution").replace("_", " ").title()
        tally[t] = tally.get(t, 0) + c
    colors_map = {"Substitution": "#1565c0", "Stop Gain": "#c62828", "Insertion": "#2e7d32",
                  "Deletion": "#ef6c00", "Frameshift Insertion": "#6a1b9a",
                  "Frameshift Deletion": "#4527a0"}
    ks = list(tally.keys())
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(ks, [tally[k] for k in ks], color=[colors_map.get(k, "#607d8b") for k in ks], edgecolor="white")
    ax.set_ylabel("Isolate-mutation events")
    plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_type_donut(counts, type_map, path, title):
    if not counts:
        return _ext_no_data(path, title)
    tally = {}
    for m, c in counts.items():
        t = type_map.get(m, "substitution").replace("_", " ").title()
        tally[t] = tally.get(t, 0) + c
    colors_map = {"Substitution": "#1565c0", "Stop Gain": "#c62828", "Insertion": "#2e7d32",
                  "Deletion": "#ef6c00", "Frameshift Insertion": "#6a1b9a",
                  "Frameshift Deletion": "#4527a0"}
    cs = [colors_map.get(k, "#607d8b") for k in tally]
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.pie(tally.values(), labels=tally.keys(), autopct="%1.0f%%", colors=cs,
           startangle=90, wedgeprops=dict(width=0.42, edgecolor="white"))
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_type_stacked_dummy(counts, type_map, path, title):
    if not counts:
        return _ext_no_data(path, title)
    tally = {}
    for m, c in counts.items():
        t = type_map.get(m, "substitution").replace("_", " ").title()
        tally[t] = tally.get(t, 0) + c
    total = sum(tally.values()) or 1
    fig, ax = plt.subplots(figsize=(8, 2.2))
    left = 0
    colors = ["#1565c0", "#c62828", "#2e7d32", "#ef6c00", "#6a1b9a", "#4527a0", "#607d8b"]
    for i, (k, v) in enumerate(tally.items()):
        w = 100.0 * v / total
        ax.barh(0, w, left=left, color=colors[i % len(colors)], edgecolor="white", height=0.6)
        if w >= 6:
            ax.text(left + w / 2, 0, f"{k}\n{w:.0f}%", ha="center", va="center", fontsize=7, color="white")
        left += w
    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.set_xlabel("Share of events (%)")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_position_scatter(counts, position_map, path, title):
    pts = [(position_map.get(m, 0), c, m) for m, c in counts.items() if position_map.get(m, 0)]
    if not pts:
        return _ext_no_data(path, title, "No positional mutations")
    pts.sort()
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    fig, ax = plt.subplots(figsize=(max(8, 0.12 * len(pts) + 4), 4.8))
    ax.scatter(xs, ys, s=36, c="#6a1b9a", alpha=0.8, edgecolors="white")
    ax.set_xlabel("Gene position")
    ax.set_ylabel("Isolate count")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_position_lollipop_local(counts, position_map, path, title):
    return _chart_position_lollipop(counts, position_map, path, title)

def _chart_position_hist_local(counts, position_map, path, title):
    return _chart_position_histogram(counts, position_map, path, title)

def _chart_snp_freq_by_position(counts, position_map, path, title):
    pts = [(position_map.get(m, 0), c) for m, c in counts.items() if position_map.get(m, 0)]
    if not pts:
        return _ext_no_data(path, title, "No positional mutations")
    pts.sort()
    xs, ys = zip(*pts)
    fig, ax = plt.subplots(figsize=(max(8, 0.1 * len(pts) + 4), 4.5))
    ax.plot(xs, ys, "-", color="#455a64", linewidth=1.0, alpha=0.6)
    ax.scatter(xs, ys, s=22, c="#c62828", zorder=3)
    ax.set_xlabel("Gene position")
    ax.set_ylabel("Isolate count")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_burden_hist(all_mutation_data, path, title):
    if not all_mutation_data:
        return _ext_no_data(path, title)
    burdens = [len(d.get("mutations", [])) for d in all_mutation_data]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    maxb = max(burdens) if burdens else 0
    ax.hist(burdens, bins=range(0, maxb + 2), color="#1565c0", edgecolor="white", align="left")
    ax.set_xlabel("Mutations per isolate")
    ax.set_ylabel("Number of isolates")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_burden_box(all_mutation_data, path, title):
    if not all_mutation_data:
        return _ext_no_data(path, title)
    burdens = [len(d.get("mutations", [])) for d in all_mutation_data]
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.boxplot(burdens, orientation="vertical", patch_artist=True,
               boxprops=dict(facecolor="#bbdefb", color="#1565c0"),
               medianprops=dict(color="#c62828"))
    ax.set_ylabel("Mutations per isolate")
    ax.set_xticks([])
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_burden_violin(all_mutation_data, path, title):
    if not all_mutation_data:
        return _ext_no_data(path, title)
    burdens = [len(d.get("mutations", [])) for d in all_mutation_data]
    fig, ax = plt.subplots(figsize=(5.5, 5))
    try:
        parts = ax.violinplot(burdens, showmeans=True, showmedians=True)
        for pc in parts.get("bodies", []):
            pc.set_facecolor("#c8e6c9")
            pc.set_edgecolor("#2e7d32")
            pc.set_alpha(0.85)
    except Exception:
        ax.hist(burdens, bins=max(5, min(15, (max(burdens) or 0) + 1)), color="#c8e6c9", edgecolor="#2e7d32")
    ax.set_ylabel("Mutations per isolate")
    ax.set_xticks([])
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_burden_ecdf(all_mutation_data, path, title):
    if not all_mutation_data:
        return _ext_no_data(path, title)
    burdens = sorted(len(d.get("mutations", [])) for d in all_mutation_data)
    n = len(burdens)
    y = [i / n for i in range(1, n + 1)]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.step(burdens, y, where="post", color="#6a1b9a", linewidth=1.6)
    ax.set_xlabel("Mutations per isolate")
    ax.set_ylabel("ECDF")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_burden_rank(all_mutation_data, path, title):
    if not all_mutation_data:
        return _ext_no_data(path, title)
    burdens = sorted((len(d.get("mutations", [])), d["isolate"]) for d in all_mutation_data)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(range(1, len(burdens) + 1), [b for b, _ in burdens], "o-", color="#e65100", markersize=3)
    ax.set_xlabel("Isolate rank (low → high burden)")
    ax.set_ylabel("Mutations per isolate")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_burden_waterfall(all_mutation_data, path, title, top_n=30):
    if not all_mutation_data:
        return _ext_no_data(path, title)
    rows = sorted(((len(d.get("mutations", [])), d["isolate"]) for d in all_mutation_data), reverse=True)[:top_n]
    labels = [short_isolate_label(iso, 14) for _, iso in rows]
    vals = [b for b, _ in rows]
    fig, ax = plt.subplots(figsize=(max(8, 0.28 * len(labels)), 4.8))
    ax.bar(range(len(labels)), vals, color="#1565c0", edgecolor="white")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=75, ha="right", fontsize=6)
    ax.set_ylabel("Mutations")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_ci_bar(counts, total, path, title, top_n=15):
    if not counts or not total:
        return _ext_no_data(path, title)
    items = sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]
    labels, mids, los, his = [], [], [], []
    for m, c in items:
        p, lo, hi = _wilson_ci(c, total)
        labels.append(m)
        mids.append(100 * p)
        los.append(100 * (p - lo))
        his.append(100 * (hi - p))
    fig, ax = plt.subplots(figsize=(8, max(3, 0.32 * len(labels) + 1)))
    y = range(len(labels))
    ax.errorbar(mids, y, xerr=[los, his], fmt="o", color="#1565c0",
                ecolor="#90caf9", elinewidth=1.4, capsize=3, markersize=5)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Prevalence (%) ± 95% Wilson CI")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_radar(counts, total, path, title, top_n=8):
    items = sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]
    if len(items) < 3:
        return _ext_no_data(path, title, "Fewer than 3 mutations")
    labels = [k for k, _ in items]
    values = [100.0 * v / total if total else 0 for _, v in items]
    values += values[:1]
    angles = [n / float(len(labels)) * 2 * math.pi for n in range(len(labels))]
    angles += angles[:1]
    fig, ax = plt.subplots(figsize=(6.5, 6.5), subplot_kw=dict(polar=True))
    ax.plot(angles, values, color="#6a1b9a", linewidth=1.5)
    ax.fill(angles, values, color="#e1bee7", alpha=0.45)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_title(title, pad=18)
    _safe_savefig(fig, path)

def _chart_waterfall(counts, total, path, title, top_n=12):
    if not counts or not total:
        return _ext_no_data(path, title)
    items = sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    running = 0
    cmap = plt.get_cmap("viridis")
    for i, (m, c) in enumerate(items):
        pct = 100.0 * c / total
        ax.bar(i, pct, bottom=running, color=cmap(i / max(1, len(items) - 1)),
               edgecolor="white", width=0.7)
        running += pct
    ax.set_xticks(range(len(items)))
    ax.set_xticklabels([m for m, _ in items], rotation=55, ha="right", fontsize=7)
    ax.set_ylabel("Cumulative isolate coverage (%)")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_pareto(counts, total, path, title, top_n=15):
    if not counts or not total:
        return _ext_no_data(path, title)
    items = sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]
    labels = [m for m, _ in items]
    values = [c for _, c in items]
    grand = sum(counts.values()) or 1
    cum, running = [], 0
    for c in values:
        running += c
        cum.append(100.0 * running / grand)
    fig, ax1 = plt.subplots(figsize=(8, 4.8))
    ax1.bar(range(len(labels)), values, color="#1565c0", edgecolor="white")
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, rotation=55, ha="right", fontsize=7)
    ax1.set_ylabel("Isolate count", color="#1565c0")
    ax2 = ax1.twinx()
    ax2.plot(range(len(labels)), cum, color="#c62828", marker="o", markersize=4)
    ax2.set_ylabel("Cumulative share of events (%)", color="#c62828")
    ax2.set_ylim(0, 105)
    ax2.grid(False)
    ax1.set_title(title)
    _safe_savefig(fig, path)

def _chart_cumulative(counts, total, path, title):
    if not counts or not total:
        return _ext_no_data(path, title)
    items = sorted(counts.items(), key=lambda kv: -kv[1])
    cum, running = [], 0
    for _, c in items:
        running = min(total, running + c)
        cum.append(100.0 * running / total)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(range(1, len(cum) + 1), cum, marker="o", markersize=3, color="#1565c0")
    ax.fill_between(range(1, len(cum) + 1), cum, alpha=0.18, color="#bbdefb")
    ax.set_xlabel("Ranked mutations (most prevalent first)")
    ax.set_ylabel("Cumulative isolate coverage (%)")
    ax.set_ylim(0, 105)
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_rank_freq(counts, path, title):
    if not counts:
        return _ext_no_data(path, title)
    values = sorted(counts.values(), reverse=True)
    ranks = range(1, len(values) + 1)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.plot(ranks, values, "o-", color="#6a1b9a", markersize=4, linewidth=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Mutation rank (log)")
    ax.set_ylabel("Isolate count (log)")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_category_bar(counts, total, path, title):
    if not counts or not total:
        return _ext_no_data(path, title)
    bins = {"Rare\n(<10%)": 0, "Low\n(10–25%)": 0, "Moderate\n(25–50%)": 0,
            "High\n(50–75%)": 0, "Majority\n(>75%)": 0}
    for _, c in counts.items():
        pct = 100.0 * c / total
        if pct < 10:
            bins["Rare\n(<10%)"] += 1
        elif pct < 25:
            bins["Low\n(10–25%)"] += 1
        elif pct < 50:
            bins["Moderate\n(25–50%)"] += 1
        elif pct < 75:
            bins["High\n(50–75%)"] += 1
        else:
            bins["Majority\n(>75%)"] += 1
    colors = ["#90a4ae", "#4fc3f7", "#ffb74d", "#ef5350", "#8e24aa"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(list(bins.keys()), list(bins.values()), color=colors, edgecolor="white")
    ax.set_ylabel("Number of distinct mutations")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_cooc(all_mutation_data, counts, path, title, top_n=12, cmap="Purples"):
    import numpy as np
    mut_order = [m for m, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]]
    if len(mut_order) < 2 or not all_mutation_data:
        return _ext_no_data(path, title, "Need ≥2 mutations")
    n = len(mut_order)
    mat = np.zeros((n, n))
    for d in all_mutation_data:
        present = [m for m in mut_order if m in set(d.get("mutations", []))]
        for a in present:
            for b in present:
                mat[mut_order.index(a)][mut_order.index(b)] += 1
    fig, ax = plt.subplots(figsize=(max(5.5, 0.4 * n), max(5.5, 0.4 * n)))
    im = ax.imshow(mat, cmap=cmap)
    ax.set_xticks(range(n))
    ax.set_xticklabels(mut_order, rotation=90, fontsize=6.5)
    ax.set_yticks(range(n))
    ax.set_yticklabels(mut_order, fontsize=6.5)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.grid(False)
    _safe_savefig(fig, path)

def _chart_pca(all_mutation_data, counts, path, title, top_n=20):
    import numpy as np
    mut_order = [m for m, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]]
    if len(mut_order) < 2 or len(all_mutation_data) < 3:
        return _ext_no_data(path, title, "Insufficient data for PCA")
    X = np.array([[1.0 if m in set(d.get("mutations", [])) else 0.0 for m in mut_order]
                   for d in all_mutation_data])
    Xc = X - X.mean(axis=0, keepdims=True)
    try:
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        coords = U[:, :2] * S[:2]
    except Exception:
        return _ext_no_data(path, title, "PCA failed")
    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.scatter(coords[:, 0], coords[:, 1], s=42, c="#1565c0", alpha=0.75, edgecolors="white")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_dendro(all_mutation_data, counts, path, title, top_n=20):
    import numpy as np
    mut_order = [m for m, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]]
    if len(mut_order) < 2 or len(all_mutation_data) < 3:
        return _ext_no_data(path, title, "Insufficient data for clustering")
    X = np.array([[1.0 if m in set(d.get("mutations", [])) else 0.0 for m in mut_order]
                   for d in all_mutation_data])
    try:
        from scipy.cluster.hierarchy import dendrogram, linkage
        Z = linkage(X, method="average", metric="hamming")
        fig, ax = plt.subplots(figsize=(max(8, 0.22 * len(all_mutation_data)), 5))
        dendrogram(Z, labels=[short_isolate_label(d["isolate"], 14) for d in all_mutation_data],
                   leaf_rotation=90, leaf_font_size=5.5, ax=ax, color_threshold=0)
        ax.set_ylabel("Hamming distance")
        ax.set_title(title)
        _safe_savefig(fig, path)
    except Exception:
        return _ext_no_data(path, title, "SciPy unavailable — dendrogram skipped")

def _chart_status_pie(all_data, path, title, mode="gene_status"):
    if not all_data:
        return _ext_no_data(path, title)
    if mode == "gene_status":
        found = sum(1 for d in all_data if d.get("gene_found"))
        disrupted = sum(1 for d in all_data if d.get("gene_found") and d.get("is_disrupted"))
        missing = sum(1 for d in all_data if not d.get("gene_found"))
        intact = max(0, found - disrupted)
        labels = ["Intact", "Disrupted", "Missing"]
        values = [intact, disrupted, missing]
        colors = ["#2e7d32", "#ef6c00", "#b0bec5"]
    else:
        with_mut = sum(1 for d in all_data if d.get("gene_found") and d.get("mutations"))
        without = sum(1 for d in all_data if d.get("gene_found") and not d.get("mutations"))
        labels = ["≥1 mutation", "No mutation"]
        values = [with_mut, without]
        colors = ["#1565c0", "#cfd8dc"]
    pairs = [(v, l, c) for v, l, c in zip(values, labels, colors) if v > 0]
    if not pairs:
        return _ext_no_data(path, title)
    values2, labels2, colors2 = zip(*pairs)
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.pie(values2, labels=labels2, autopct=lambda p: f"{p:.0f}%",
           colors=list(colors2), startangle=90, wedgeprops=dict(edgecolor="white"))
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_qc_scatter(all_data, path, title):
    pts = [(d.get("pident"), d.get("qcov"), d.get("is_disrupted", False))
           for d in all_data if d.get("gene_found") and isinstance(d.get("pident"), (int, float))]
    if not pts:
        return _ext_no_data(path, title, "No QC metrics")
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for pident, qcov, dis in pts:
        ax.scatter(pident, qcov, s=38, color="#c62828" if dis else "#1565c0",
                   alpha=0.75, edgecolor="white")
    ax.set_xlabel("Percent identity (%)")
    ax.set_ylabel("Query coverage (%)")
    ax.set_title(title)
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#1565c0", markersize=8, label="Intact"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#c62828", markersize=8, label="Disrupted"),
    ]
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=8)
    _safe_savefig(fig, path)

def _chart_metric_bar(all_data, metric, path, title):
    rows = [(d["isolate"], d.get(metric)) for d in all_data
            if d.get("gene_found") and isinstance(d.get(metric), (int, float))]
    if not rows:
        return _ext_no_data(path, title)
    rows.sort(key=lambda r: -r[1])
    labels = [short_isolate_label(r[0], 16) for r in rows]
    values = [r[1] for r in rows]
    fig, ax = plt.subplots(figsize=(max(8, 0.16 * len(rows)), 4.5))
    ax.bar(range(len(rows)), values, color="#1565c0" if metric == "pident" else "#2e7d32", width=0.85)
    if len(rows) <= 40:
        ax.set_xticks(range(len(rows)))
        ax.set_xticklabels(labels, rotation=90, fontsize=5.5)
    else:
        ax.set_xticks([])
    ax.set_ylabel("Percent identity (%)" if metric == "pident" else "Query coverage (%)")
    ax.set_xlabel(f"{len(rows)} isolates")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_threshold_bar(counts, total, thr, path, title):
    if not counts or not total:
        return _ext_no_data(path, title)
    min_c = int(total * thr / 100.0)
    common = {m: c for m, c in counts.items() if c >= min_c}
    if not common:
        return _ext_no_data(path, title, f"No mutation ≥{thr}%")
    items = sorted(common.items(), key=lambda kv: -kv[1])
    labels = [k for k, _ in items]
    vals = [100.0 * v / total for _, v in items]
    fig, ax = plt.subplots(figsize=(8, max(3, 0.28 * len(labels) + 1)))
    ax.barh(range(len(labels)), vals, color="#2e7d32", edgecolor="white")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Prevalence (%)")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_threshold_trend(counts, total, thresholds, path, title):
    if not counts or not total:
        return _ext_no_data(path, title)
    xs, ys = [], []
    for thr in thresholds:
        min_c = int(total * thr / 100.0)
        xs.append(thr)
        ys.append(sum(1 for c in counts.values() if c >= min_c))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(xs, ys, "o-", color="#6a1b9a", markersize=7, linewidth=2)
    for x, y in zip(xs, ys):
        ax.annotate(str(y), (x, y), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8)
    ax.set_xlabel("Prevalence threshold (%)")
    ax.set_ylabel("Mutations meeting threshold")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_aa_matrix(counts, path, title):
    import numpy as np
    aas = list("ACDEFGHIKLMNPQRSTVWY")
    idx = {a: i for i, a in enumerate(aas)}
    mat = [[0] * len(aas) for _ in aas]
    found = False
    for notation, c in counts.items():
        parsed = _classify_substitution(notation)
        if not parsed:
            continue
        ref_aa, _, alt_aa, _ = parsed
        if ref_aa in idx and alt_aa in idx:
            mat[idx[ref_aa]][idx[alt_aa]] += c
            found = True
    if not found:
        return _ext_no_data(path, title, "No missense substitutions")
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(mat, cmap="OrRd")
    ax.set_xticks(range(len(aas)))
    ax.set_xticklabels(aas, fontsize=8)
    ax.set_yticks(range(len(aas)))
    ax.set_yticklabels(aas, fontsize=8)
    ax.set_xlabel("Variant residue")
    ax.set_ylabel("Reference residue")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.grid(False)
    _safe_savefig(fig, path)

def _chart_physchem(counts, path, title):
    tally = {}
    any_found = False
    for notation, c in counts.items():
        parsed = _classify_substitution(notation)
        if not parsed:
            continue
        any_found = True
        change = parsed[3]
        tally[change] = tally.get(change, 0) + c
    if not any_found:
        return _ext_no_data(path, title, "No missense substitutions")
    colors = {"Conservative (same class)": "#2e7d32", "Radical (class change)": "#c62828",
              "Unclassified": "#90a4ae"}
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ks = list(tally.keys())
    ax.bar(ks, [tally[k] for k in ks], color=[colors.get(k, "#607d8b") for k in ks], edgecolor="white")
    ax.set_ylabel("Isolate-mutation events")
    ax.set_title(title)
    plt.setp(ax.get_xticklabels(), rotation=12, ha="right")
    _safe_savefig(fig, path)

def _chart_summary_table(gene_name, counts, total, all_data, path, title, label):
    found = sum(1 for d in all_data if d.get("gene_found"))
    disrupted = sum(1 for d in all_data if d.get("gene_found") and d.get("is_disrupted"))
    missing = sum(1 for d in all_data if not d.get("gene_found"))
    n_muts = len(counts)
    n_events = sum(counts.values())
    top3 = sorted(counts.items(), key=lambda kv: -kv[1])[:3]
    top3_str = ", ".join(f"{m} ({c})" for m, c in top3) if top3 else "None"
    rows = [
        ["Gene", gene_name],
        ["Category", label],
        ["Isolates analysed", str(len(all_data))],
        ["Gene detected", str(found)],
        ["Gene disrupted", str(disrupted)],
        ["Gene missing", str(missing)],
        ["Distinct mutations", str(n_muts)],
        ["Total events", str(n_events)],
        ["Top 3 mutations", top3_str],
    ]
    fig, ax = plt.subplots(figsize=(8, 0.48 * len(rows) + 1.2))
    ax.axis("off")
    tbl = ax.table(cellText=rows, colLabels=["Metric", "Value"], cellLoc="left", loc="center",
                    colWidths=[0.4, 0.6])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.55)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor(PUB_PALETTE["primary"])
            cell.set_text_props(color="white", fontweight="bold")
        else:
            cell.set_facecolor("#f5f5f5" if r % 2 == 0 else "white")
    ax.set_title(title, pad=18)
    _safe_savefig(fig, path)

def _chart_summary_panel(gene_name, counts, total, all_data, all_mutation_data, path, title):
    if not counts:
        return _ext_no_data(path, title)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    ax = axes[0, 0]
    top = sorted(counts.items(), key=lambda kv: -kv[1])[:10]
    ax.barh([m for m, _ in top][::-1], [c for _, c in top][::-1], color=PUB_PALETTE["primary"])
    ax.set_title("Top 10", fontsize=10)
    ax = axes[0, 1]
    bins = {"Rare": 0, "Low": 0, "Mod": 0, "High": 0, "Maj": 0}
    for _, c in counts.items():
        pct = 100.0 * c / total if total else 0
        key = ("Rare" if pct < 10 else "Low" if pct < 25 else "Mod" if pct < 50
               else "High" if pct < 75 else "Maj")
        bins[key] += 1
    ax.bar(list(bins.keys()), list(bins.values()), color="#4fc3f7")
    ax.set_title("Prevalence categories", fontsize=10)
    ax = axes[0, 2]
    burdens = [len(d.get("mutations", [])) for d in all_mutation_data] or [0]
    ax.hist(burdens, bins=range(0, max(burdens) + 2), color="#8e24aa", align="left")
    ax.set_title("Burden / isolate", fontsize=10)
    ax = axes[1, 0]
    found = sum(1 for d in all_data if d.get("gene_found"))
    missing = sum(1 for d in all_data if not d.get("gene_found"))
    disrupted = sum(1 for d in all_data if d.get("gene_found") and d.get("is_disrupted"))
    vals = [max(0, found - disrupted), disrupted, missing]
    labs = ["Intact", "Disrupted", "Missing"]
    nz = [(v, l) for v, l in zip(vals, labs) if v > 0]
    if nz:
        ax.pie([v for v, _ in nz], labels=[l for _, l in nz], autopct="%1.0f%%",
               colors=["#2e7d32", "#ef6c00", "#b0bec5"])
    ax.set_title("Gene status", fontsize=10)
    ax = axes[1, 1]
    items = sorted(counts.items(), key=lambda kv: -kv[1])
    cum, running = [], 0
    for _, c in items:
        running += c
        cum.append(100.0 * running / total if total else 0)
    ax.plot(range(1, len(cum) + 1), cum, color=PUB_PALETTE["warn"])
    ax.set_title("Cumulative coverage", fontsize=10)
    ax.set_ylim(0, 105)
    ax = axes[1, 2]
    ax.axis("off")
    ax.text(0.05, 0.85, f"Gene: {gene_name}", fontweight="bold", fontsize=11, transform=ax.transAxes)
    ax.text(0.05, 0.65, f"Isolates: {len(all_data)}", fontsize=9, transform=ax.transAxes)
    ax.text(0.05, 0.50, f"Distinct muts: {len(counts)}", fontsize=9, transform=ax.transAxes)
    ax.text(0.05, 0.35, f"Events: {sum(counts.values())}", fontsize=9, transform=ax.transAxes)
    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _safe_savefig(fig, path)

def _chart_area(counts, total, path, title, top_n=20, mode="percent"):
    items = _ext_items(counts, top_n)
    if not items:
        return _ext_no_data(path, title)
    labels = [k for k, _ in items]
    values = [100.0 * v / total if (mode == "percent" and total) else v for _, v in items]
    fig, ax = plt.subplots(figsize=(max(7, 0.35 * len(labels)), 4.5))
    x = range(len(labels))
    ax.fill_between(x, values, color="#b2dfdb", alpha=0.85)
    ax.plot(x, values, color="#00695c", linewidth=1.5)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Isolates (%)" if mode == "percent" else "Count")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_step(counts, total, path, title, top_n=20, mode="percent"):
    items = _ext_items(counts, top_n)
    if not items:
        return _ext_no_data(path, title)
    labels = [k for k, _ in items]
    values = [100.0 * v / total if (mode == "percent" and total) else v for _, v in items]
    fig, ax = plt.subplots(figsize=(max(7, 0.35 * len(labels)), 4.5))
    ax.step(range(len(labels)), values, where="mid", color="#283593", linewidth=1.8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Isolates (%)" if mode == "percent" else "Count")
    ax.set_title(title)
    _safe_savefig(fig, path)

def _chart_pie_top(counts, path, title, top_n=8, donut=False, explode_first=False):
    if not counts:
        return _ext_no_data(path, title)
    items = sorted(counts.items(), key=lambda kv: -kv[1])
    top = items[:top_n]
    other = sum(v for _, v in items[top_n:])
    labels = [k for k, _ in top] + (["Other"] if other > 0 else [])
    values = [v for _, v in top] + ([other] if other > 0 else [])
    if not values or sum(values) == 0:
        return _ext_no_data(path, title)
    cmap = plt.get_cmap("tab20")
    colors = [cmap(i % 20) for i in range(len(values))]
    explode = [0.06] + [0] * (len(values) - 1) if explode_first and len(values) > 1 else None
    fig, ax = plt.subplots(figsize=(7, 7))
    wp = dict(width=0.42, edgecolor="white") if donut else dict(edgecolor="white")
    wedges, _, _ = ax.pie(values, labels=None, autopct=lambda p: f"{p:.1f}%" if p >= 3 else "",
                          colors=colors, startangle=90, pctdistance=0.78,
                          explode=explode, wedgeprops=wp)
    ax.legend(wedges, labels, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7.5, frameon=False)
    ax.set_title(title)
    _safe_savefig(fig, path)

def _ext_chart_manifest(gene_name, counts, total, all_data, all_mutation_data,
                         type_map, position_map, thresholds, label):
    thr_list = thresholds or [25, 30, 35, 40, 50]
    m = []

    m.append(("001_Mutation_Frequency_Bar_Plot.png",
              lambda p: _chart_mut_freq_bar(counts, total, p, f"{gene_name} — Mutation Frequency Bar [{label}]")))
    m.append(("002_Horizontal_Mutation_Bar_Plot.png",
              lambda p: _chart_mut_freq_hbar(counts, total, p, f"{gene_name} — Horizontal Mutation Bar [{label}]")))
    m.append(("003_Sorted_Mutation_Frequency_Plot.png",
              lambda p: _chart_sorted_freq(counts, total, p, f"{gene_name} — Sorted Mutation Frequency [{label}]")))
    m.append(("004_Percentage_Mutation_Bar_Plot.png",
              lambda p: _chart_mut_pct_bar(counts, total, p, f"{gene_name} — Percentage Mutation Bar [{label}]")))
    m.append(("005_Percentage_Mutation_Horizontal.png",
              lambda p: _chart_mut_pct_hbar(counts, total, p, f"{gene_name} — Percentage Horizontal [{label}]")))
    m.append(("006_Mutation_Count_Dot_Plot.png",
              lambda p: _chart_dot_freq(counts, total, p, f"{gene_name} — Mutation Count Dot Plot [{label}]")))
    m.append(("007_Mutation_Frequency_Lollipop.png",
              lambda p: _chart_lollipop(counts, total, p, f"{gene_name} — Mutation Frequency Lollipop [{label}]", mode="count")))
    m.append(("008_Mutation_Prevalence_Lollipop.png",
              lambda p: _chart_lollipop(counts, total, p, f"{gene_name} — Mutation Prevalence Lollipop [{label}]", mode="percent")))
    m.append(("009_Mutation_Stem_Plot.png",
              lambda p: _chart_stem(counts, total, p, f"{gene_name} — Mutation Stem Plot [{label}]", mode="count")))
    m.append(("010_Mutation_Stem_Percent.png",
              lambda p: _chart_stem(counts, total, p, f"{gene_name} — Mutation Stem % [{label}]", mode="percent")))
    m.append(("011_Mutation_Needle_Plot.png",
              lambda p: _chart_needle(counts, total, p, f"{gene_name} — Mutation Needle Plot [{label}]")))
    m.append(("012_Mutation_Frequency_Histogram.png",
              lambda p: _chart_hist_counts(counts, p, f"{gene_name} — Mutation Frequency Histogram [{label}]")))
    m.append(("013_Mutation_Density_Plot.png",
              lambda p: _chart_density_counts(counts, p, f"{gene_name} — Mutation Density Plot [{label}]")))
    m.append(("014_Rug_Plot_of_Mutations.png",
              lambda p: _chart_rug(counts, p, f"{gene_name} — Rug Plot of Mutations [{label}]")))
    m.append(("015_Mutation_Strip_Plot.png",
              lambda p: _chart_strip(counts, p, f"{gene_name} — Mutation Strip Plot [{label}]")))
    m.append(("016_Bubble_Mutation_Plot.png",
              lambda p: _chart_bubble(counts, total, p, f"{gene_name} — Bubble Mutation Plot [{label}]")))
    m.append(("017_Sized_Point_Mutation_Plot.png",
              lambda p: _chart_sized_point(counts, total, p, f"{gene_name} — Sized-Point Mutation Plot [{label}]")))
    m.append(("018_Area_Chart_Top20.png",
              lambda p: _chart_area(counts, total, p, f"{gene_name} — Area Chart Top 20 [{label}]")))
    m.append(("019_Step_Chart_Top20.png",
              lambda p: _chart_step(counts, total, p, f"{gene_name} — Step Chart Top 20 [{label}]")))

    for n in (10, 15, 20, 25, 30):
        m.append((f"02{n}_Top{n}_Frequency_Bar.png",
                  lambda p, n=n: _chart_mut_freq_bar(counts, total, p, f"{gene_name} — Top {n} Frequency [{label}]", top_n=n)))
        m.append((f"03{n}_Top{n}_Percent_HBar.png",
                  lambda p, n=n: _chart_mut_pct_hbar(counts, total, p, f"{gene_name} — Top {n} % Horizontal [{label}]", top_n=n)))
        m.append((f"04{n}_Top{n}_Lollipop.png",
                  lambda p, n=n: _chart_lollipop(counts, total, p, f"{gene_name} — Top {n} Lollipop [{label}]", top_n=n)))
        m.append((f"05{n}_Top{n}_Stem.png",
                  lambda p, n=n: _chart_stem(counts, total, p, f"{gene_name} — Top {n} Stem [{label}]", top_n=n)))
        m.append((f"06{n}_Top{n}_Bubble.png",
                  lambda p, n=n: _chart_bubble(counts, total, p, f"{gene_name} — Top {n} Bubble [{label}]", top_n=n)))

    m.append(("070_Binary_Mutation_Matrix.png",
              lambda p: _chart_binary_matrix(all_mutation_data, counts, p, f"{gene_name} — Binary Mutation Matrix [{label}]")))
    m.append(("071_Gene_x_Isolate_Heatmap.png",
              lambda p: _chart_presence_heatmap(all_mutation_data, counts, p, f"{gene_name} — Gene × Isolate Heatmap [{label}]", cmap="Blues")))
    m.append(("072_Mutation_Burden_Heatmap.png",
              lambda p: _chart_burden_heatmap(all_mutation_data, counts, p, f"{gene_name} — Mutation Burden Heatmap [{label}]")))
    m.append(("073_Mutation_Intensity_Heatmap.png",
              lambda p: _chart_intensity_heatmap(all_mutation_data, counts, p, f"{gene_name} — Mutation Intensity Heatmap [{label}]")))
    m.append(("074_Dot_Matrix_Top15.png",
              lambda p: _chart_dotmatrix(all_mutation_data, counts, p, f"{gene_name} — Dot Matrix Top 15 [{label}]", top_n=15)))
    m.append(("075_Dot_Matrix_Top20.png",
              lambda p: _chart_dotmatrix(all_mutation_data, counts, p, f"{gene_name} — Dot Matrix Top 20 [{label}]", top_n=20)))
    m.append(("076_Dot_Matrix_Top25.png",
              lambda p: _chart_dotmatrix(all_mutation_data, counts, p, f"{gene_name} — Dot Matrix Top 25 [{label}]", top_n=25)))
    for n, cmap in [(15, "Greens"), (20, "Oranges"), (25, "Purples"), (18, "YlOrRd")]:
        m.append((f"08{n}_Presence_Heatmap_{cmap}.png",
                  lambda p, n=n, c=cmap: _chart_presence_heatmap(
                      all_mutation_data, counts, p,
                      f"{gene_name} — Presence Heatmap top {n} ({c}) [{label}]", top_n=n, cmap=c)))

    m.append(("090_Mutation_Type_Pie.png",
              lambda p: _chart_type_pie(counts, type_map, p, f"{gene_name} — Mutation Type Pie [{label}]")))
    m.append(("091_Mutation_Type_Bar.png",
              lambda p: _chart_type_bar(counts, type_map, p, f"{gene_name} — Mutation Type Bar [{label}]")))
    m.append(("092_Mutation_Type_Donut.png",
              lambda p: _chart_type_donut(counts, type_map, p, f"{gene_name} — Mutation Type Donut [{label}]")))
    m.append(("093_Mutation_Type_Stacked_Share.png",
              lambda p: _chart_type_stacked_dummy(counts, type_map, p, f"{gene_name} — Mutation Type Stacked Share [{label}]")))
    m.append(("094_Pie_Top8.png",
              lambda p: _chart_pie_top(counts, p, f"{gene_name} — Pie Top 8 [{label}]", top_n=8)))
    m.append(("095_Donut_Top8.png",
              lambda p: _chart_pie_top(counts, p, f"{gene_name} — Donut Top 8 [{label}]", top_n=8, donut=True)))
    m.append(("096_Pie_Exploded_Top8.png",
              lambda p: _chart_pie_top(counts, p, f"{gene_name} — Pie Exploded Top 8 [{label}]", top_n=8, explode_first=True)))
    m.append(("097_Pie_Top12.png",
              lambda p: _chart_pie_top(counts, p, f"{gene_name} — Pie Top 12 [{label}]", top_n=12)))
    m.append(("098_Donut_Top12.png",
              lambda p: _chart_pie_top(counts, p, f"{gene_name} — Donut Top 12 [{label}]", top_n=12, donut=True)))

    m.append(("100_SNP_Position_Scatter.png",
              lambda p: _chart_position_scatter(counts, position_map, p, f"{gene_name} — SNP Position Scatter [{label}]")))
    m.append(("101_SNP_Position_Lollipop.png",
              lambda p: _chart_position_lollipop_local(counts, position_map, p, f"{gene_name} — SNP Position Lollipop [{label}]")))
    m.append(("102_SNP_Position_Histogram.png",
              lambda p: _chart_position_hist_local(counts, position_map, p, f"{gene_name} — SNP Position Histogram [{label}]")))
    m.append(("103_SNP_Frequency_by_Position.png",
              lambda p: _chart_snp_freq_by_position(counts, position_map, p, f"{gene_name} — SNP Frequency-by-Position [{label}]")))

    m.append(("110_Mutation_Burden_Histogram.png",
              lambda p: _chart_burden_hist(all_mutation_data, p, f"{gene_name} — Mutation Burden Histogram [{label}]")))
    m.append(("111_Mutation_Burden_Boxplot.png",
              lambda p: _chart_burden_box(all_mutation_data, p, f"{gene_name} — Mutation Burden Boxplot [{label}]")))
    m.append(("112_Mutation_Burden_Violin.png",
              lambda p: _chart_burden_violin(all_mutation_data, p, f"{gene_name} — Mutation Burden Violin [{label}]")))
    m.append(("113_Mutation_Burden_ECDF.png",
              lambda p: _chart_burden_ecdf(all_mutation_data, p, f"{gene_name} — Mutation Burden ECDF [{label}]")))
    m.append(("114_Mutation_Burden_Rank_Plot.png",
              lambda p: _chart_burden_rank(all_mutation_data, p, f"{gene_name} — Mutation Burden Rank Plot [{label}]")))
    m.append(("115_Mutation_Burden_Waterfall.png",
              lambda p: _chart_burden_waterfall(all_mutation_data, p, f"{gene_name} — Mutation Burden Waterfall [{label}]")))

    m.append(("120_Wilson_CI_Top15.png",
              lambda p: _chart_ci_bar(counts, total, p, f"{gene_name} — Wilson CI Top 15 [{label}]", top_n=15)))
    m.append(("121_Wilson_CI_Top20.png",
              lambda p: _chart_ci_bar(counts, total, p, f"{gene_name} — Wilson CI Top 20 [{label}]", top_n=20)))
    m.append(("122_Radar_Top8.png",
              lambda p: _chart_radar(counts, total, p, f"{gene_name} — Radar Top 8 [{label}]", top_n=8)))
    m.append(("123_Radar_Top10.png",
              lambda p: _chart_radar(counts, total, p, f"{gene_name} — Radar Top 10 [{label}]", top_n=10)))
    m.append(("124_Waterfall_Top12.png",
              lambda p: _chart_waterfall(counts, total, p, f"{gene_name} — Waterfall Top 12 [{label}]", top_n=12)))
    m.append(("125_Pareto_Top15.png",
              lambda p: _chart_pareto(counts, total, p, f"{gene_name} — Pareto Top 15 [{label}]", top_n=15)))
    m.append(("126_Cumulative_Coverage.png",
              lambda p: _chart_cumulative(counts, total, p, f"{gene_name} — Cumulative Coverage [{label}]")))
    m.append(("127_Rank_Frequency_LogLog.png",
              lambda p: _chart_rank_freq(counts, p, f"{gene_name} — Rank–Frequency (log-log) [{label}]")))
    m.append(("128_Prevalence_Category_Bar.png",
              lambda p: _chart_category_bar(counts, total, p, f"{gene_name} — Prevalence Category Bar [{label}]")))

    m.append(("130_CoOccurrence_Top10.png",
              lambda p: _chart_cooc(all_mutation_data, counts, p, f"{gene_name} — Co-Occurrence Top 10 [{label}]", top_n=10)))
    m.append(("131_CoOccurrence_Top12.png",
              lambda p: _chart_cooc(all_mutation_data, counts, p, f"{gene_name} — Co-Occurrence Top 12 [{label}]", top_n=12, cmap="Blues")))
    m.append(("132_CoOccurrence_Top15.png",
              lambda p: _chart_cooc(all_mutation_data, counts, p, f"{gene_name} — Co-Occurrence Top 15 [{label}]", top_n=15, cmap="Greens")))
    m.append(("133_Isolate_PCA.png",
              lambda p: _chart_pca(all_mutation_data, counts, p, f"{gene_name} — Isolate PCA [{label}]")))
    m.append(("134_Isolate_PCA_Top30.png",
              lambda p: _chart_pca(all_mutation_data, counts, p, f"{gene_name} — Isolate PCA Top 30 [{label}]", top_n=30)))
    m.append(("135_Isolate_Dendrogram.png",
              lambda p: _chart_dendro(all_mutation_data, counts, p, f"{gene_name} — Isolate Dendrogram [{label}]")))

    m.append(("140_Gene_Status_Pie.png",
              lambda p: _chart_status_pie(all_data, p, f"{gene_name} — Gene Status Pie [{label}]", mode="gene_status")))
    m.append(("141_Mutation_Presence_Pie.png",
              lambda p: _chart_status_pie(all_data, p, f"{gene_name} — Mutation Presence Pie [{label}]", mode="mutation_presence")))
    m.append(("142_QC_pident_vs_qcov.png",
              lambda p: _chart_qc_scatter(all_data, p, f"{gene_name} — QC pident vs qcov [{label}]")))
    m.append(("143_pident_by_Isolate.png",
              lambda p: _chart_metric_bar(all_data, "pident", p, f"{gene_name} — pident by Isolate [{label}]")))
    m.append(("144_qcov_by_Isolate.png",
              lambda p: _chart_metric_bar(all_data, "qcov", p, f"{gene_name} — qcov by Isolate [{label}]")))

    for thr in thr_list[:5]:
        m.append((f"15{thr}_Threshold_ge{thr}pct.png",
                  lambda p, t=thr: _chart_threshold_bar(
                      counts, total, t, p, f"{gene_name} — Mutations ≥{t}% [{label}]")))
    m.append(("160_Threshold_Trend.png",
              lambda p: _chart_threshold_trend(counts, total, thr_list, p, f"{gene_name} — Threshold Trend [{label}]")))

    m.append(("170_AA_Substitution_Matrix.png",
              lambda p: _chart_aa_matrix(counts, p, f"{gene_name} — AA Substitution Matrix [{label}]")))
    m.append(("171_Conservative_vs_Radical.png",
              lambda p: _chart_physchem(counts, p, f"{gene_name} — Conservative vs Radical [{label}]")))

    m.append(("180_Summary_Panel.png",
              lambda p: _chart_summary_panel(gene_name, counts, total, all_data, all_mutation_data, p,
                  f"{gene_name} — Summary Panel [{label}]")))
    m.append(("181_Summary_Table.png",
              lambda p: _chart_summary_table(gene_name, counts, total, all_data, p,
                  f"{gene_name} — Summary Table [{label}]", label)))

    for n in (8, 12, 16, 18, 22, 24):
        m.append((f"19{n}_Top{n}_Needle.png",
                  lambda p, n=n: _chart_needle(counts, total, p, f"{gene_name} — Top {n} Needle [{label}]", top_n=n)))
        m.append((f"20{n}_Top{n}_Area.png",
                  lambda p, n=n: _chart_area(counts, total, p, f"{gene_name} — Top {n} Area [{label}]", top_n=n)))
        m.append((f"21{n}_Top{n}_Step.png",
                  lambda p, n=n: _chart_step(counts, total, p, f"{gene_name} — Top {n} Step [{label}]", top_n=n)))
        m.append((f"22{n}_Top{n}_Strip.png",
                  lambda p, n=n: _chart_strip(counts, p, f"{gene_name} — Top {n} Strip [{label}]", top_n=n)))
        m.append((f"23{n}_Top{n}_SizedPoint.png",
                  lambda p, n=n: _chart_sized_point(counts, total, p, f"{gene_name} — Top {n} Sized-Point [{label}]", top_n=n)))
        m.append((f"24{n}_Top{n}_Waterfall.png",
                  lambda p, n=n: _chart_waterfall(counts, total, p, f"{gene_name} — Top {n} Waterfall [{label}]", top_n=n)))
        m.append((f"25{n}_Top{n}_Pareto.png",
                  lambda p, n=n: _chart_pareto(counts, total, p, f"{gene_name} — Top {n} Pareto [{label}]", top_n=n)))

    return m

def _ext_run_battery(gene_name, mutation_counts, total_isolates, all_data, all_mutation_data,
                      curated_set, out_dir, thresholds, relevant_only: bool):
    all_counts, relevant_counts = _split_relevant(gene_name, mutation_counts, curated_set)
    counts = relevant_counts if relevant_only else all_counts
    label = "NG AMR-relevant only" if relevant_only else "all detected mutations"
    if relevant_only:
        filtered_amd = []
        for d in all_mutation_data or []:
            keep = [x for x in d.get("mutations", []) if x in counts]
            filtered_amd.append({"isolate": d["isolate"], "mutations": keep})
    else:
        filtered_amd = all_mutation_data or []
    position_map, type_map = _notation_metadata(all_data)
    os.makedirs(out_dir, exist_ok=True)
    manifest = _ext_chart_manifest(gene_name, counts, total_isolates, all_data, filtered_amd,
                                    type_map, position_map, thresholds, label)
    ok, failed = 0, 0
    for fname, fn in manifest:
        out_path = os.path.join(out_dir, fname)
        try:
            fn(out_path)
            ok += 1
        except Exception as e:
            failed += 1
            logger.warning(f"    chart failed [{gene_name}] {fname}: {e}")
            try:
                _ext_no_data(out_path, f"{gene_name} — {fname}", f"Error: {e}")
            except Exception as _exc:
                logger.debug("Ignored non-critical exception: %s", _exc)
                pass
    return ok, failed, len(manifest)

def generate_extended_gene_wise_images(gene_results: List[dict], output_root: str,
                                        curated_master: Optional[dict] = None,
                                        thresholds: Optional[List[int]] = None) -> str:
    if not MATPLOTLIB_OK:
        logger.info("  [Extended_Gene_Figures] matplotlib missing — skipped")
        return ""
    _pub_rc()
    root = os.path.join(output_root, "Extended_Gene_Figures")
    os.makedirs(root, exist_ok=True)
    logger.info("\n" + "=" * 70)
    logger.info("EXTENDED GENE-WISE FIGURES (≥120 distinct styles × All / AMR-relevant)")
    logger.info("=" * 70)
    for res in gene_results:
        if not res:
            continue
        gene_name = res.get("gene")
        counts = res.get("mutation_counts") or {}
        total = res.get("total_isolates") or 0
        all_data = res.get("all_data") or []
        all_mutation_data = res.get("all_mutation_data") or []
        curated_set = res.get("curated_set") or set()
        if curated_master:
            curated_set = set(curated_set) | set(curated_master.get(gene_name, []))
        folder_name = GENE_FOLDER_NAMES.get(gene_name, gene_name)
        gene_root = os.path.join(root, folder_name)
        all_dir = os.path.join(gene_root, "All_Mutations")
        rel_dir = os.path.join(gene_root, "NG_AMR_Relevant")
        ok1, fail1, n1 = _ext_run_battery(gene_name, counts, total, all_data, all_mutation_data,
                                           curated_set, all_dir, thresholds, relevant_only=False)
        ok2, fail2, n2 = _ext_run_battery(gene_name, counts, total, all_data, all_mutation_data,
                                           curated_set, rel_dir, thresholds, relevant_only=True)
        logger.info(f"  {gene_name:16s}  All_Mutations: {ok1}/{n1}   NG_AMR_Relevant: {ok2}/{n2}"
              + (f"  (failures: {fail1 + fail2})" if (fail1 + fail2) else ""))
    logger.info(f"\nExtended gene-wise galleries → {root}")
    return root

def generate_extended_master_publication_figures(gene_results: List[dict], output_root: str,
                                                  curated_master: Optional[dict] = None) -> str:
    if not MATPLOTLIB_OK:
        logger.info("  [Extended_Publication_Figures] matplotlib missing — skipped")
        return ""
    _pub_rc()
    root = os.path.join(output_root, "Extended_Publication_Figures")
    all_dir = os.path.join(root, "All_Mutations")
    rel_dir = os.path.join(root, "NG_AMR_Relevant")
    os.makedirs(all_dir, exist_ok=True)
    os.makedirs(rel_dir, exist_ok=True)

    valid = [r for r in gene_results if r and r.get("mutation_counts") is not None]

    def _gene_rows(relevant_only):
        rows = []
        for r in valid:
            gene = r["gene"]
            cset = set(r.get("curated_set") or set())
            if curated_master:
                cset |= set(curated_master.get(gene, []))
            all_c, rel_c = _split_relevant(gene, r["mutation_counts"], cset)
            c = rel_c if relevant_only else all_c
            rows.append((gene, c, r.get("total_isolates") or 0, r))
        return rows

    def _build(relevant_only, out_dir, label):
        rows = _gene_rows(relevant_only)
        gene_names = [g for g, _, _, _ in rows]
        n_distinct = [len(c) for _, c, _, _ in rows]
        n_events = [sum(c.values()) for _, c, _, _ in rows]
        mean_burden = [(sum(c.values()) / t) if t else 0 for _, c, t, _ in rows]

        fig, ax = plt.subplots(figsize=(max(9, 0.38 * len(gene_names)), 5))
        ax.bar(gene_names, n_distinct, color="#1565c0", edgecolor="white")
        ax.set_ylabel("Distinct mutations")
        ax.set_title(f"Distinct mutations per gene — {label}")
        plt.setp(ax.get_xticklabels(), rotation=55, ha="right", fontsize=7.5)
        _safe_savefig(fig, os.path.join(out_dir, "001_distinct_mutations_bar.png"))

        fig, ax = plt.subplots(figsize=(8, max(4, 0.28 * len(gene_names) + 1)))
        ax.barh(gene_names, n_distinct, color="#2e7d32", edgecolor="white")
        ax.set_xlabel("Distinct mutations")
        ax.set_title(f"Distinct mutations (horizontal) — {label}")
        _safe_savefig(fig, os.path.join(out_dir, "002_distinct_mutations_hbar.png"))

        fig, ax = plt.subplots(figsize=(max(9, 0.38 * len(gene_names)), 5))
        ax.bar(gene_names, n_events, color="#e65100", edgecolor="white")
        ax.set_ylabel("Isolate-mutation events")
        ax.set_title(f"Mutation events per gene — {label}")
        plt.setp(ax.get_xticklabels(), rotation=55, ha="right", fontsize=7.5)
        _safe_savefig(fig, os.path.join(out_dir, "003_events_bar.png"))

        fig, ax = plt.subplots(figsize=(8, max(4, 0.28 * len(gene_names) + 1)))
        ax.barh(gene_names, n_events, color="#6a1b9a", edgecolor="white")
        ax.set_xlabel("Events")
        ax.set_title(f"Events (horizontal) — {label}")
        _safe_savefig(fig, os.path.join(out_dir, "004_events_hbar.png"))

        fig, ax = plt.subplots(figsize=(max(9, 0.38 * len(gene_names)), 5))
        ax.bar(gene_names, mean_burden, color="#00897b", edgecolor="white")
        ax.set_ylabel("Mean mutations / isolate")
        ax.set_title(f"Mean burden per gene — {label}")
        plt.setp(ax.get_xticklabels(), rotation=55, ha="right", fontsize=7.5)
        _safe_savefig(fig, os.path.join(out_dir, "005_mean_burden_bar.png"))

        items = sorted(zip(gene_names, n_distinct), key=lambda x: x[1])
        labs = [x[0] for x in items]
        vals = [x[1] for x in items]
        fig, ax = plt.subplots(figsize=(8, max(3.5, 0.28 * len(labs) + 1)))
        ax.hlines(range(len(labs)), 0, vals, color="#455a64", linewidth=1.2)
        ax.plot(vals, range(len(labs)), "o", color="#c62828", markersize=6)
        ax.set_yticks(range(len(labs)))
        ax.set_yticklabels(labs, fontsize=7.5)
        ax.set_xlabel("Distinct mutations")
        ax.set_title(f"Lollipop — distinct mutations — {label}")
        _safe_savefig(fig, os.path.join(out_dir, "006_distinct_lollipop.png"))

        items = sorted(zip(gene_names, n_events), key=lambda x: -x[1])
        labs = [x[0] for x in items]
        vals = [x[1] for x in items]
        fig, ax = plt.subplots(figsize=(max(8, 0.35 * len(labs)), 4.8))
        ax.stem(range(len(labs)), vals, linefmt="#6a1b9a", markerfmt="o", basefmt=" ")
        ax.set_xticks(range(len(labs)))
        ax.set_xticklabels(labs, rotation=55, ha="right", fontsize=7)
        ax.set_ylabel("Events")
        ax.set_title(f"Stem — events by gene — {label}")
        _safe_savefig(fig, os.path.join(out_dir, "007_events_stem.png"))

        items = sorted(zip(gene_names, n_distinct), key=lambda x: -x[1])
        labs = [x[0] for x in items]
        vals = [x[1] for x in items]
        fig, ax = plt.subplots(figsize=(max(8, 0.35 * len(labs)), 4.5))
        ax.fill_between(range(len(labs)), vals, color="#b2dfdb", alpha=0.85)
        ax.plot(range(len(labs)), vals, color="#00695c", linewidth=1.4)
        ax.set_xticks(range(len(labs)))
        ax.set_xticklabels(labs, rotation=55, ha="right", fontsize=7)
        ax.set_ylabel("Distinct mutations")
        ax.set_title(f"Area — distinct mutations — {label}")
        _safe_savefig(fig, os.path.join(out_dir, "008_distinct_area.png"))

        fig, ax = plt.subplots(figsize=(max(8, 0.35 * len(labs)), 4.5))
        ax.step(range(len(labs)), vals, where="mid", color="#283593", linewidth=1.6)
        ax.set_xticks(range(len(labs)))
        ax.set_xticklabels(labs, rotation=55, ha="right", fontsize=7)
        ax.set_ylabel("Distinct mutations")
        ax.set_title(f"Step — distinct mutations — {label}")
        _safe_savefig(fig, os.path.join(out_dir, "009_distinct_step.png"))

        pooled, gene_of = {}, {}
        for gene, c, _, _ in rows:
            for m, cnt in c.items():
                key = f"{gene}:{m}"
                pooled[key] = cnt
                gene_of[key] = gene
        for top_n in (15, 20, 25, 30):
            top = sorted(pooled.items(), key=lambda kv: -kv[1])[:top_n]
            if not top:
                _ext_no_data(os.path.join(out_dir, f"01{top_n}_top{top_n}_pooled.png"),
                             f"Top {top_n} pooled — {label}")
                continue
            cmap = plt.get_cmap("tab20")
            genes_u = sorted(set(gene_of[k] for k, _ in top))
            gcolor = {g: cmap(i % 20) for i, g in enumerate(genes_u)}
            fig, ax = plt.subplots(figsize=(9, max(3.5, 0.28 * len(top) + 1)))
            labels = [k for k, _ in top]
            vals = [v for _, v in top]
            colors = [gcolor[gene_of[k]] for k, _ in top]
            ax.barh(range(len(labels)), vals, color=colors)
            ax.set_yticks(range(len(labels)))
            ax.set_yticklabels(labels, fontsize=6.5)
            ax.invert_yaxis()
            ax.set_xlabel("Isolates")
            ax.set_title(f"Top {top_n} mutations across genes — {label}")
            from matplotlib.patches import Patch
            handles = [Patch(color=gcolor[g], label=g) for g in genes_u]
            ax.legend(handles=handles, loc="lower right", fontsize=6.5, ncol=2, frameon=False)
            _safe_savefig(fig, os.path.join(out_dir, f"01{top_n}_top{top_n}_pooled.png"))

        pct_with = []
        for gene, c, t, r in rows:
            if t == 0:
                pct_with.append(0)
                continue
            amd = r.get("all_mutation_data") or []
            keep = set(c.keys())
            n_with = sum(1 for d in amd if set(d.get("mutations", [])) & keep)
            pct_with.append(100.0 * n_with / t)
        fig, ax = plt.subplots(figsize=(max(9, 0.38 * len(gene_names)), 5))
        ax.bar(gene_names, pct_with, color="#ef6c00", edgecolor="white")
        ax.set_ylabel("% isolates with ≥1 mutation")
        ax.set_ylim(0, 105)
        ax.set_title(f"Isolate-level prevalence — {label}")
        plt.setp(ax.get_xticklabels(), rotation=55, ha="right", fontsize=7.5)
        _safe_savefig(fig, os.path.join(out_dir, "020_pct_with_mutation.png"))

        found, disrupted, missing = [], [], []
        for r in valid:
            ad = r.get("all_data") or []
            f = sum(1 for d in ad if d.get("gene_found"))
            dch = sum(1 for d in ad if d.get("gene_found") and d.get("is_disrupted"))
            mi = sum(1 for d in ad if not d.get("gene_found"))
            found.append(max(0, f - dch))
            disrupted.append(dch)
            missing.append(mi)
        gnames = [r["gene"] for r in valid]
        fig, ax = plt.subplots(figsize=(max(9, 0.38 * len(gnames)), 5))
        ax.bar(gnames, found, label="Intact", color="#2e7d32")
        ax.bar(gnames, disrupted, bottom=found, label="Disrupted", color="#ef6c00")
        bottom2 = [a + b for a, b in zip(found, disrupted)]
        ax.bar(gnames, missing, bottom=bottom2, label="Missing", color="#b0bec5")
        ax.set_ylabel("Isolates")
        ax.set_title("Gene status (stacked) — all genes")
        ax.legend(frameon=False, fontsize=8)
        plt.setp(ax.get_xticklabels(), rotation=55, ha="right", fontsize=7.5)
        _safe_savefig(fig, os.path.join(out_dir, "021_gene_status_stacked.png"))

        import numpy as np
        fig, ax = plt.subplots(figsize=(max(9, 0.38 * len(gnames)), 5))
        x = np.arange(len(gnames))
        w = 0.25
        ax.bar(x - w, found, w, label="Intact", color="#2e7d32")
        ax.bar(x, disrupted, w, label="Disrupted", color="#ef6c00")
        ax.bar(x + w, missing, w, label="Missing", color="#b0bec5")
        ax.set_xticks(x)
        ax.set_xticklabels(gnames)
        ax.set_ylabel("Isolates")
        ax.set_title("Gene status (grouped) — all genes")
        ax.legend(frameon=False, fontsize=8)
        plt.setp(ax.get_xticklabels(), rotation=55, ha="right", fontsize=7.5)
        _safe_savefig(fig, os.path.join(out_dir, "022_gene_status_grouped.png"))

        gene_share = {g: sum(c.values()) for g, c, _, _ in rows}
        gene_share = {g: v for g, v in gene_share.items() if v > 0}
        for donut, fname in [(False, "023_gene_share_pie.png"), (True, "024_gene_share_donut.png")]:
            if not gene_share:
                _ext_no_data(os.path.join(out_dir, fname), f"Gene share — {label}")
                continue
            fig, ax = plt.subplots(figsize=(7.5, 7.5))
            cmap = plt.get_cmap("tab20")
            colors = [cmap(i % 20) for i in range(len(gene_share))]
            wp = dict(width=0.42, edgecolor="white") if donut else dict(edgecolor="white")
            ax.pie(gene_share.values(), labels=gene_share.keys(), autopct="%1.0f%%",
                   colors=colors, startangle=90, textprops={"fontsize": 7}, wedgeprops=wp)
            ax.set_title(f"Share of mutation events by gene — {label}")
            _safe_savefig(fig, os.path.join(out_dir, fname))

        for max_per, cmap_name in [(4, "YlOrRd"), (6, "YlGnBu"), (8, "PuRd"), (5, "OrRd")]:
            cell_data = []
            for gene, c, t, _ in rows:
                topk = sorted(c.items(), key=lambda kv: -kv[1])[:max_per]
                cell_data.append((gene, t, topk))
            max_cols = max((len(x[2]) for x in cell_data), default=0)
            if not max_cols:
                _ext_no_data(os.path.join(out_dir, f"03{max_per}_top{max_per}_heatmap.png"),
                             f"Top{max_per} heatmap — {label}")
                continue
            mat = np.zeros((len(rows), max_cols))
            annot = [["" for _ in range(max_cols)] for _ in range(len(rows))]
            for i, (gene, t, topk) in enumerate(cell_data):
                for j, (m, cnt) in enumerate(topk):
                    mat[i, j] = 100.0 * cnt / t if t else 0
                    annot[i][j] = m
            fig, ax = plt.subplots(figsize=(max(7, 0.85 * max_cols), max(5.5, 0.32 * len(rows) + 1)))
            im = ax.imshow(mat, cmap=cmap_name, aspect="auto")
            ax.set_yticks(range(len(gene_names)))
            ax.set_yticklabels(gene_names, fontsize=7.5)
            ax.set_xticks(range(max_cols))
            ax.set_xticklabels([f"R{i+1}" for i in range(max_cols)], fontsize=8)
            for i in range(len(rows)):
                for j in range(max_cols):
                    if annot[i][j]:
                        ax.text(j, i, annot[i][j], ha="center", va="center", fontsize=5.5)
            fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Prevalence (%)")
            ax.set_title(f"Top-{max_per} mutation per gene — {label}")
            ax.grid(False)
            _safe_savefig(fig, os.path.join(out_dir, f"03{max_per}_top{max_per}_heatmap_{cmap_name}.png"))

        rows_tbl = [["Gene", "Isolates", "Distinct", "Events", "% w/ mut"]]
        for (gene, c, t, _), pw in zip(rows, pct_with):
            rows_tbl.append([gene, str(t), str(len(c)), str(sum(c.values())), f"{pw:.1f}%"])
        fig, ax = plt.subplots(figsize=(9, 0.38 * len(rows_tbl) + 1))
        ax.axis("off")
        tbl = ax.table(cellText=rows_tbl[1:], colLabels=rows_tbl[0], cellLoc="center", loc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8)
        tbl.scale(1, 1.35)
        for (r, c_), cell in tbl.get_celld().items():
            if r == 0:
                cell.set_facecolor(PUB_PALETTE["primary"])
                cell.set_text_props(color="white", fontweight="bold")
        ax.set_title(f"Cross-gene summary — {label}", pad=14)
        _safe_savefig(fig, os.path.join(out_dir, "040_cross_gene_summary_table.png"))

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes[0, 0].bar(gene_names, n_distinct, color=PUB_PALETTE["primary"])
        axes[0, 0].set_title("Distinct mutations")
        plt.setp(axes[0, 0].get_xticklabels(), rotation=55, ha="right", fontsize=6)
        axes[0, 1].bar(gene_names, n_events, color=PUB_PALETTE["accent"])
        axes[0, 1].set_title("Events")
        plt.setp(axes[0, 1].get_xticklabels(), rotation=55, ha="right", fontsize=6)
        axes[1, 0].bar(gene_names, mean_burden, color=PUB_PALETTE["curated"])
        axes[1, 0].set_title("Mean burden")
        plt.setp(axes[1, 0].get_xticklabels(), rotation=55, ha="right", fontsize=6)
        axes[1, 1].bar(gene_names, pct_with, color=PUB_PALETTE["warn"])
        axes[1, 1].set_title("% isolates with mut")
        axes[1, 1].set_ylim(0, 105)
        plt.setp(axes[1, 1].get_xticklabels(), rotation=55, ha="right", fontsize=6)
        fig.suptitle(f"Cross-gene overview panel — {label}", fontsize=13, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        _safe_savefig(fig, os.path.join(out_dir, "041_overview_panel.png"))

        for vals, ylab, tname in [
            (n_distinct, "Distinct muts", "distinct"),
            (n_events, "Events", "events"),
            (mean_burden, "Mean burden", "burden"),
            (pct_with, "% with mut", "pct"),
        ]:
            fig, ax = plt.subplots(figsize=(max(9, 0.38 * len(gene_names)), 4.5))
            ax.plot(range(len(gene_names)), vals, "o-", color="#1565c0", markersize=5)
            ax.fill_between(range(len(gene_names)), vals, alpha=0.15, color="#1565c0")
            ax.set_xticks(range(len(gene_names)))
            ax.set_xticklabels(gene_names, rotation=55, ha="right", fontsize=7)
            ax.set_ylabel(ylab)
            ax.set_title(f"Line view — {tname} — {label}")
            _safe_savefig(fig, os.path.join(out_dir, f"05{tname}_line.png"))

        for vals, name in [(n_distinct, "distinct"), (n_events, "events")]:
            if max(vals or [0]) <= 0:
                continue
            fig, ax = plt.subplots(figsize=(max(9, 0.38 * len(gene_names)), 4.8))
            ax.bar(gene_names, [max(v, 0.1) for v in vals], color="#e65100", edgecolor="white")
            ax.set_yscale("log")
            ax.set_ylabel(name + " (log)")
            ax.set_title(f"Log-scale {name} — {label}")
            plt.setp(ax.get_xticklabels(), rotation=55, ha="right", fontsize=7.5)
            _safe_savefig(fig, os.path.join(out_dir, f"06{name}_log.png"))

        if len(gene_names) >= 3:
            for metric_vals, mname in [(n_distinct, "distinct"), (pct_with, "pct_with")]:
                n = min(12, len(gene_names))
                labs = gene_names[:n]
                vals = list(metric_vals[:n])
                mx = max(vals) or 1
                vals_n = [100.0 * v / mx for v in vals]
                vals_n += vals_n[:1]
                angles = [k / float(n) * 2 * math.pi for k in range(n)]
                angles += angles[:1]
                fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))
                ax.plot(angles, vals_n, color="#6a1b9a", linewidth=1.5)
                ax.fill(angles, vals_n, color="#e1bee7", alpha=0.25)
                ax.set_xticks(angles[:-1])
                ax.set_xticklabels(labs, fontsize=7)
                ax.set_title(f"Gene radar — {mname} — {label}", pad=18)
                _safe_savefig(fig, os.path.join(out_dir, f"07{mname}_radar.png"))

        extra_i = 0
        for top_n in (8, 10, 12, 15, 18, 20, 25, 30):
            top = sorted(pooled.items(), key=lambda kv: -kv[1])[:top_n]
            if not top:
                continue
            fig, ax = plt.subplots(figsize=(8, max(3.2, 0.26 * len(top) + 1)))
            ax.barh([k for k, _ in top][::-1], [v for _, v in top][::-1],
                    color="#1565c0", edgecolor="white")
            ax.set_xlabel("Isolates")
            ax.set_title(f"Top {top_n} pooled — {label}")
            _safe_savefig(fig, os.path.join(out_dir, f"08{extra_i:02d}_top{top_n}_pooled_hbar.png"))
            extra_i += 1

        fig, ax = plt.subplots(figsize=(max(8, 0.35 * len(gene_names)), 4.8))
        ax.vlines(range(len(gene_names)), 0, n_distinct, color="#1565c0", linewidth=1.5)
        ax.scatter(range(len(gene_names)), n_distinct, s=28, color="#c62828", zorder=3)
        ax.set_xticks(range(len(gene_names)))
        ax.set_xticklabels(gene_names, rotation=55, ha="right", fontsize=7)
        ax.set_ylabel("Distinct mutations")
        ax.set_title(f"Needle — distinct by gene — {label}")
        _safe_savefig(fig, os.path.join(out_dir, "090_distinct_needle.png"))

        sizes = [max(40, 8 * v) for v in n_events]
        fig, ax = plt.subplots(figsize=(max(8, 0.35 * len(gene_names)), 4.8))
        ax.scatter(range(len(gene_names)), n_events, s=sizes, c=n_events, cmap="YlOrRd",
                   alpha=0.8, edgecolors="white")
        ax.set_xticks(range(len(gene_names)))
        ax.set_xticklabels(gene_names, rotation=55, ha="right", fontsize=7)
        ax.set_ylabel("Events")
        ax.set_title(f"Bubble — events by gene — {label}")
        _safe_savefig(fig, os.path.join(out_dir, "091_events_bubble.png"))

        for i, (vals, name) in enumerate([(n_distinct, "distinct"), (n_events, "events"),
                                           (mean_burden, "burden"), (pct_with, "pct")]):
            items = sorted(zip(gene_names, vals), key=lambda x: x[1])
            labs = [x[0] for x in items]
            v = [x[1] for x in items]
            fig, ax = plt.subplots(figsize=(8, max(3.5, 0.28 * len(labs) + 1)))
            ax.plot(v, range(len(labs)), "o", color="#1565c0", markersize=7)
            ax.set_yticks(range(len(labs)))
            ax.set_yticklabels(labs, fontsize=7.5)
            ax.set_xlabel(name)
            ax.set_title(f"Cleveland — {name} — {label}")
            _safe_savefig(fig, os.path.join(out_dir, f"10{i}_{name}_cleveland.png"))

            fig, ax = plt.subplots(figsize=(8, max(3.5, 0.28 * len(labs) + 1)))
            ax.barh(range(len(labs)), v, color="#00897b", edgecolor="white")
            ax.set_yticks(range(len(labs)))
            ax.set_yticklabels(labs, fontsize=7.5)
            ax.set_xlabel(name)
            ax.set_title(f"Ascending bar — {name} — {label}")
            _safe_savefig(fig, os.path.join(out_dir, f"11{i}_{name}_asc_bar.png"))

    _build(False, all_dir, "all detected mutations, whole dataset")
    _build(True, rel_dir, "NG AMR-relevant only, whole dataset")
    logger.info(f"\nExtended master publication figures → {root}")
    return root

def build_master_figures_pptx(output_root: str,
                               gene_results: Optional[List[dict]] = None) -> Optional[str]:
    try:
        from pptx import Presentation
        from pptx.util import Inches, Pt
        from pptx.dml.color import RGBColor
        from pptx.enum.text import PP_ALIGN
    except ImportError as e:
        logger.error("  [PPTX] python-pptx not installed (or import failed) — "
              "skip NG_AMR_Master_Mutation_Analysis_Figures.pptx")
        logger.info(f"         Detail: {e}")
        logger.info("         Install with:  pip install python-pptx")
        return None

    search_roots = [
        os.path.join(output_root, "Extended_Gene_Figures"),
        os.path.join(output_root, "Extended_Publication_Figures"),
        os.path.join(output_root, "genes_wise_images"),
        os.path.join(output_root, "Publication_Figures"),
        os.path.join(output_root, "master_analysis_figures"),
        os.path.join(output_root, "gene_phenotype_figures"),
        os.path.join(output_root, "Exclusive_Figures"),
        os.path.join(output_root, "NG_MIC_EST_DD_vs_AD", "Figures"),
        os.path.join(output_root, "NG_AMR_MIC_Comparison_Illustration"),
        os.path.join(output_root, "Final_for_Thesis"),
        os.path.join(output_root, "Thesis_Publication_Package"),
    ]

    images: List[Tuple[str, str]] = []
    for root in search_roots:
        if not os.path.isdir(root):
            continue
        for dirpath, _, files in os.walk(root):
            for f in sorted(files):
                if f.lower().endswith((".png", ".jpg", ".jpeg", ".tif", ".tiff")):
                    abs_p = os.path.join(dirpath, f)
                    rel = os.path.relpath(abs_p, output_root)
                    images.append((rel, abs_p))

    if not images:
        logger.info("  [PPTX] No figure images found under output root — PPTX skipped")
        return None

    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    blank = prs.slide_layouts[6]

    slide = prs.slides.add_slide(blank)
    box = slide.shapes.add_textbox(Inches(0.8), Inches(2.2), Inches(11.7), Inches(2.5))
    tf = box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = "NG AMR Master Mutation Analysis — Figure Gallery"
    p.font.size = Pt(28)
    p.font.bold = True
    p.font.color.rgb = RGBColor(0x15, 0x65, 0xC0)
    p.alignment = PP_ALIGN.CENTER
    p2 = tf.add_paragraph()
    p2.text = (f"{len(images)} publication-quality figures\n"
               "Images can be moved / resized / replaced inside PowerPoint\n"
               "(raster PNGs — individual text/icons are not natively editable)")
    p2.font.size = Pt(14)
    p2.font.color.rgb = RGBColor(0x42, 0x42, 0x42)
    p2.alignment = PP_ALIGN.CENTER

    for idx, (rel, abs_p) in enumerate(images, 1):
        slide = prs.slides.add_slide(blank)
        cap = slide.shapes.add_textbox(Inches(0.3), Inches(0.15), Inches(12.7), Inches(0.45))
        ctf = cap.text_frame
        ctf.word_wrap = True
        cp = ctf.paragraphs[0]
        cp.text = f"[{idx}/{len(images)}]  {rel}"
        cp.font.size = Pt(10)
        cp.font.color.rgb = RGBColor(0x37, 0x47, 0x4F)

        try:
            slide.shapes.add_picture(abs_p, Inches(0.3), Inches(0.7),
                                     width=Inches(12.7))
        except Exception as e:
            note = slide.shapes.add_textbox(Inches(1), Inches(3), Inches(11), Inches(1))
            note.text_frame.paragraphs[0].text = f"Could not embed: {rel}\n{e}"
            note.text_frame.paragraphs[0].font.size = Pt(12)

    out_path = os.path.join(output_root, "NG_AMR_Master_Mutation_Analysis_Figures.pptx")
    try:
        prs.save(out_path)
        logger.info(f"\n  [PPTX] Figure gallery written → {out_path}")
        logger.info(f"         ({len(images)} images embedded; open in PowerPoint to move/resize)")
        return out_path
    except Exception as e:
        logger.error(f"  [PPTX] Failed to save: {e}")
        return None

def build_master_figures_pptx_native_editable(
    output_root: str,
    gene_results: Optional[List[dict]] = None,
    curated_master: Optional[dict] = None,
) -> Optional[str]:
    if not gene_results:
        logger.info("  [PPTX-Editable] No gene_results — skip native editable PPTX")
        return None

    try:
        from pptx import Presentation
        from pptx.util import Inches, Pt
        from pptx.dml.color import RGBColor
        from pptx.enum.text import PP_ALIGN
        from pptx.enum.shapes import MSO_SHAPE
        from pptx.chart.data import CategoryChartData
        from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION
    except ImportError as e:
        logger.info("  [PPTX-Editable] python-pptx (or dependency) missing — skip")
        logger.info(f"         Detail: {e}")
        logger.info("         Install with:  pip install python-pptx")
        return None

    def _ensure_run(paragraph):
        if paragraph.runs:
            return paragraph.runs[0]
        return paragraph.add_run()

    def _set_run_font(run, size=12, bold=False, color=None, name="Calibri"):
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.name = name
        if color is not None:
            run.font.color.rgb = color

    def _style_paragraph(paragraph, text, size=12, bold=False, color=None, align=None):
        paragraph.text = text if text is not None else ""
        if align is not None:
            paragraph.alignment = align
        run = _ensure_run(paragraph)
        if text and not run.text:
            run.text = str(text)
        _set_run_font(run, size=size, bold=bold, color=color)
        return run

    def _add_title(slide, text, top=0.2, height=0.55, size=20, color=RGBColor(0x15, 0x65, 0xC0)):
        box = slide.shapes.add_textbox(Inches(0.35), Inches(top), Inches(12.6), Inches(height))
        tf = box.text_frame
        tf.word_wrap = True
        _style_paragraph(tf.paragraphs[0], text, size=size, bold=True, color=color, align=PP_ALIGN.LEFT)
        return box

    def _add_subtitle(slide, text, top=0.7, size=11, color=RGBColor(0x42, 0x42, 0x42)):
        box = slide.shapes.add_textbox(Inches(0.35), Inches(top), Inches(12.6), Inches(0.4))
        tf = box.text_frame
        tf.word_wrap = True
        _style_paragraph(tf.paragraphs[0], text, size=size, bold=False, color=color)
        return box

    def _fill_table_cell(cell, text, bold=False, size=9, fill_rgb=None, font_rgb=None, align=PP_ALIGN.CENTER):
        cell.text = str(text) if text is not None else ""
        for p in cell.text_frame.paragraphs:
            p.alignment = align
            if p.runs:
                for r in p.runs:
                    _set_run_font(r, size=size, bold=bold, color=font_rgb or RGBColor(0x21, 0x21, 0x21))
            else:
                r = p.add_run()
                r.text = str(text) if text is not None else ""
                _set_run_font(r, size=size, bold=bold, color=font_rgb or RGBColor(0x21, 0x21, 0x21))
        if fill_rgb is not None:
            cell.fill.solid()
            cell.fill.fore_color.rgb = fill_rgb

    def _style_chart(chart, has_legend=True):
        try:
            chart.has_legend = has_legend
            if has_legend and chart.legend is not None:
                chart.legend.position = XL_LEGEND_POSITION.BOTTOM
                chart.legend.include_in_layout = False
        except Exception as _exc:
            logger.debug("Ignored non-critical exception: %s", _exc)
            pass

    def _add_col(slide, title, cats, series_name, values, ctype=XL_CHART_TYPE.COLUMN_CLUSTERED, legend=False):
        if not cats or not values:
            return
        _add_title(slide, title, size=16)
        _add_subtitle(slide, "Native chart — double-click to edit data / colours")
        cd = CategoryChartData()
        cd.categories = list(cats)
        cd.add_series(series_name, list(values))
        ch = slide.shapes.add_chart(ctype, Inches(0.4), Inches(1.15), Inches(12.5), Inches(5.8), cd).chart
        _style_chart(ch, has_legend=legend)

    def _add_multi_col(slide, title, cats, series_dict, ctype=XL_CHART_TYPE.COLUMN_CLUSTERED):
        if not cats or not series_dict:
            return
        _add_title(slide, title, size=16)
        _add_subtitle(slide, "Native multi-series chart — fully editable")
        cd = CategoryChartData()
        cd.categories = list(cats)
        for name, vals in series_dict.items():
            cd.add_series(name, list(vals))
        ch = slide.shapes.add_chart(ctype, Inches(0.4), Inches(1.15), Inches(12.5), Inches(5.8), cd).chart
        _style_chart(ch, has_legend=True)

    def _add_pie(slide, title, cats, values):
        if not cats or not values or sum(values) <= 0:
            return
        _add_title(slide, title, size=16)
        _add_subtitle(slide, "Native pie — every slice is editable")
        cd = CategoryChartData()
        cd.categories = list(cats)
        cd.add_series("Share", list(values))
        ch = slide.shapes.add_chart(XL_CHART_TYPE.PIE, Inches(2.2), Inches(1.15), Inches(8.8), Inches(5.8), cd).chart
        _style_chart(ch, has_legend=True)

    BLUE = RGBColor(0x15, 0x65, 0xC0)
    DARK = RGBColor(0x21, 0x21, 0x21)
    GREY = RGBColor(0x61, 0x61, 0x61)
    GREEN = RGBColor(0x2E, 0x7D, 0x32)
    HEADER_FILL = RGBColor(0x15, 0x65, 0xC0)
    ALT_ROW = RGBColor(0xE3, 0xF2, 0xFD)
    WHITE = RGBColor(0xFF, 0xFF, 0xFF)

    rows = []
    for gr in gene_results:
        if not gr:
            continue
        gene = gr.get("gene") or "?"
        counts = dict(gr.get("mutation_counts") or {})
        if curated_master:
            cset = set(gr.get("curated_set") or set()) | set(curated_master.get(gene, []))
        else:
            cset = set(gr.get("curated_set") or set())
        total = int(gr.get("total_isolates") or 0)
        all_data = gr.get("all_data") or []
        amd = gr.get("all_mutation_data") or []
        distinct = len(counts)
        events = sum(counts.values())
        if total > 0 and amd:
            keep = set(counts.keys())
            n_with = sum(1 for d in amd if set(d.get("mutations") or []) & keep)
            pct_with = 100.0 * n_with / total
        else:
            n_with = 0
            pct_with = 0.0
        found = sum(1 for d in all_data if d.get("gene_found"))
        disrupted = sum(1 for d in all_data if d.get("gene_found") and d.get("is_disrupted"))
        missing = sum(1 for d in all_data if not d.get("gene_found"))
        intact = max(0, found - disrupted)
        mean_burden = (events / total) if total else 0.0
        top_muts = sorted(counts.items(), key=lambda kv: -kv[1])
        rel_counts = {m: c for m, c in counts.items() if m in cset} if cset else {}
        non_counts = {m: c for m, c in counts.items() if m not in cset}
        type_tally: Dict[str, int] = {}
        for d in all_data:
            for m in d.get("mutations") or []:
                notation = m.get("mutation") or m.get("notation")
                if not notation or notation not in counts:
                    continue
                t = str(m.get("type") or "substitution").replace("_", " ").title()
                type_tally[t] = type_tally.get(t, 0) + 1
        if not type_tally and counts:
            type_tally = {"Substitution / other": events}
        burdens = [len(d.get("mutations") or []) for d in amd] if amd else []
        burden_hist: Dict[str, int] = {}
        for b in burdens:
            k = str(b)
            burden_hist[k] = burden_hist.get(k, 0) + 1
        cat_bins = {"Rare (<10%)": 0, "Low (10–25%)": 0, "Moderate (25–50%)": 0,
                    "High (50–75%)": 0, "Majority (>75%)": 0}
        if total:
            for _, c in counts.items():
                pct = 100.0 * c / total
                if pct < 10:
                    cat_bins["Rare (<10%)"] += 1
                elif pct < 25:
                    cat_bins["Low (10–25%)"] += 1
                elif pct < 50:
                    cat_bins["Moderate (25–50%)"] += 1
                elif pct < 75:
                    cat_bins["High (50–75%)"] += 1
                else:
                    cat_bins["Majority (>75%)"] += 1
        rows.append({
            "gene": gene, "distinct": distinct, "events": events, "total": total,
            "pct_with": pct_with, "n_with": n_with, "intact": intact,
            "disrupted": disrupted, "missing": missing, "mean_burden": mean_burden,
            "counts": counts, "top_muts": top_muts, "curated_set": cset,
            "rel_counts": rel_counts, "non_counts": non_counts,
            "type_tally": type_tally, "burden_hist": burden_hist, "cat_bins": cat_bins,
            "all_data": all_data, "amd": amd,
        })

    if not rows:
        logger.info("  [PPTX-Editable] No usable gene rows — skip")
        return None

    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    blank = prs.slide_layouts[6]
    gene_names = [r["gene"] for r in rows]

    slide = prs.slides.add_slide(blank)
    shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0), Inches(0), Inches(13.333), Inches(1.4))
    shape.fill.solid()
    shape.fill.fore_color.rgb = BLUE
    shape.line.fill.background()
    box = slide.shapes.add_textbox(Inches(0.6), Inches(0.3), Inches(12), Inches(0.95))
    tf = box.text_frame
    _style_paragraph(tf.paragraphs[0], "NG AMR Master Mutation Analysis",
                     size=26, bold=True, color=WHITE)
    p2 = tf.add_paragraph()
    _style_paragraph(p2, "OPTION A — Maximum native editable coverage  •  no PNG embeds",
                     size=13, bold=False, color=RGBColor(0xBB, 0xDE, 0xFB))
    info = slide.shapes.add_textbox(Inches(0.7), Inches(1.9), Inches(11.8), Inches(4.8))
    tf = info.text_frame
    tf.word_wrap = True
    for i, line in enumerate([
        f"Genes: {len(rows)}   |   Max isolates: {max((r['total'] for r in rows), default=0)}",
        " ",
        "This deck contains ONLY native PowerPoint objects derived from pipeline data:",
        "  • Tables — every cell editable",
        "  • Charts — series, categories, colours, axis titles fully editable",
        "  • Multiple encodings of the same metrics (column / bar / line / pie / stacked)",
        "  • Cross-gene overview + rich per-gene sections",
        " ",
        "Existing PNG generators are untouched. Companion PNG gallery PPTX is separate.",
        "Open in PowerPoint → click any label, bar, slice, or cell to edit.",
    ]):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        _style_paragraph(para, line, size=13, bold=(i == 0), color=DARK)

    slide = prs.slides.add_slide(blank)
    _add_title(slide, "Cross-gene summary (editable table)")
    _add_subtitle(slide, "All values are live text cells")
    headers = ["Gene", "Isolates", "Distinct", "Events", "% with ≥1", "Mean burden", "Intact", "Disrupted", "Missing"]
    n_rows = len(rows) + 1
    table_shape = slide.shapes.add_table(n_rows, len(headers), Inches(0.3), Inches(1.15), Inches(12.7),
                                         Inches(min(5.9, 0.30 * n_rows + 0.4)))
    table = table_shape.table
    for c, h in enumerate(headers):
        _fill_table_cell(table.cell(0, c), h, bold=True, size=9, fill_rgb=HEADER_FILL, font_rgb=WHITE)
    for i, r in enumerate(rows):
        vals = [r["gene"], r["total"], r["distinct"], r["events"], f"{r['pct_with']:.1f}%",
                f"{r['mean_burden']:.2f}", r["intact"], r["disrupted"], r["missing"]]
        fill = ALT_ROW if i % 2 == 0 else WHITE
        for c, v in enumerate(vals):
            _fill_table_cell(table.cell(i + 1, c), v, bold=(c == 0), size=8, fill_rgb=fill)

    metrics = [
        ("distinct", "Distinct mutations", [r["distinct"] for r in rows]),
        ("events", "Isolate–mutation events", [r["events"] for r in rows]),
        ("pct", "% isolates with ≥1 mutation", [round(r["pct_with"], 1) for r in rows]),
        ("burden", "Mean mutations / isolate", [round(r["mean_burden"], 2) for r in rows]),
    ]
    for key, label, vals in metrics:
        slide = prs.slides.add_slide(blank)
        _add_col(slide, f"{label} — column", gene_names, label, vals)
        items = sorted(zip(gene_names, vals), key=lambda x: x[1])
        slide = prs.slides.add_slide(blank)
        _add_col(slide, f"{label} — horizontal (sorted)", [x[0] for x in items], label, [x[1] for x in items],
                 ctype=XL_CHART_TYPE.BAR_CLUSTERED)
        slide = prs.slides.add_slide(blank)
        _add_col(slide, f"{label} — line", gene_names, label, vals, ctype=XL_CHART_TYPE.LINE)

    slide = prs.slides.add_slide(blank)
    _add_multi_col(slide, "Gene status — grouped (Intact / Disrupted / Missing)",
                   gene_names,
                   {"Intact": [r["intact"] for r in rows],
                    "Disrupted": [r["disrupted"] for r in rows],
                    "Missing": [r["missing"] for r in rows]})
    slide = prs.slides.add_slide(blank)
    _add_multi_col(slide, "Gene status — stacked",
                   gene_names,
                   {"Intact": [r["intact"] for r in rows],
                    "Disrupted": [r["disrupted"] for r in rows],
                    "Missing": [r["missing"] for r in rows]},
                   ctype=XL_CHART_TYPE.COLUMN_STACKED)

    pie_rows = [r for r in rows if r["events"] > 0]
    if pie_rows:
        slide = prs.slides.add_slide(blank)
        _add_pie(slide, "Share of mutation events by gene",
                 [r["gene"] for r in pie_rows], [r["events"] for r in pie_rows])

    pooled = []
    for r in rows:
        for m, c in r["top_muts"][:15]:
            pooled.append((f"{r['gene']}:{m}", c, r["gene"]))
    pooled.sort(key=lambda x: -x[1])
    for top_n in (15, 25, 40):
        top = pooled[:top_n]
        if not top:
            continue
        slide = prs.slides.add_slide(blank)
        _add_title(slide, f"Top {top_n} mutations across genes (editable table)", size=16)
        _add_subtitle(slide, "Ranked by isolate count")
        headers2 = ["Rank", "Gene:Mutation", "Count"]
        n_r = len(top) + 1
        tbl = slide.shapes.add_table(n_r, 3, Inches(1.5), Inches(1.2), Inches(10), Inches(min(5.7, 0.26 * n_r + 0.3))).table
        for c, h in enumerate(headers2):
            _fill_table_cell(tbl.cell(0, c), h, bold=True, size=10, fill_rgb=HEADER_FILL, font_rgb=WHITE)
        for i, (lab, cnt, _) in enumerate(top):
            fill = ALT_ROW if i % 2 == 0 else WHITE
            for c, v in enumerate([i + 1, lab, cnt]):
                _fill_table_cell(tbl.cell(i + 1, c), v, bold=(c == 1), size=9, fill_rgb=fill)
        slide = prs.slides.add_slide(blank)
        _add_col(slide, f"Top {top_n} mutations across genes — bar",
                 [x[0] for x in top], "Isolates", [x[1] for x in top],
                 ctype=XL_CHART_TYPE.BAR_CLUSTERED)

    logger.info(f"  [PPTX-Editable] Building per-gene native sections for {len(rows)} genes…")
    for ri, r in enumerate(rows, 1):
        gene = r["gene"]
        total = r["total"]
        counts = r["counts"]
        top_all = r["top_muts"]
        cset = r["curated_set"]
        if ri % 5 == 0 or ri == len(rows):
            logger.info(f"    … gene {ri}/{len(rows)} ({gene})")

        slide = prs.slides.add_slide(blank)
        shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0), Inches(2.5), Inches(13.333), Inches(1.6))
        shape.fill.solid()
        shape.fill.fore_color.rgb = BLUE
        shape.line.fill.background()
        box = slide.shapes.add_textbox(Inches(0.5), Inches(2.7), Inches(12.3), Inches(1.2))
        tf = box.text_frame
        _style_paragraph(tf.paragraphs[0], f"{gene}", size=28, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
        p2 = tf.add_paragraph()
        _style_paragraph(p2,
            f"n={total}  •  {r['distinct']} distinct  •  {r['events']} events  •  {r['pct_with']:.1f}% with ≥1 mut",
            size=13, bold=False, color=RGBColor(0xBB, 0xDE, 0xFB), align=PP_ALIGN.CENTER)

        slide = prs.slides.add_slide(blank)
        _add_title(slide, f"{gene} — summary metrics", size=16)
        metrics_rows = [
            ["Metric", "Value"],
            ["Gene", gene],
            ["Isolates", str(total)],
            ["Detected / Intact / Disrupted / Missing",
             f"{r['intact']+r['disrupted']} / {r['intact']} / {r['disrupted']} / {r['missing']}"],
            ["Distinct mutations", str(r["distinct"])],
            ["Events", str(r["events"])],
            ["Isolates with ≥1 mut", f"{r['n_with']} ({r['pct_with']:.1f}%)"],
            ["Mean burden", f"{r['mean_burden']:.2f}"],
            ["AMR-relevant (curated)", str(len(r["rel_counts"]))],
            ["Non-curated", str(len(r["non_counts"]))],
        ]
        tbl = slide.shapes.add_table(len(metrics_rows), 2, Inches(1.8), Inches(1.2), Inches(9.5),
                                     Inches(min(5.8, 0.38 * len(metrics_rows)))).table
        for ri2, row in enumerate(metrics_rows):
            for ci, val in enumerate(row):
                fill = HEADER_FILL if ri2 == 0 else (ALT_ROW if ri2 % 2 == 0 else WHITE)
                font_c = WHITE if ri2 == 0 else DARK
                _fill_table_cell(tbl.cell(ri2, ci), val, bold=(ri2 == 0 or ci == 0), size=11,
                                 fill_rgb=fill, font_rgb=font_c,
                                 align=PP_ALIGN.LEFT if ci == 0 else PP_ALIGN.CENTER)

        def _mut_table(title, items, note, max_rows=20):
            if not items:
                return
            slide = prs.slides.add_slide(blank)
            _add_title(slide, title, size=15)
            _add_subtitle(slide, note)
            headers = ["Rank", "Mutation", "Count", "Prevalence %", "Curated"]
            show = items[:max_rows]
            n_r = len(show) + 1
            tbl = slide.shapes.add_table(n_r, 5, Inches(0.4), Inches(1.15), Inches(12.5),
                                         Inches(min(5.8, 0.27 * n_r + 0.25))).table
            for c, h in enumerate(headers):
                _fill_table_cell(tbl.cell(0, c), h, bold=True, size=10, fill_rgb=HEADER_FILL, font_rgb=WHITE)
            for i, (mut, cnt) in enumerate(show):
                pct = (100.0 * cnt / total) if total else 0.0
                is_cur = "Yes" if mut in cset else "No"
                fill = ALT_ROW if i % 2 == 0 else WHITE
                cur_c = GREEN if is_cur == "Yes" else GREY
                for c, v in enumerate([i + 1, mut, cnt, f"{pct:.1f}%", is_cur]):
                    _fill_table_cell(tbl.cell(i + 1, c), v, bold=(c == 1), size=9, fill_rgb=fill,
                                     font_rgb=(cur_c if c == 4 else DARK))

        if top_all:
            _mut_table(f"{gene} — top mutations (all)", top_all, f"{len(top_all)} distinct — editable", 25)
            for start in range(25, len(top_all), 25):
                part = top_all[start:start + 25]
                _mut_table(f"{gene} — mutation catalogue {start+1}–{start+len(part)} of {len(top_all)}",
                           part, "Complete ranked list", 25)

        top_rel = sorted(r["rel_counts"].items(), key=lambda kv: -kv[1])
        if top_rel:
            _mut_table(f"{gene} — top AMR-relevant (curated)", top_rel, f"{len(top_rel)} curated", 20)

        for n, ctype, suffix in [
            (10, XL_CHART_TYPE.COLUMN_CLUSTERED, "top10 column"),
            (15, XL_CHART_TYPE.COLUMN_CLUSTERED, "top15 column"),
            (15, XL_CHART_TYPE.BAR_CLUSTERED, "top15 bar"),
            (20, XL_CHART_TYPE.BAR_CLUSTERED, "top20 bar"),
            (30, XL_CHART_TYPE.BAR_CLUSTERED, "top30 bar"),
        ]:
            show = top_all[:n]
            if not show:
                continue
            slide = prs.slides.add_slide(blank)
            _add_col(slide, f"{gene} — {suffix} (count)",
                     [m for m, _ in show], "Isolates", [c for _, c in show], ctype=ctype)
            if total:
                slide = prs.slides.add_slide(blank)
                _add_col(slide, f"{gene} — {suffix} (prevalence %)",
                         [m for m, _ in show], "Prevalence %",
                         [round(100.0 * c / total, 1) for _, c in show], ctype=ctype)

        if top_rel:
            show = top_rel[:15]
            slide = prs.slides.add_slide(blank)
            _add_col(slide, f"{gene} — top AMR-relevant (count)",
                     [m for m, _ in show], "Isolates", [c for _, c in show],
                     ctype=XL_CHART_TYPE.BAR_CLUSTERED)

        if top_all:
            top8 = top_all[:8]
            other = sum(c for _, c in top_all[8:])
            cats = [m for m, _ in top8] + (["Other"] if other else [])
            vals = [c for _, c in top8] + ([other] if other else [])
            slide = prs.slides.add_slide(blank)
            _add_pie(slide, f"{gene} — event share (top 8 + Other)", cats, vals)

        status_pairs = [(r["intact"], "Intact"), (r["disrupted"], "Disrupted"), (r["missing"], "Missing")]
        status_pairs = [(v, l) for v, l in status_pairs if v > 0]
        if status_pairs:
            slide = prs.slides.add_slide(blank)
            _add_pie(slide, f"{gene} — gene detection status",
                     [l for _, l in status_pairs], [v for v, _ in status_pairs])
        if total:
            slide = prs.slides.add_slide(blank)
            _add_pie(slide, f"{gene} — isolates with vs without mutation",
                     ["≥1 mutation", "No mutation"], [r["n_with"], max(0, total - r["n_with"])])

        if r["type_tally"]:
            items_t = sorted(r["type_tally"].items(), key=lambda kv: -kv[1])
            slide = prs.slides.add_slide(blank)
            _add_pie(slide, f"{gene} — mutation type distribution",
                     [k for k, _ in items_t], [v for _, v in items_t])
            slide = prs.slides.add_slide(blank)
            _add_col(slide, f"{gene} — mutation type counts",
                     [k for k, _ in items_t], "Events", [v for _, v in items_t])

        cat_items = [(k, v) for k, v in r["cat_bins"].items() if v > 0]
        if cat_items:
            slide = prs.slides.add_slide(blank)
            _add_col(slide, f"{gene} — mutations by prevalence category",
                     [k for k, _ in cat_items], "Distinct mutations", [v for _, v in cat_items])

        cur_e = sum(r["rel_counts"].values())
        non_e = sum(r["non_counts"].values())
        if cur_e or non_e:
            slide = prs.slides.add_slide(blank)
            _add_pie(slide, f"{gene} — curated vs non-curated events",
                     ["AMR-relevant", "Non-curated"], [cur_e, non_e])
            slide = prs.slides.add_slide(blank)
            _add_col(slide, f"{gene} — curated vs non-curated (distinct + events)",
                     ["Distinct curated", "Distinct non-curated", "Events curated", "Events non-curated"],
                     "Count", [len(r["rel_counts"]), len(r["non_counts"]), cur_e, non_e])

        if r["burden_hist"]:
            keys = sorted(r["burden_hist"].keys(), key=lambda x: int(x) if str(x).isdigit() else 0)
            slide = prs.slides.add_slide(blank)
            _add_col(slide, f"{gene} — mutation burden per isolate",
                     keys, "Isolates", [r["burden_hist"][k] for k in keys])

        thr_list = [25, 30, 35, 40, 50]
        if total and counts:
            thr_counts = []
            for thr in thr_list:
                min_c = int(total * thr / 100.0)
                common = sorted(((m, c) for m, c in counts.items() if c >= min_c), key=lambda kv: -kv[1])
                thr_counts.append(len(common))
                if common:
                    _mut_table(f"{gene} — mutations ≥ {thr}%", common,
                               f"min count {min_c}  •  {len(common)} mutations", 20)
                    show = common[:20]
                    slide = prs.slides.add_slide(blank)
                    _add_col(slide, f"{gene} — ≥{thr}% (top {len(show)}) prevalence %",
                             [m for m, _ in show], "Prevalence %",
                             [round(100.0 * c / total, 1) for _, c in show],
                             ctype=XL_CHART_TYPE.BAR_CLUSTERED)
            slide = prs.slides.add_slide(blank)
            _add_col(slide, f"{gene} — mutations meeting prevalence thresholds",
                     [f"≥{t}%" for t in thr_list], "Distinct mutations", thr_counts)

    slide = prs.slides.add_slide(blank)
    _add_title(slide, "How to edit this deck (OPTION A)", size=22)
    tips = slide.shapes.add_textbox(Inches(0.7), Inches(1.15), Inches(11.8), Inches(5.5))
    tf = tips.text_frame
    tf.word_wrap = True
    for i, line in enumerate([
        "Tables: click any cell → type to change values or labels.",
        "Charts: double-click → Edit Data opens the embedded spreadsheet.",
        "        Right-click a series → Format Data Series to change colours.",
        "Titles: click any title or subtitle text box and type.",
        " ",
        "This is OPTION A (maximum native coverage). All slides are native objects",
        "built from pipeline data. Existing PNG generators were not modified.",
        "Companion PNG gallery: NG_AMR_Master_Mutation_Analysis_Figures.pptx",
    ]):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        _style_paragraph(para, line, size=13, bold=False, color=DARK)

    out_path = os.path.join(output_root, "NG_AMR_Master_Mutation_Analysis_Figures_Editable.pptx")
    try:
        prs.save(out_path)
        logger.info(f"\n  [PPTX-Editable] OPTION A deck written → {out_path}")
        logger.info(f"         ({len(prs.slides)} native editable slides; no PNG embeds)")
        return out_path
    except Exception as e:
        logger.error(f"  [PPTX-Editable] Failed to save: {e}")
        return None

def build_gene_specific_editable_pptx_folder(
    output_root: str,
    gene_results: Optional[List[dict]] = None,
    curated_master: Optional[dict] = None,
    thresholds: Optional[List[int]] = None,
) -> Optional[str]:
    if not gene_results:
        logger.info("  [Gene-PPTX] No gene_results — skip gene-specific editable PPTX folder")
        return None

    try:
        from pptx import Presentation
        from pptx.util import Inches, Pt
        from pptx.dml.color import RGBColor
        from pptx.enum.text import PP_ALIGN
        from pptx.enum.shapes import MSO_SHAPE
        from pptx.chart.data import CategoryChartData
        from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION
    except ImportError as e:
        logger.info("  [Gene-PPTX] python-pptx missing — skip gene-specific editable PPTX")
        logger.info(f"         Detail: {e}")
        return None

    thr_list = thresholds or [25, 30, 35, 40, 50]
    out_dir = os.path.join(output_root, "Gene_Editable_PPTX")
    os.makedirs(out_dir, exist_ok=True)

    def _ensure_run(paragraph):
        if paragraph.runs:
            return paragraph.runs[0]
        return paragraph.add_run()

    def _set_run_font(run, size=12, bold=False, color=None, name="Calibri"):
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.name = name
        if color is not None:
            run.font.color.rgb = color

    def _style_paragraph(paragraph, text, size=12, bold=False, color=None, align=None):
        paragraph.text = text if text is not None else ""
        if align is not None:
            paragraph.alignment = align
        run = _ensure_run(paragraph)
        if text and not run.text:
            run.text = str(text)
        _set_run_font(run, size=size, bold=bold, color=color)
        return run

    def _add_title(slide, text, top=0.2, height=0.55, size=16, color=RGBColor(0x15, 0x65, 0xC0)):
        box = slide.shapes.add_textbox(Inches(0.35), Inches(top), Inches(12.6), Inches(height))
        tf = box.text_frame
        tf.word_wrap = True
        _style_paragraph(tf.paragraphs[0], text, size=size, bold=True, color=color, align=PP_ALIGN.LEFT)
        return box

    def _add_subtitle(slide, text, top=0.7, size=11, color=RGBColor(0x42, 0x42, 0x42)):
        box = slide.shapes.add_textbox(Inches(0.35), Inches(top), Inches(12.6), Inches(0.4))
        tf = box.text_frame
        tf.word_wrap = True
        _style_paragraph(tf.paragraphs[0], text, size=size, bold=False, color=color)
        return box

    def _fill_table_cell(cell, text, bold=False, size=9, fill_rgb=None, font_rgb=None, align=PP_ALIGN.CENTER):
        cell.text = str(text) if text is not None else ""
        for p in cell.text_frame.paragraphs:
            p.alignment = align
            if p.runs:
                for r in p.runs:
                    _set_run_font(r, size=size, bold=bold, color=font_rgb or RGBColor(0x21, 0x21, 0x21))
            else:
                r = p.add_run()
                r.text = str(text) if text is not None else ""
                _set_run_font(r, size=size, bold=bold, color=font_rgb or RGBColor(0x21, 0x21, 0x21))
        if fill_rgb is not None:
            cell.fill.solid()
            cell.fill.fore_color.rgb = fill_rgb

    def _style_chart(chart, has_legend=True):
        try:
            chart.has_legend = has_legend
            if has_legend and chart.legend is not None:
                chart.legend.position = XL_LEGEND_POSITION.BOTTOM
                chart.legend.include_in_layout = False
        except Exception as _exc:
            logger.debug("Ignored non-critical exception: %s", _exc)
            pass

    def _add_col(slide, title, cats, series_name, values, ctype=XL_CHART_TYPE.COLUMN_CLUSTERED):
        if not cats or not values:
            return
        _add_title(slide, title, size=15)
        _add_subtitle(slide, "Native chart — double-click to edit data / colours")
        cd = CategoryChartData()
        cd.categories = list(cats)
        cd.add_series(series_name, list(values))
        ch = slide.shapes.add_chart(ctype, Inches(0.4), Inches(1.15), Inches(12.5), Inches(5.8), cd).chart
        _style_chart(ch, has_legend=False)

    def _add_pie(slide, title, cats, values):
        if not cats or not values or sum(values) <= 0:
            return
        _add_title(slide, title, size=15)
        _add_subtitle(slide, "Native pie — every slice is editable")
        cd = CategoryChartData()
        cd.categories = list(cats)
        cd.add_series("Share", list(values))
        ch = slide.shapes.add_chart(XL_CHART_TYPE.PIE, Inches(2.2), Inches(1.15), Inches(8.8), Inches(5.8), cd).chart
        _style_chart(ch, has_legend=True)

    BLUE = RGBColor(0x15, 0x65, 0xC0)
    DARK = RGBColor(0x21, 0x21, 0x21)
    GREY = RGBColor(0x61, 0x61, 0x61)
    GREEN = RGBColor(0x2E, 0x7D, 0x32)
    HEADER_FILL = RGBColor(0x15, 0x65, 0xC0)
    ALT_ROW = RGBColor(0xE3, 0xF2, 0xFD)
    WHITE = RGBColor(0xFF, 0xFF, 0xFF)

    written = 0
    failed = 0

    for gr in gene_results:
        if not gr:
            continue
        gene = gr.get("gene") or "gene"
        folder_label = GENE_FOLDER_NAMES.get(gene, gene)
        safe_name = re.sub(r'[<>:"/\\\\|?*]', "_", str(folder_label)).strip() or gene

        counts = dict(gr.get("mutation_counts") or {})
        if curated_master:
            cset = set(gr.get("curated_set") or set()) | set(curated_master.get(gene, []))
        else:
            cset = set(gr.get("curated_set") or set())
        total = int(gr.get("total_isolates") or 0)
        all_data = gr.get("all_data") or []
        amd = gr.get("all_mutation_data") or []
        rel_counts = {m: c for m, c in counts.items() if m in cset} if cset else {}
        non_counts = {m: c for m, c in counts.items() if m not in cset}
        found = sum(1 for d in all_data if d.get("gene_found"))
        disrupted = sum(1 for d in all_data if d.get("gene_found") and d.get("is_disrupted"))
        missing = sum(1 for d in all_data if not d.get("gene_found"))
        intact = max(0, found - disrupted)
        if total > 0 and amd:
            keep = set(counts.keys())
            n_with = sum(1 for d in amd if set(d.get("mutations") or []) & keep)
            pct_with = 100.0 * n_with / total
        else:
            n_with = 0
            pct_with = 0.0
        top_all = sorted(counts.items(), key=lambda kv: -kv[1])
        top_rel = sorted(rel_counts.items(), key=lambda kv: -kv[1])

        type_tally: Dict[str, int] = {}
        for d in all_data:
            for m in d.get("mutations") or []:
                notation = m.get("mutation") or m.get("notation")
                if not notation or notation not in counts:
                    continue
                t = str(m.get("type") or "substitution").replace("_", " ").title()
                type_tally[t] = type_tally.get(t, 0) + 1
        if not type_tally and counts:
            type_tally = {"Substitution / other": sum(counts.values())}

        cat_bins = {"Rare (<10%)": 0, "Low (10–25%)": 0, "Moderate (25–50%)": 0,
                    "High (50–75%)": 0, "Majority (>75%)": 0}
        if total:
            for _, c in counts.items():
                pct = 100.0 * c / total
                if pct < 10: cat_bins["Rare (<10%)"] += 1
                elif pct < 25: cat_bins["Low (10–25%)"] += 1
                elif pct < 50: cat_bins["Moderate (25–50%)"] += 1
                elif pct < 75: cat_bins["High (50–75%)"] += 1
                else: cat_bins["Majority (>75%)"] += 1

        burdens = [len(d.get("mutations") or []) for d in amd] if amd else []
        burden_hist: Dict[str, int] = {}
        for b in burdens:
            k = str(b)
            burden_hist[k] = burden_hist.get(k, 0) + 1

        try:
            prs = Presentation()
            prs.slide_width = Inches(13.333)
            prs.slide_height = Inches(7.5)
            blank = prs.slide_layouts[6]

            slide = prs.slides.add_slide(blank)
            shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0), Inches(0), Inches(13.333), Inches(1.3))
            shape.fill.solid()
            shape.fill.fore_color.rgb = BLUE
            shape.line.fill.background()
            box = slide.shapes.add_textbox(Inches(0.5), Inches(0.25), Inches(12.3), Inches(0.9))
            tf = box.text_frame
            _style_paragraph(tf.paragraphs[0], f"{gene} — Editable Mutation Figures (OPTION A)",
                             size=22, bold=True, color=WHITE)
            p2 = tf.add_paragraph()
            _style_paragraph(p2, "Native PowerPoint tables & charts only  •  fully editable",
                             size=12, bold=False, color=RGBColor(0xBB, 0xDE, 0xFB))
            info = slide.shapes.add_textbox(Inches(0.7), Inches(1.8), Inches(11.8), Inches(4.5))
            tf = info.text_frame
            tf.word_wrap = True
            for i, line in enumerate([
                f"Gene: {gene}",
                f"Isolates: {total}  |  Distinct: {len(counts)}  |  Events: {sum(counts.values())}",
                f"With ≥1 mut: {n_with} ({pct_with:.1f}%)  |  Curated: {len(rel_counts)}",
                f"Detected/Disrupted/Missing: {found}/{disrupted}/{missing}",
                " ",
                "All slides use native PowerPoint objects. Double-click charts to edit data.",
            ]):
                para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
                _style_paragraph(para, line, size=13, bold=(i == 0), color=DARK)

            def _mut_table(title, items, note, max_rows=25):
                if not items:
                    return
                slide = prs.slides.add_slide(blank)
                _add_title(slide, title)
                _add_subtitle(slide, note)
                headers = ["Rank", "Mutation", "Count", "Prevalence %", "Curated"]
                show = items[:max_rows]
                n_r = len(show) + 1
                tbl = slide.shapes.add_table(n_r, 5, Inches(0.4), Inches(1.15), Inches(12.5),
                                             Inches(min(5.8, 0.27 * n_r + 0.25))).table
                for c, h in enumerate(headers):
                    _fill_table_cell(tbl.cell(0, c), h, bold=True, size=10, fill_rgb=HEADER_FILL, font_rgb=WHITE)
                for i, (mut, cnt) in enumerate(show):
                    pct = (100.0 * cnt / total) if total else 0.0
                    is_cur = "Yes" if mut in cset else "No"
                    fill = ALT_ROW if i % 2 == 0 else WHITE
                    cur_c = GREEN if is_cur == "Yes" else GREY
                    for c, v in enumerate([i + 1, mut, cnt, f"{pct:.1f}%", is_cur]):
                        _fill_table_cell(tbl.cell(i + 1, c), v, bold=(c == 1), size=9, fill_rgb=fill,
                                         font_rgb=(cur_c if c == 4 else DARK))

            slide = prs.slides.add_slide(blank)
            _add_title(slide, f"{gene} — summary metrics")
            metrics = [
                ["Metric", "Value"], ["Gene", gene], ["Isolates", str(total)],
                ["Intact / Disrupted / Missing", f"{intact} / {disrupted} / {missing}"],
                ["Distinct / Events", f"{len(counts)} / {sum(counts.values())}"],
                ["With ≥1 mut", f"{n_with} ({pct_with:.1f}%)"],
                ["Curated / Non-curated", f"{len(rel_counts)} / {len(non_counts)}"],
            ]
            tbl = slide.shapes.add_table(len(metrics), 2, Inches(1.8), Inches(1.2), Inches(9.5),
                                         Inches(min(5.5, 0.4 * len(metrics)))).table
            for ri, row in enumerate(metrics):
                for ci, val in enumerate(row):
                    fill = HEADER_FILL if ri == 0 else (ALT_ROW if ri % 2 == 0 else WHITE)
                    _fill_table_cell(tbl.cell(ri, ci), val, bold=(ri == 0 or ci == 0), size=11,
                                     fill_rgb=fill, font_rgb=(WHITE if ri == 0 else DARK),
                                     align=PP_ALIGN.LEFT if ci == 0 else PP_ALIGN.CENTER)

            if top_all:
                _mut_table(f"{gene} — top mutations", top_all, f"{len(top_all)} distinct", 25)
                for start in range(25, len(top_all), 25):
                    part = top_all[start:start + 25]
                    _mut_table(f"{gene} — catalogue {start+1}–{start+len(part)} of {len(top_all)}", part, "ranked", 25)
            if top_rel:
                _mut_table(f"{gene} — top AMR-relevant", top_rel, f"{len(top_rel)} curated", 20)

            for n, ctype in [(10, XL_CHART_TYPE.COLUMN_CLUSTERED), (15, XL_CHART_TYPE.COLUMN_CLUSTERED),
                             (15, XL_CHART_TYPE.BAR_CLUSTERED), (20, XL_CHART_TYPE.BAR_CLUSTERED),
                             (30, XL_CHART_TYPE.BAR_CLUSTERED)]:
                show = top_all[:n]
                if not show:
                    continue
                slide = prs.slides.add_slide(blank)
                _add_col(slide, f"{gene} — top {n} (count)", [m for m, _ in show], "Isolates",
                         [c for _, c in show], ctype=ctype)
                if total:
                    slide = prs.slides.add_slide(blank)
                    _add_col(slide, f"{gene} — top {n} (prevalence %)", [m for m, _ in show], "Prevalence %",
                             [round(100.0 * c / total, 1) for _, c in show], ctype=ctype)

            if top_all:
                top8 = top_all[:8]
                other = sum(c for _, c in top_all[8:])
                cats = [m for m, _ in top8] + (["Other"] if other else [])
                vals = [c for _, c in top8] + ([other] if other else [])
                slide = prs.slides.add_slide(blank)
                _add_pie(slide, f"{gene} — share top 8 + Other", cats, vals)

            status_pairs = [(v, l) for v, l in [(intact, "Intact"), (disrupted, "Disrupted"), (missing, "Missing")] if v > 0]
            if status_pairs:
                slide = prs.slides.add_slide(blank)
                _add_pie(slide, f"{gene} — detection status", [l for _, l in status_pairs], [v for v, _ in status_pairs])
            if total:
                slide = prs.slides.add_slide(blank)
                _add_pie(slide, f"{gene} — with vs without mutation",
                         ["≥1 mutation", "No mutation"], [n_with, max(0, total - n_with)])

            if type_tally:
                items_t = sorted(type_tally.items(), key=lambda kv: -kv[1])
                slide = prs.slides.add_slide(blank)
                _add_pie(slide, f"{gene} — mutation types", [k for k, _ in items_t], [v for _, v in items_t])
                slide = prs.slides.add_slide(blank)
                _add_col(slide, f"{gene} — type counts", [k for k, _ in items_t], "Events", [v for _, v in items_t])

            cat_items = [(k, v) for k, v in cat_bins.items() if v > 0]
            if cat_items:
                slide = prs.slides.add_slide(blank)
                _add_col(slide, f"{gene} — prevalence categories", [k for k, _ in cat_items],
                         "Distinct", [v for _, v in cat_items])

            cur_e = sum(rel_counts.values())
            non_e = sum(non_counts.values())
            if cur_e or non_e:
                slide = prs.slides.add_slide(blank)
                _add_pie(slide, f"{gene} — curated vs non-curated", ["AMR-relevant", "Non-curated"], [cur_e, non_e])

            if burden_hist:
                keys = sorted(burden_hist.keys(), key=lambda x: int(x) if str(x).isdigit() else 0)
                slide = prs.slides.add_slide(blank)
                _add_col(slide, f"{gene} — burden per isolate", keys, "Isolates", [burden_hist[k] for k in keys])

            if total and counts:
                ys = []
                for thr in thr_list:
                    min_c = int(total * thr / 100.0)
                    common = sorted(((m, c) for m, c in counts.items() if c >= min_c), key=lambda kv: -kv[1])
                    ys.append(len(common))
                    if common:
                        _mut_table(f"{gene} — ≥{thr}%", common, f"min count {min_c}", 20)
                        show = common[:15]
                        slide = prs.slides.add_slide(blank)
                        _add_col(slide, f"{gene} — ≥{thr}% prevalence", [m for m, _ in show], "Prevalence %",
                                 [round(100.0 * c / total, 1) for _, c in show], ctype=XL_CHART_TYPE.BAR_CLUSTERED)
                slide = prs.slides.add_slide(blank)
                _add_col(slide, f"{gene} — threshold trend", [f"≥{t}%" for t in thr_list], "Distinct", ys)

            slide = prs.slides.add_slide(blank)
            _add_title(slide, f"{gene} — how to edit", size=20)
            tips = slide.shapes.add_textbox(Inches(0.7), Inches(1.2), Inches(11.8), Inches(5))
            tf = tips.text_frame
            tf.word_wrap = True
            for i, line in enumerate([
                "Tables: click any cell and type.",
                "Charts: double-click → Edit Data; right-click series for colours.",
                "Titles: click text boxes and type.",
                f"Gene-specific native deck for: {gene}",
            ]):
                para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
                _style_paragraph(para, line, size=13, bold=False, color=DARK)

            out_path = os.path.join(out_dir, f"{safe_name}_Editable.pptx")
            prs.save(out_path)
            written += 1
            logger.info(f"  [Gene-PPTX] {gene:16s} → {os.path.basename(out_path)}  ({len(prs.slides)} slides)")
        except Exception as e:
            failed += 1
            logger.error(f"  [Gene-PPTX] {gene}: failed — {e}")

    logger.info(f"\n  [Gene-PPTX] Folder: {out_dir}")
    logger.error(f"         Written: {written} gene decks  |  Failed: {failed}")
    return out_dir if written else None

_DEFAULT_THESIS_DRUG_CLASS = {
    "penA": "β-lactams (cephalosporins / penicillin)",
    "ponA1": "β-lactams (penicillin)",
    "porB1b": "β-lactams / permeability",
    "blaTEM": "β-lactams (penicillinase)",
    "gyrA": "Fluoroquinolones",
    "parC": "Fluoroquinolones",
    "parE": "Fluoroquinolones",
    "mtrR": "Efflux (multi-drug)",
    "mtrR_promoter": "Efflux (multi-drug)",
    "mtrC": "Efflux (multi-drug)",
    "mtrD": "Efflux (multi-drug)",
    "macAB_promoter": "Efflux (macrolides)",
    "23S_rDNA": "Macrolides (azithromycin)",
    "rpsJ": "Tetracycline",
    "rpsE": "Spectinomycin",
    "folP": "Sulphonamides",
    "ermA": "Macrolides",
    "ermB": "Macrolides",
    "ermC": "Macrolides",
    "ermF": "Macrolides",
    "ereA": "Macrolides",
    "ereB": "Macrolides",
    "mefA": "Macrolides",
    "16S_rDNA": "Aminoglycosides / spectinomycin context",
    "rpoB": "Rifampicin context",
    "rpoD": "Transcription context",
    "rplD": "Macrolides context",
    "rplV": "Macrolides context",
}

_THESIS_DRUG_CLASS = (
    _DOMAIN_CONFIG["drug_class_map"] if _DOMAIN_CONFIG else _DEFAULT_THESIS_DRUG_CLASS
)

def _thesis_safe_savefig(fig, path, dpi=200):
    try:
        fig.tight_layout()
    except Exception as _exc:
        logger.debug("Ignored non-critical exception: %s", _exc)
        pass
    try:
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    except Exception as e:
        logger.error(f"    [Final_for_Thesis] savefig failed {os.path.basename(path)}: {e}")
    try:
        plt.close(fig)
    except Exception as _exc:
        logger.debug("Ignored non-critical exception: %s", _exc)
        pass

def _thesis_pub_rc():
    if not MATPLOTLIB_OK:
        return
    try:
        plt.rcParams.update({
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "figure.dpi": 120,
            "savefig.dpi": 200,
            "axes.spines.top": False,
            "axes.spines.right": False,
        })
    except Exception as _exc:
        logger.debug("Ignored non-critical exception: %s", _exc)
        pass

def _thesis_copy_master_matrix(output_root: str, thesis_root: str) -> Optional[str]:
    src = os.path.join(output_root, "NG_AMR_Mutation_MIC_Matrix.xlsx")
    if not os.path.isfile(src):return None
    dst = os.path.join(thesis_root, "00_NG_AMR_Master_Mutation_Matrix.xlsx")
    try:
        import shutil
        shutil.copy2(src, dst)
        logger.info(f"  [Final_for_Thesis] Copied master matrix → {os.path.basename(dst)}")
        return dst
    except Exception as e:
        logger.error(f"  [Final_for_Thesis] Copy master matrix failed: {e}")
        return None

def _thesis_aggregate_rows(gene_results, curated_master=None):
    rows = []
    for gr in gene_results or []:
        if not gr:
            continue
        gene = gr.get("gene") or "?"
        counts = dict(gr.get("mutation_counts") or {})
        if curated_master:
            cset = set(gr.get("curated_set") or set()) | set(curated_master.get(gene, []))
        else:
            cset = set(gr.get("curated_set") or set())
        total = int(gr.get("total_isolates") or 0)
        all_data = gr.get("all_data") or []
        amd = gr.get("all_mutation_data") or []
        events = sum(counts.values())
        distinct = len(counts)
        if total > 0 and amd:
            keep = set(counts.keys())
            n_with = sum(1 for d in amd if set(d.get("mutations") or []) & keep)
            pct_with = 100.0 * n_with / total
        else:
            n_with, pct_with = 0, 0.0
        found = sum(1 for d in all_data if d.get("gene_found"))
        disrupted = sum(1 for d in all_data if d.get("gene_found") and d.get("is_disrupted"))
        missing = sum(1 for d in all_data if not d.get("gene_found"))
        intact = max(0, found - disrupted)
        rel = {m: c for m, c in counts.items() if m in cset}
        non = {m: c for m, c in counts.items() if m not in cset}
        top = sorted(counts.items(), key=lambda kv: -kv[1])
        burdens = [len(d.get("mutations") or []) for d in amd] if amd else []
        rows.append({
            "gene": gene, "counts": counts, "total": total, "events": events,
            "distinct": distinct, "n_with": n_with, "pct_with": pct_with,
            "intact": intact, "disrupted": disrupted, "missing": missing,
            "curated_set": cset, "rel_counts": rel, "non_counts": non,
            "top": top, "amd": amd, "all_data": all_data, "burdens": burdens,
            "mean_burden": (events / total) if total else 0.0,
            "drug_class": _THESIS_DRUG_CLASS.get(gene, "Other / context"),
        })
    return rows

def _thesis_write_mutation_catalogue(rows, thesis_root, thresholds):
    out_dir = os.path.join(thesis_root, "01_Mutation_Catalogue")
    os.makedirs(out_dir, exist_ok=True)
    cat_rows = []
    for r in rows:
        gene, total, cset = r["gene"], r["total"], r["curated_set"]
        for rank, (mut, cnt) in enumerate(r["top"], 1):
            pct = (100.0 * cnt / total) if total else 0.0
            cat_rows.append({
                "Gene": gene,
                "Drug_class": r["drug_class"],
                "Rank": rank,
                "Mutation": mut,
                "Count": cnt,
                "Prevalence_pct": round(pct, 2),
                "Curated_AMR_relevant": "Yes" if mut in cset else "No",
                "Isolates_analysed": total,
            })
    csv_path = os.path.join(out_dir, "Mutation_Catalogue_All.csv")
    if cat_rows:
        keys = list(cat_rows[0].keys())
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(cat_rows)
    xlsx_path = os.path.join(out_dir, "Mutation_Catalogue_All.xlsx")
    if pd is not None and cat_rows:
        try:
            df = pd.DataFrame(cat_rows)
            with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
                df.to_excel(writer, sheet_name="All_mutations", index=False)
                df[df["Curated_AMR_relevant"] == "Yes"].to_excel(
                    writer, sheet_name="AMR_relevant_only", index=False)
                summ = []
                for r in rows:
                    summ.append({
                        "Gene": r["gene"], "Drug_class": r["drug_class"],
                        "Isolates": r["total"], "Distinct_mutations": r["distinct"],
                        "Events": r["events"], "Pct_with_ge1": round(r["pct_with"], 2),
                        "Curated_distinct": len(r["rel_counts"]),
                        "Non_curated_distinct": len(r["non_counts"]),
                        "Intact": r["intact"], "Disrupted": r["disrupted"], "Missing": r["missing"],
                    })
                pd.DataFrame(summ).to_excel(writer, sheet_name="Gene_summary", index=False)
        except Exception as e:
            logger.error(f"    catalogue xlsx failed: {e}")
    logger.info(f"  [Final_for_Thesis] Mutation catalogue → {out_dir}")
    return out_dir

def _thesis_burden_and_carriage(rows, thesis_root):
    out_dir = os.path.join(thesis_root, "02_Mutation_Burden")
    os.makedirs(out_dir, exist_ok=True)
    if not MATPLOTLIB_OK:
        return out_dir
    _thesis_pub_rc()

    genes = [r["gene"] for r in rows]
    mean_b = [r["mean_burden"] for r in rows]
    pct_w = [r["pct_with"] for r in rows]
    fig, ax = plt.subplots(figsize=(max(8, 0.35 * len(genes)), 4.8))
    ax.bar(genes, mean_b, color="#1565c0", edgecolor="white")
    ax.set_ylabel("Mean mutations / isolate")
    ax.set_title("Mutation burden by gene (mean mutations per isolate)")
    plt.setp(ax.get_xticklabels(), rotation=55, ha="right", fontsize=7)
    _thesis_safe_savefig(fig, os.path.join(out_dir, "01_mean_burden_by_gene.png"))

    fig, ax = plt.subplots(figsize=(max(8, 0.35 * len(genes)), 4.8))
    ax.bar(genes, pct_w, color="#2e7d32", edgecolor="white")
    ax.set_ylabel("% isolates with ≥1 mutation")
    ax.set_ylim(0, 105)
    ax.set_title("Isolate-level carriage of ≥1 mutation by gene")
    plt.setp(ax.get_xticklabels(), rotation=55, ha="right", fontsize=7)
    _thesis_safe_savefig(fig, os.path.join(out_dir, "02_pct_isolates_with_mutation.png"))

    all_burdens = []
    for r in rows:
        all_burdens.extend(r["burdens"])
    if all_burdens:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        maxb = max(all_burdens)
        ax.hist(all_burdens, bins=range(0, maxb + 2), color="#6a1b9a", edgecolor="white", align="left")
        ax.set_xlabel("Mutations per isolate (per gene analysis)")
        ax.set_ylabel("Count")
        ax.set_title("Pooled mutation-burden distribution")
        _thesis_safe_savefig(fig, os.path.join(out_dir, "03_pooled_burden_histogram.png"))

    summ_path = os.path.join(out_dir, "Burden_summary.csv")
    with open(summ_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Gene", "Isolates", "Mean_burden", "Pct_with_ge1", "Distinct", "Events"])
        for r in rows:
            w.writerow([r["gene"], r["total"], f"{r['mean_burden']:.3f}",
                        f"{r['pct_with']:.2f}", r["distinct"], r["events"]])
    logger.info(f"  [Final_for_Thesis] Burden analysis → {out_dir}")
    return out_dir

def _thesis_cooccurrence(rows, thesis_root):
    out_dir = os.path.join(thesis_root, "03_CoOccurrence")
    os.makedirs(out_dir, exist_ok=True)
    if not MATPLOTLIB_OK:
        return out_dir
    _thesis_pub_rc()

    genes = [r["gene"] for r in rows if r["events"] > 0]
    events = [r["events"] for r in rows if r["events"] > 0]
    if genes:
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        colors = plt.cm.tab20.colors
        ax.pie(events, labels=genes, autopct="%1.0f%%", colors=colors[:len(genes)],
               textprops={"fontsize": 7}, startangle=90)
        ax.set_title("Share of isolate–mutation events by gene")
        _thesis_safe_savefig(fig, os.path.join(out_dir, "01_events_share_pie.png"))

    for r in rows:
        top = r["top"][:8]
        amd = r["amd"] or []
        if len(top) < 2 or not amd:
            continue
        labels = [m for m, _ in top]
        n = len(labels)
        mat = [[0] * n for _ in range(n)]
        idx = {m: i for i, m in enumerate(labels)}
        for d in amd:
            present = [m for m in d.get("mutations") or [] if m in idx]
            for i, a in enumerate(present):
                for b in present[i:]:
                    mat[idx[a]][idx[b]] += 1
                    if a != b:
                        mat[idx[b]][idx[a]] += 1
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        im = ax.imshow(mat, cmap="YlOrRd", aspect="auto")
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=7)
        ax.set_yticklabels(labels, fontsize=7)
        for i in range(n):
            for j in range(n):
                ax.text(j, i, str(mat[i][j]), ha="center", va="center", fontsize=7)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Co-occurrence count")
        ax.set_title(f"{r['gene']} — top-mutation co-occurrence")
        _thesis_safe_savefig(fig, os.path.join(out_dir, f"02_{r['gene']}_cooccurrence.png"))
    logger.info(f"  [Final_for_Thesis] Co-occurrence → {out_dir}")
    return out_dir

def _thesis_thresholds(rows, thesis_root, thresholds):
    out_dir = os.path.join(thesis_root, "04_Threshold_Stratification")
    os.makedirs(out_dir, exist_ok=True)
    thr_list = thresholds or [25, 30, 35, 40, 50]
    table_rows = []
    for r in rows:
        total = r["total"]
        if not total:
            continue
        for thr in thr_list:
            min_c = int(total * thr / 100.0)
            common = [(m, c) for m, c in r["counts"].items() if c >= min_c]
            table_rows.append({
                "Gene": r["gene"], "Threshold_pct": thr, "Min_count": min_c,
                "N_mutations_meeting": len(common),
                "Mutations": "; ".join(f"{m}({c})" for m, c in sorted(common, key=lambda x: -x[1])[:20]),
            })
    csv_path = os.path.join(out_dir, "Threshold_mutations.csv")
    if table_rows:
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(table_rows[0].keys()))
            w.writeheader()
            w.writerows(table_rows)
    if MATPLOTLIB_OK and rows:
        _thesis_pub_rc()
        genes = [r["gene"] for r in rows]
        for thr in thr_list:
            vals = []
            for r in rows:
                total = r["total"]
                if not total:
                    vals.append(0)
                    continue
                min_c = int(total * thr / 100.0)
                vals.append(sum(1 for c in r["counts"].values() if c >= min_c))
            fig, ax = plt.subplots(figsize=(max(8, 0.35 * len(genes)), 4.5))
            ax.bar(genes, vals, color="#e65100", edgecolor="white")
            ax.set_ylabel(f"Distinct mutations ≥ {thr}%")
            ax.set_title(f"Mutations meeting ≥{thr}% prevalence threshold")
            plt.setp(ax.get_xticklabels(), rotation=55, ha="right", fontsize=7)
            _thesis_safe_savefig(fig, os.path.join(out_dir, f"threshold_ge{thr}pct.png"))
    logger.info(f"  [Final_for_Thesis] Threshold stratification → {out_dir}")
    return out_dir

def _thesis_curated_vs_non(rows, thesis_root):
    out_dir = os.path.join(thesis_root, "05_Curated_vs_NonCurated")
    os.makedirs(out_dir, exist_ok=True)
    summ = []
    for r in rows:
        ce = sum(r["rel_counts"].values())
        ne = sum(r["non_counts"].values())
        summ.append({
            "Gene": r["gene"], "Curated_distinct": len(r["rel_counts"]),
            "Non_curated_distinct": len(r["non_counts"]),
            "Curated_events": ce, "Non_curated_events": ne,
        })
    with open(os.path.join(out_dir, "Curated_vs_NonCurated.csv"), "w", newline="", encoding="utf-8") as fh:
        if summ:
            w = csv.DictWriter(fh, fieldnames=list(summ[0].keys()))
            w.writeheader()
            w.writerows(summ)
    if MATPLOTLIB_OK and rows:
        _thesis_pub_rc()
        genes = [r["gene"] for r in rows]
        cur_d = [len(r["rel_counts"]) for r in rows]
        non_d = [len(r["non_counts"]) for r in rows]
        import numpy as np
        x = np.arange(len(genes))
        w = 0.38
        fig, ax = plt.subplots(figsize=(max(9, 0.4 * len(genes)), 5))
        ax.bar(x - w / 2, cur_d, w, label="AMR-relevant (curated)", color="#2e7d32")
        ax.bar(x + w / 2, non_d, w, label="Non-curated", color="#e65100")
        ax.set_xticks(x)
        ax.set_xticklabels(genes, rotation=55, ha="right", fontsize=7)
        ax.set_ylabel("Distinct mutations")
        ax.set_title("Curated (AMR-relevant) vs non-curated mutations")
        ax.legend(frameon=False)
        _thesis_safe_savefig(fig, os.path.join(out_dir, "01_curated_vs_non_distinct.png"))
    logger.info(f"  [Final_for_Thesis] Curated vs non-curated → {out_dir}")
    return out_dir

def _thesis_drug_class_summary(rows, thesis_root):
    out_dir = os.path.join(thesis_root, "06_DrugClass_Summary")
    os.makedirs(out_dir, exist_ok=True)
    by_class: Dict[str, List] = {}
    for r in rows:
        by_class.setdefault(r["drug_class"], []).append(r)
    table = []
    for cls, lst in sorted(by_class.items()):
        table.append({
            "Drug_class": cls,
            "Genes": ", ".join(x["gene"] for x in lst),
            "N_genes": len(lst),
            "Total_distinct_mutations": sum(x["distinct"] for x in lst),
            "Total_events": sum(x["events"] for x in lst),
            "Mean_pct_with_mutation": round(sum(x["pct_with"] for x in lst) / len(lst), 2) if lst else 0,
        })
    with open(os.path.join(out_dir, "DrugClass_summary.csv"), "w", newline="", encoding="utf-8") as fh:
        if table:
            w = csv.DictWriter(fh, fieldnames=list(table[0].keys()))
            w.writeheader()
            w.writerows(table)
    if MATPLOTLIB_OK and table:
        _thesis_pub_rc()
        labels = [t["Drug_class"][:40] for t in table]
        vals = [t["Total_events"] for t in table]
        fig, ax = plt.subplots(figsize=(8, max(3.5, 0.35 * len(labels) + 1)))
        ax.barh(range(len(labels)), vals, color="#1565c0")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Isolate–mutation events")
        ax.set_title("Mutation events by antibiotic drug class")
        _thesis_safe_savefig(fig, os.path.join(out_dir, "01_events_by_drug_class.png"))
    logger.info(f"  [Final_for_Thesis] Drug-class summary → {out_dir}")
    return out_dir

def _thesis_cross_gene_panels(rows, thesis_root):
    out_dir = os.path.join(thesis_root, "07_CrossGene_Panels")
    os.makedirs(out_dir, exist_ok=True)
    if not MATPLOTLIB_OK or not rows:
        return out_dir
    _thesis_pub_rc()
    genes = [r["gene"] for r in rows]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes[0, 0].bar(genes, [r["distinct"] for r in rows], color="#1565c0")
    axes[0, 0].set_title("Distinct mutations")
    plt.setp(axes[0, 0].get_xticklabels(), rotation=55, ha="right", fontsize=6)
    axes[0, 1].bar(genes, [r["events"] for r in rows], color="#e65100")
    axes[0, 1].set_title("Events")
    plt.setp(axes[0, 1].get_xticklabels(), rotation=55, ha="right", fontsize=6)
    axes[1, 0].bar(genes, [r["mean_burden"] for r in rows], color="#2e7d32")
    axes[1, 0].set_title("Mean burden")
    plt.setp(axes[1, 0].get_xticklabels(), rotation=55, ha="right", fontsize=6)
    axes[1, 1].bar(genes, [r["pct_with"] for r in rows], color="#6a1b9a")
    axes[1, 1].set_title("% isolates with ≥1 mut")
    axes[1, 1].set_ylim(0, 105)
    plt.setp(axes[1, 1].get_xticklabels(), rotation=55, ha="right", fontsize=6)
    fig.suptitle("Cross-gene overview panel", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _thesis_safe_savefig(fig, os.path.join(out_dir, "01_cross_gene_overview_panel.png"))

    fig, ax = plt.subplots(figsize=(max(9, 0.38 * len(genes)), 5))
    intact = [r["intact"] for r in rows]
    disrupted = [r["disrupted"] for r in rows]
    missing = [r["missing"] for r in rows]
    ax.bar(genes, intact, label="Intact", color="#2e7d32")
    ax.bar(genes, disrupted, bottom=intact, label="Disrupted", color="#ef6c00")
    bottom2 = [a + b for a, b in zip(intact, disrupted)]
    ax.bar(genes, missing, bottom=bottom2, label="Missing", color="#b0bec5")
    ax.set_ylabel("Isolates")
    ax.set_title("Gene detection status (stacked)")
    ax.legend(frameon=False, fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=55, ha="right", fontsize=7)
    _thesis_safe_savefig(fig, os.path.join(out_dir, "02_gene_status_stacked.png"))
    logger.info(f"  [Final_for_Thesis] Cross-gene panels → {out_dir}")
    return out_dir

def _thesis_key_findings(rows, thesis_root, thresholds):
    out_dir = os.path.join(thesis_root, "08_Key_Findings_and_Index")
    os.makedirs(out_dir, exist_ok=True)
    lines = []
    lines.append("# Final_for_Thesis — Key findings summary")
    lines.append("")
    lines.append("Generated automatically from the NG AMR Master Mutation pipeline.")
    lines.append("Use these statements as a starting point for Results / Discussion; cite your exact n and methods.")
    lines.append("")
    n_genes = len(rows)
    max_iso = max((r["total"] for r in rows), default=0)
    lines.append(f"- Genes analysed: **{n_genes}**")
    lines.append(f"- Maximum isolates across genes: **{max_iso}**")
    lines.append("")
    lines.append("## Gene-level highlights")
    for r in sorted(rows, key=lambda x: -x["events"]):
        lines.append(
            f"- **{r['gene']}** ({r['drug_class']}): "
            f"{r['distinct']} distinct mutations, {r['events']} events, "
            f"{r['pct_with']:.1f}% isolates with ≥1 mutation; "
            f"curated AMR-relevant distinct = {len(r['rel_counts'])}."
        )
        if r["top"][:3]:
            top3 = ", ".join(f"{m} (n={c})" for m, c in r["top"][:3])
            lines.append(f"  - Top mutations: {top3}")
    lines.append("")
    lines.append("## Threshold notes")
    thr_list = thresholds or [25, 30, 35, 40, 50]
    lines.append(f"Prevalence thresholds used: {', '.join(str(t) + '%' for t in thr_list)}.")
    lines.append("See `04_Threshold_Stratification/` for mutations meeting each threshold.")
    lines.append("")
    lines.append("## Folder index")
    lines.append("| Folder | Content |")
    lines.append("|--------|---------|")
    lines.append("| `00_NG_AMR_Master_Mutation_Matrix.xlsx` | Copy of master matrix |")
    lines.append("| `01_Mutation_Catalogue` | Full mutation catalogue (CSV/XLSX) |")
    lines.append("| `02_Mutation_Burden` | Burden distributions and summaries |")
    lines.append("| `03_CoOccurrence` | Within-gene co-occurrence + event share |")
    lines.append("| `04_Threshold_Stratification` | ≥25–50% prevalence strata |")
    lines.append("| `05_Curated_vs_NonCurated` | AMR-relevant vs other mutations |")
    lines.append("| `06_DrugClass_Summary` | Grouped by antibiotic class |")
    lines.append("| `07_CrossGene_Panels` | Overview multi-panel figures |")
    lines.append("| `08_Key_Findings_and_Index` | This narrative + index |")
    lines.append("Method-specific MIC figures are in the run folder when matched measurements exist.")
    lines.append("")
    lines.append("## Limitations (n and scope)")
    lines.append("- Analyses are descriptive / exploratory relative to sample size.")
    lines.append("- Genome-based trees are recorded under 03_Phylogeny when available.")
    lines.append("- MLST, NG-STAR, plasmid typing, and structural modelling")
    lines.append("  require assemblies and external tools; they are not produced by this package.")
    lines.append("- Genotype–phenotype concordance modules run only when a completed MIC workbook is supplied.")
    lines.append("")
    path = os.path.join(out_dir, "KEY_FINDINGS.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    logger.info(f"  [Final_for_Thesis] Key findings → {path}")
    return out_dir

def _thesis_build_editable_pptx(rows, thesis_root, thresholds):
    try:
        from pptx import Presentation
        from pptx.util import Inches, Pt
        from pptx.dml.color import RGBColor
        from pptx.enum.text import PP_ALIGN
        from pptx.chart.data import CategoryChartData
        from pptx.enum.chart import XL_CHART_TYPE
    except ImportError:
        logger.info("  [Final_for_Thesis] python-pptx missing – skip editable PPTX")
        return None

    BLUE = RGBColor(0x15, 0x65, 0xC0)
    GREEN = RGBColor(0x2E, 0x7D, 0x32)
    ORANGE = RGBColor(0xE6, 0x51, 0x00)
    DARK = RGBColor(0x21, 0x21, 0x21)
    WHITE = RGBColor(0xFF, 0xFF, 0xFF)
    GREY = RGBColor(0x75, 0x75, 0x75)
    HEADER_FILL = BLUE
    ALT_ROW = RGBColor(0xE3, 0xF2, 0xFD)

    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    blank = prs.slide_layouts[6]

    def _style_paragraph(para, text, size=14, bold=False, color=DARK):
        para.text = text
        para.font.size = Pt(size)
        para.font.bold = bold
        para.font.color.rgb = color
        para.font.name = "Calibri"

    def _add_title(slide, text, top=0.2, height=0.5, size=18, color=BLUE):
        box = slide.shapes.add_textbox(Inches(0.4), Inches(top), Inches(12.5), Inches(height))
        tf = box.text_frame
        tf.word_wrap = True
        _style_paragraph(tf.paragraphs[0], text, size=size, bold=True, color=color)

    def _add_subtitle(slide, text, top=0.65):
        box = slide.shapes.add_textbox(Inches(0.4), Inches(top), Inches(12.5), Inches(0.35))
        _style_paragraph(box.text_frame.paragraphs[0], text, size=11, bold=False, color=GREY)

    def _fill_table_cell(cell, text, bold=False, size=9, fill_rgb=None, font_rgb=None, align=PP_ALIGN.CENTER):
        cell.text = str(text)
        for p in cell.text_frame.paragraphs:
            p.font.size = Pt(size)
            p.font.bold = bold
            p.font.name = "Calibri"
            p.alignment = align
            if font_rgb is not None:
                p.font.color.rgb = font_rgb
        if fill_rgb is not None:
            from pptx.oxml.ns import nsmap
            from pptx.oxml import parse_xml
            tc = cell._tc
            tcPr = tc.get_or_add_tcPr()
            solid = parse_xml(
                f'<a:solidFill xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
                f'<a:srgbClr val="{fill_rgb[0]:02X}{fill_rgb[1]:02X}{fill_rgb[2]:02X}"/></a:solidFill>'
            )
            try:
                cell.fill.solid()
                cell.fill.fore_color.rgb = fill_rgb
            except Exception as _exc:
                logger.debug("Ignored non-critical exception: %s", _exc)
                pass

    def _add_col(slide, title, cats, series_name, values, ctype=XL_CHART_TYPE.COLUMN_CLUSTERED):
        _add_title(slide, title)
        cd = CategoryChartData()
        cd.categories = cats
        cd.add_series(series_name, values)
        slide.shapes.add_chart(ctype, Inches(0.5), Inches(1.0), Inches(12.3), Inches(5.8), cd)

    def _add_pie(slide, title, cats, vals):
        _add_title(slide, title)
        cd = CategoryChartData()
        cd.categories = cats
        cd.add_series("Share", vals)
        slide.shapes.add_chart(XL_CHART_TYPE.PIE, Inches(2.2), Inches(1.0), Inches(8.8), Inches(5.8), cd)

    slide = prs.slides.add_slide(blank)
    _add_title(slide, "Final_for_Thesis — Editable analysis pack", size=24)
    _add_subtitle(slide, "Native tables & charts from NG AMR Master Mutation pipeline (editable in PowerPoint)")
    tip = slide.shapes.add_textbox(Inches(0.6), Inches(2.0), Inches(12), Inches(4))
    tf = tip.text_frame
    tf.word_wrap = True
    for i, line in enumerate([
        "How to edit: click any table cell to type; double-click charts → Edit Data;",
        "right-click series for colours; click title text boxes to rename.",
        f"Genes included: {len(rows)}",
        "Contents: gene summary • burden • curated vs non-curated • thresholds • drug class • top mutations",
    ]):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        _style_paragraph(p, line, size=13, bold=False, color=DARK)

    slide = prs.slides.add_slide(blank)
    _add_title(slide, "Gene summary metrics")
    headers = ["Gene", "Isolates", "Distinct", "Events", "% ≥1 mut", "Curated", "Non-cur", "Intact", "Disrupted", "Missing"]
    n_r = min(len(rows), 28) + 1
    tbl = slide.shapes.add_table(n_r, len(headers), Inches(0.3), Inches(0.9), Inches(12.7),
                                 Inches(min(6.2, 0.22 * n_r + 0.3))).table
    for c, h in enumerate(headers):
        _fill_table_cell(tbl.cell(0, c), h, bold=True, size=8, fill_rgb=HEADER_FILL, font_rgb=WHITE)
    for i, r in enumerate(rows[:28]):
        vals = [r["gene"], r["total"], r["distinct"], r["events"],
                f"{r['pct_with']:.1f}", len(r["rel_counts"]), len(r["non_counts"]),
                r["intact"], r["disrupted"], r["missing"]]
        fill = ALT_ROW if i % 2 == 0 else WHITE
        for c, v in enumerate(vals):
            _fill_table_cell(tbl.cell(i + 1, c), v, bold=(c == 0), size=8, fill_rgb=fill, font_rgb=DARK)

    genes = [r["gene"] for r in rows]
    slide = prs.slides.add_slide(blank)
    _add_col(slide, "Distinct mutations by gene", genes, "Distinct", [r["distinct"] for r in rows])
    slide = prs.slides.add_slide(blank)
    _add_col(slide, "Mutation events by gene", genes, "Events", [r["events"] for r in rows],
             ctype=XL_CHART_TYPE.BAR_CLUSTERED)
    slide = prs.slides.add_slide(blank)
    _add_col(slide, "Mean mutation burden by gene", genes, "Mean burden",
             [round(r["mean_burden"], 2) for r in rows])
    slide = prs.slides.add_slide(blank)
    _add_col(slide, "% isolates with ≥1 mutation", genes, "Prevalence %",
             [round(r["pct_with"], 1) for r in rows])

    slide = prs.slides.add_slide(blank)
    _add_title(slide, "Curated (AMR-relevant) vs non-curated — distinct counts")
    cd = CategoryChartData()
    cd.categories = genes
    cd.add_series("AMR-relevant", [len(r["rel_counts"]) for r in rows])
    cd.add_series("Non-curated", [len(r["non_counts"]) for r in rows])
    slide.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(0.5), Inches(1.0),
                           Inches(12.3), Inches(5.8), cd)

    by_class = {}
    for r in rows:
        by_class.setdefault(r["drug_class"], []).append(r)
    class_labels = sorted(by_class.keys())
    class_events = [sum(x["events"] for x in by_class[c]) for c in class_labels]
    if class_labels:
        slide = prs.slides.add_slide(blank)
        short = [c[:45] for c in class_labels]
        _add_col(slide, "Mutation events by drug class", short, "Events", class_events,
                 ctype=XL_CHART_TYPE.BAR_CLUSTERED)
        slide = prs.slides.add_slide(blank)
        _add_pie(slide, "Share of events by drug class", short, class_events)

    slide = prs.slides.add_slide(blank)
    _add_title(slide, "Gene detection status (Intact / Disrupted / Missing)")
    cd = CategoryChartData()
    cd.categories = genes
    cd.add_series("Intact", [r["intact"] for r in rows])
    cd.add_series("Disrupted", [r["disrupted"] for r in rows])
    cd.add_series("Missing", [r["missing"] for r in rows])
    slide.shapes.add_chart(XL_CHART_TYPE.COLUMN_STACKED, Inches(0.5), Inches(1.0),
                           Inches(12.3), Inches(5.8), cd)

    thr_list = thresholds or [25, 30, 35, 40, 50]
    slide = prs.slides.add_slide(blank)
    _add_title(slide, "Distinct mutations meeting prevalence thresholds (all genes)")
    headers = ["Gene"] + [f"≥{t}%" for t in thr_list]
    n_show = min(len(rows), 25)
    tbl = slide.shapes.add_table(n_show + 1, len(headers), Inches(0.4), Inches(0.95),
                                 Inches(12.5), Inches(min(6.0, 0.22 * (n_show + 1) + 0.3))).table
    for c, h in enumerate(headers):
        _fill_table_cell(tbl.cell(0, c), h, bold=True, size=9, fill_rgb=HEADER_FILL, font_rgb=WHITE)
    for i, r in enumerate(rows[:n_show]):
        total = r["total"] or 1
        rowv = [r["gene"]]
        for thr in thr_list:
            min_c = int(total * thr / 100.0)
            rowv.append(sum(1 for c in r["counts"].values() if c >= min_c))
        fill = ALT_ROW if i % 2 == 0 else WHITE
        for c, v in enumerate(rowv):
            _fill_table_cell(tbl.cell(i + 1, c), v, bold=(c == 0), size=8, fill_rgb=fill)

    flat = []
    for r in rows:
        for mut, cnt in r["top"][:5]:
            flat.append((r["gene"], mut, cnt, r["total"], mut in r["curated_set"]))
    flat.sort(key=lambda x: -x[2])
    flat = flat[:30]
    if flat:
        slide = prs.slides.add_slide(blank)
        _add_title(slide, "Top mutations across genes (up to 5 per gene, ranked by count)")
        headers = ["Gene", "Mutation", "Count", "Prevalence %", "Curated"]
        tbl = slide.shapes.add_table(len(flat) + 1, 5, Inches(0.4), Inches(0.95),
                                     Inches(12.5), Inches(min(6.0, 0.2 * len(flat) + 0.35))).table
        for c, h in enumerate(headers):
            _fill_table_cell(tbl.cell(0, c), h, bold=True, size=9, fill_rgb=HEADER_FILL, font_rgb=WHITE)
        for i, (gene, mut, cnt, total, is_c) in enumerate(flat):
            pct = round(100.0 * cnt / total, 1) if total else 0
            fill = ALT_ROW if i % 2 == 0 else WHITE
            for c, v in enumerate([gene, mut, cnt, f"{pct}%", "Yes" if is_c else "No"]):
                _fill_table_cell(tbl.cell(i + 1, c), v, bold=(c == 1), size=8, fill_rgb=fill,
                                 font_rgb=(GREEN if c == 4 and is_c else DARK))

    slide = prs.slides.add_slide(blank)
    _add_title(slide, "How to use this deck in your thesis", size=20)
    box = slide.shapes.add_textbox(Inches(0.6), Inches(1.2), Inches(12), Inches(5.5))
    tf = box.text_frame
    tf.word_wrap = True
    for i, line in enumerate([
        "1. Tables: click a cell and type; change fonts via Home ribbon.",
        "2. Charts: double-click → Edit Data in Excel sheet that opens; close to update chart.",
        "3. Colours: right-click a series → Format Data Series → Fill.",
        "4. Titles / labels: click text boxes and edit.",
        "5. Complementary raster figures live in Final_for_Thesis/01–08 folders (PNG).",
        "6. Advanced genomics modules (typing, phylogeny, plasmid, structure) are documented under 10–13.",
    ]):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        _style_paragraph(p, line, size=14, bold=False, color=DARK)

    out_path = os.path.join(thesis_root, "Final_for_Thesis_Editable.pptx")
    prs.save(out_path)
    logger.info(f"  [Final_for_Thesis] Editable PPTX → {os.path.basename(out_path)}  ({len(prs.slides)} slides)")
    return out_path

def _thesis_advanced_genomics_modules(rows, thesis_root, gene_results, output_root):
    d10 = os.path.join(thesis_root, "10_Typing_MLST_NGSTAR")
    os.makedirs(d10, exist_ok=True)
    readme10 = os.path.join(d10, "README_Typing_requirements.md")
    with open(readme10, "w", encoding="utf-8") as fh:
        fh.write("""# Molecular typing (MLST / NG-MAST / NG-STAR)

## Why this folder exists
WHO EGASP and published NG WGS theses report MLST, NG-MAST and NG-STAR alongside AMR determinants.

## What this pipeline cannot invent
Allele calls require assemblies + tools (e.g. mlst, pyngSTar, Pathogenwatch, PubMLST).

## Recommended workflow (run outside or extend this script)
1. Assemblies: same FASTAs used for gene extraction.
2. MLST: `mlst --scheme neisseria your.fasta` or Pathogenwatch.
3. NG-STAR: pyngSTar or Pathogenwatch (penA, mtrR, porB, ponA, gyrA, parC, 23S).
4. NG-MAST: porB + tbpB alleles.

## Outputs to place here after external typing
- `MLST_profiles.csv` — Isolate, ST, allele profile
- `NGSTAR_profiles.csv` — Isolate, NG-STAR ST / CC
- `NGMAST_profiles.csv` — Isolate, NG-MAST ST
- Figures: ST diversity bar, ST vs resistance phenotype overlay

## Template table
See `Typing_template.csv` in this folder.
""")
    with open(os.path.join(d10, "Typing_template.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Isolate_ID", "MLST_ST", "NG_STAR_ST", "NG_MAST_ST", "Notes"])
        isolates = []
        for r in rows:
            for d in (r.get("all_data") or []):
                iso = d.get("isolate") or d.get("isolate_id") or d.get("sample")
                if iso and iso not in isolates:
                    isolates.append(iso)
            if isolates:
                break
        for iso in isolates:
            w.writerow([iso, "", "", "", "fill after external typing"])
    logger.info(f"  [Final_for_Thesis] Typing module (templates) → {d10}")

    d11 = os.path.join(thesis_root, "11_Phylogeny")
    os.makedirs(d11, exist_ok=True)
    with open(os.path.join(d11, "README_Phylogeny_requirements.md"), "w", encoding="utf-8") as fh:
        fh.write("""# Core-genome / SNP phylogeny

## Scientific role
Annotate tree with MIC phenotype, key determinants, and typing STs so the tree answers a biological question (not decoration).

## Recommended tools
- Core alignment: Roary / Panaroo / Snippy / Parsnp
- Recombination filter (optional): Gubbins
- Tree: IQ-TREE or RAxML-NG
- Visualisation: iTOL / ggtree / Microreact

## Outputs to place here
- `core_snp.tree` or Newick
- Annotated figure PNG/PDF with tip colours for CIP-R, mosaic penA, ST, etc.

## Minimal question checklist
- Do CIP-R isolates cluster?
- Does an NG-STAR lineage match a phenotype?
- Are genotype–phenotype discordant isolates phylogenetically distinct?
""")
    logger.info(f"  [Final_for_Thesis] Phylogeny module (checklist) → {d11}")

    d12 = os.path.join(thesis_root, "12_Plasmid_MGE")
    os.makedirs(d12, exist_ok=True)
    with open(os.path.join(d12, "README_Plasmid_requirements.md"), "w", encoding="utf-8") as fh:
        fh.write("""# Plasmid / mobile genetic element analysis

## Relevance for NG AMR
blaTEM (PPNG) and tetM (TRNG) are often plasmid-borne; classify Asia/Africa/Toronto/Rio blaTEM types and Dutch/American tetM when possible.

## Tools
- PlasmidFinder / ABRicate against plasmid DB
- MOB-suite (optional)
- Check gene presence from this pipeline: blaTEM, tetM (if in gene list)

## Data from this pipeline already useful
- Presence/absence of blaTEM from gene_results (see mutation catalogue / status tables)
""")
    plasmid_rows = []
    for r in rows:
        if r["gene"] in ("blaTEM", "tetM"):
            plasmid_rows.append({
                "Gene": r["gene"],
                "Isolates": r["total"],
                "Found_intact": r["intact"],
                "Disrupted": r["disrupted"],
                "Missing": r["missing"],
                "Pct_found": round(100.0 * (r["intact"] + r["disrupted"]) / r["total"], 1) if r["total"] else 0,
            })
    if plasmid_rows:
        with open(os.path.join(d12, "blaTEM_tetM_detection_from_pipeline.csv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(plasmid_rows[0].keys()))
            w.writeheader()
            w.writerows(plasmid_rows)
    logger.info(f"  [Final_for_Thesis] Plasmid module → {d12}")

    d13 = os.path.join(thesis_root, "13_Structural_penA_PBP2")
    os.makedirs(d13, exist_ok=True)
    PBP2_MOTIFS = {
        "SxxK_active_site": list(range(310, 320)),
        "SxN_motif": list(range(360, 366)),
        "KTG_motif": list(range(497, 505)),
        "beta3_beta4_loop": list(range(500, 520)),
        "alpha2_helix_cluster": list(range(480, 500)),
        "Ala501_cluster": [501],
        "Gly545_region": list(range(542, 548)),
    }
    with open(os.path.join(d13, "README_Structural_mapping.md"), "w", encoding="utf-8") as fh:
        fh.write("""# penA / PBP2 structural context (literature-anchored)

## What this module does without AlphaFold
Maps amino-acid substitutions detected by this pipeline onto published PBP2 structural motifs
(active-site SxxK, SxN, KTG, β3–β4 loop, Ala501 cluster, mosaic-associated positions).

## What full structural modelling would add
AlphaFold/ColabFold models of each isolate's PBP2 with mutations highlighted — run externally
and drop PNG/PDB into this folder.

## Key literature themes
- Ala501 substitutions near the active site increase rigidity associated with ESC reduced susceptibility
- Mosaic penA alleles accumulate changes in α2 / β3–β4 loop and active-site-interacting residues
""")
    penA_rows = [r for r in rows if r["gene"] == "penA"]
    mapped = []
    if penA_rows:
        r = penA_rows[0]
        for mut, cnt in r["top"]:
            m = re.match(r"^([A-Z*])(\d+)([A-Z*]|Ter)?", str(mut))
            pos = int(m.group(2)) if m else None
            motifs = []
            if pos is not None:
                for name, positions in PBP2_MOTIFS.items():
                    if pos in positions:
                        motifs.append(name)
            mapped.append({
                "Mutation": mut,
                "Count": cnt,
                "Prevalence_pct": round(100.0 * cnt / r["total"], 2) if r["total"] else 0,
                "AA_position": pos if pos is not None else "",
                "Structural_motif": "; ".join(motifs) if motifs else "outside_mapped_motifs",
                "Curated": "Yes" if mut in r["curated_set"] else "No",
            })
        with open(os.path.join(d13, "penA_mutations_structural_motif_map.csv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(mapped[0].keys()))
            w.writeheader()
            w.writerows(mapped)
        if MATPLOTLIB_OK and mapped:
            from collections import Counter
            motif_counts = Counter()
            for row in mapped:
                if row["Structural_motif"] == "outside_mapped_motifs":
                    motif_counts["outside_mapped_motifs"] += 1
                else:
                    for mname in row["Structural_motif"].split("; "):
                        motif_counts[mname] += 1
            labels = list(motif_counts.keys())
            vals = [motif_counts[k] for k in labels]
            fig, ax = plt.subplots(figsize=(8, 4.5))
            ax.barh(labels, vals, color="#6a1b9a")
            ax.set_xlabel("Distinct penA mutations mapping to motif")
            ax.set_title("penA mutations vs published PBP2 structural motifs")
            _thesis_safe_savefig(fig, os.path.join(d13, "penA_motif_mapping_counts.png"))
    logger.info(f"  [Final_for_Thesis] Structural penA mapping → {d13}")
    return d10, d11, d12, d13

def generate_final_for_thesis_package(
    gene_results: List[dict],
    output_root: str,
    curated_master: Optional[dict] = None,
    thresholds: Optional[List[int]] = None,
    mic_excel: Optional[str] = None,
) -> str:
    logger.info("\n" + "#" * 70)
    logger.info("# FINAL_FOR_THESIS PACKAGE (scientist-level thesis outputs + editable PPTX)")
    logger.info("#" * 70)
    thesis_root = os.path.join(output_root, "Final_for_Thesis")
    os.makedirs(thesis_root, exist_ok=True)

    rows = _thesis_aggregate_rows(gene_results, curated_master)
    if not rows:
        logger.info("  [Final_for_Thesis] No gene_results — skip")
        return thesis_root

    _thesis_copy_master_matrix(output_root, thesis_root)
    _thesis_write_mutation_catalogue(rows, thesis_root, thresholds)
    _thesis_burden_and_carriage(rows, thesis_root)
    _thesis_cooccurrence(rows, thesis_root)
    _thesis_thresholds(rows, thesis_root, thresholds)
    _thesis_curated_vs_non(rows, thesis_root)
    _thesis_drug_class_summary(rows, thesis_root)
    _thesis_cross_gene_panels(rows, thesis_root)
    _thesis_key_findings(rows, thesis_root, thresholds)


    try:
        _thesis_build_editable_pptx(rows, thesis_root, thresholds)
    except Exception as e:
        logger.error(f"  [Final_for_Thesis] Editable PPTX failed: {e}")

    try:
        _thesis_advanced_genomics_modules(rows, thesis_root, gene_results, output_root)
    except Exception as e:
        logger.error(f"  [Final_for_Thesis] Advanced genomics modules failed: {e}")

    try:
        idx_path = os.path.join(thesis_root, "08_Key_Findings_and_Index", "KEY_FINDINGS.md")
        extra = (
            "\n\n## Also in this package\n"
            "- `Final_for_Thesis_Editable.pptx` — native editable tables/charts\n"
            "- `10_Typing_MLST_NGSTAR/` — templates + checklist for external typing\n"
            "- `11_Phylogeny/` — phylogeny checklist\n"
            "- `12_Plasmid_MGE/` — plasmid module + blaTEM/tetM detection if analysed\n"
            "- `13_Structural_penA_PBP2/` — penA mutations mapped to PBP2 motifs\n"
        )
        if os.path.isfile(idx_path):
            with open(idx_path, "a", encoding="utf-8") as fh:
                fh.write(extra)
    except Exception as _exc:
        logger.debug("Ignored non-critical exception: %s", _exc)
        pass

    logger.info("\nFinal_for_Thesis complete → %s", thesis_root)
    logger.info("  Includes: analyses 01–09 + Final_for_Thesis_Editable.pptx + modules 10–13")
    return thesis_root

def generate_extended_thesis_publication_package(gene_results: List[dict], output_root: str,
                                                  curated_master: Optional[dict] = None,
                                                  thresholds: Optional[List[int]] = None) -> None:

    logger.info("\n" + "#" * 70)
    logger.info("# EXTENDED THESIS PUBLICATION PACKAGE (≥120 styles + master PPTX + gene PPTX)")
    logger.info("#" * 70)
    generate_extended_master_publication_figures(gene_results, output_root, curated_master)
    generate_extended_gene_wise_images(gene_results, output_root, curated_master, thresholds)
    build_master_figures_pptx(output_root, gene_results)
    build_master_figures_pptx_native_editable(output_root, gene_results, curated_master)
    build_gene_specific_editable_pptx_folder(
        output_root, gene_results, curated_master, thresholds
    )
    logger.info("\nExtended package complete. New folders / files:")
    logger.info(f"  {os.path.join(output_root, 'Extended_Publication_Figures', 'All_Mutations')}")
    logger.info(f"  {os.path.join(output_root, 'Extended_Publication_Figures', 'NG_AMR_Relevant')}")
    logger.info(f"  {os.path.join(output_root, 'Extended_Gene_Figures', '<gene>', 'All_Mutations')}")
    logger.info(f"  {os.path.join(output_root, 'Extended_Gene_Figures', '<gene>', 'NG_AMR_Relevant')}")
    logger.info(f"  {os.path.join(output_root, 'NG_AMR_Master_Mutation_Analysis_Figures.pptx')}")
    logger.info(f"  {os.path.join(output_root, 'NG_AMR_Master_Mutation_Analysis_Figures_Editable.pptx')}")
    logger.info(f"  {os.path.join(output_root, 'Gene_Editable_PPTX')}  (one editable PPTX per gene)")

def run_self_tests() -> Optional[bool]:
    failures = []
    checks_run = [0]

    def check(label, condition):
        checks_run[0] += 1
        if not condition:
            failures.append(label)
            logger.error("SELF-TEST FAILED: %s", label)
        else:
            logger.info("SELF-TEST OK: %s", label)

    try:
        p, lo, hi = _wilson_ci(0, 0)
        check("_wilson_ci(0,0) returns zeros", (p, lo, hi) == (0.0, 0.0, 0.0))
        p, lo, hi = _wilson_ci(50, 100)
        check("_wilson_ci(50,100) point estimate is 0.5", abs(p - 0.5) < 1e-9)
        check("_wilson_ci(50,100) CI bounds ordered and within [0,1]", 0.0 <= lo <= p <= hi <= 1.0)
        p, lo, hi = _wilson_ci(100, 100)
        check("_wilson_ci(100,100) upper bound <= 1.0", hi <= 1.0 and p == 1.0)
    except Exception as exc:
        failures.append(f"_wilson_ci raised {exc!r}")
        logger.error("SELF-TEST FAILED: _wilson_ci raised %r", exc)

    try:
        parsed = _classify_substitution("S91F")
        check("_classify_substitution parses S91F", parsed is not None and parsed[0] == "S" and parsed[1] == 91 and parsed[2] == "F")
        parsed_cons = _classify_substitution("L99I")
        check("_classify_substitution flags L99I as Conservative", parsed_cons is not None and parsed_cons[3] == "Conservative (same class)")
        parsed_rad = _classify_substitution("D95K")
        check("_classify_substitution flags D95K as Radical", parsed_rad is not None and parsed_rad[3] == "Radical (class change)")
        check("_classify_substitution rejects malformed notation", _classify_substitution("not_a_mutation") is None)
    except Exception as exc:
        failures.append(f"_classify_substitution raised {exc!r}")
        logger.error("SELF-TEST FAILED: _classify_substitution raised %r", exc)

    try:
        check("short_isolate_label truncates long names", len(short_isolate_label("A" * 50, max_len=18)) <= 18)
        check("normalize_sample_id is deterministic", normalize_sample_id("Sample_1") == normalize_sample_id("Sample_1"))
        check("normalize_sample_id lower-cases output", normalize_sample_id("SAMPLE_ABC") == normalize_sample_id("SAMPLE_ABC").lower())
    except Exception as exc:
        failures.append(f"ID normalisation raised {exc!r}")
        logger.error("SELF-TEST FAILED: ID normalisation raised %r", exc)

    try:
        muts = call_protein_mutations("MSTV", "MSTV")
        check("call_protein_mutations finds no changes for identical sequences", muts == [])

        muts = call_protein_mutations("MSTV", "MSTF")
        check("call_protein_mutations detects a simple substitution",
              any(m["type"] == "substitution" and m["notation"] == "V4F" for m in muts))

        muts = call_protein_mutations("MST-V", "MSTAV")
        check("call_protein_mutations detects an insertion", any(m["type"] == "insertion" for m in muts))

        muts = call_protein_mutations("MSTAV", "MST-V")
        check("call_protein_mutations detects a deletion", any(m["type"] == "deletion" for m in muts))

        muts = call_protein_mutations("MSTV", "MST*")
        check("call_protein_mutations detects a stop-gain", any(m["type"] == "stop_gain" for m in muts))
    except Exception as exc:
        failures.append(f"call_protein_mutations raised {exc!r}")
        logger.error("SELF-TEST FAILED: call_protein_mutations raised %r", exc)

    try:
        counts = {"A501V": 5, "X1Y": 2}
        all_counts, rel_counts = _split_relevant("penA", counts, curated_set={"X1Y"})
        check("_split_relevant returns all counts unchanged", all_counts == counts)
        check("_split_relevant keeps curated UNION literature-hotspot mutations only",
              set(rel_counts.keys()) == {"A501V", "X1Y"})
    except Exception as exc:
        failures.append(f"_split_relevant raised {exc!r}")
        logger.error("SELF-TEST FAILED: _split_relevant raised %r", exc)


    check("16S_rDNA and 23S_rDNA are configured as non-CDS",
          get_gene_type("16S_rDNA") != "CDS" and get_gene_type("23S_rDNA") != "CDS")
    check("promoter regions are configured as non-CDS",
          get_gene_type("mtrR_promoter") != "CDS"
          and get_gene_type("macAB_promoter") != "CDS")

    if Seq is not None and PairwiseAligner is not None:
        try:
            ref = "ATGAAACCCGGGTTTAAACCCGGGTAA"
            qry_fs = "ATGAACCCGGGTTTAAACCCGGGTAA"
            nt_fs = analyze_nucleotide(ref, qry_fs)
            prot_fs = analyze_protein("selftest", ref, qry_fs,
                                      nt_mutations=nt_fs["mutations"])
            check("1 bp deletion is called as a frameshift",
                  len(prot_fs["frameshifts"]) == 1
                  and prot_fs["frameshifts"][0]["type"] == "frameshift_deletion")
            check("frameshifted gene is assessed as DISRUPTED",
                  assess_disruption(prot_fs) is True)

            qry_inframe = "ATG" + "CCCGGGTTTAAACCCGGGTAA"
            nt_if = analyze_nucleotide(ref, qry_inframe)
            prot_if = analyze_protein("selftest", ref, qry_inframe,
                                     nt_mutations=nt_if["mutations"])
            check("in-frame 3 bp deletion is NOT called a frameshift",
                  prot_if["frameshifts"] == [])
        except Exception as exc:
            failures.append(f"frameshift regression raised {exc!r}")
            logger.error("SELF-TEST FAILED: frameshift regression raised %r", exc)

        try:
            ref_syn = "ATGAAACCCGGGTTTTAA"
            qry_syn = "ATGAAACCCGGGTTCTAA"
            nt_syn = analyze_nucleotide(ref_syn, qry_syn)
            cat_syn = categorize_nt_mutations(nt_syn["mutations"], ref_syn, qry_syn, True)
            check("synonymous SNP is called synonymous (no indel present)",
                  cat_syn["syn_count"] == 1 and cat_syn["missense_count"] == 0)

            ref_i = "ATGAAACCCGGGTTTTAA"
            qry_i = "ATGAAAGCCCGGGTTCTAA"
            nt_i = analyze_nucleotide(ref_i, qry_i)
            cat_i = categorize_nt_mutations(nt_i["mutations"], ref_i, qry_i, True)
            snp_effects = [m.get("effect") for m in nt_i["mutations"] if m["type"] == "SNP"]
            check("SNP downstream of an indel is not mis-called missense",
                  all(e != "Non-synonymous (missense)" for e in snp_effects))
            check("frame-affected SNPs are labelled, not silently guessed",
                  all((e or "").startswith("SNP (") for e in snp_effects))
            check("categorize_nt_mutations keeps syn + non_syn == snp_count",
                  cat_i["syn_count"] + cat_i["non_syn_count"] == cat_i["snp_count"])
        except Exception as exc:
            failures.append(f"syn/non-syn regression raised {exc!r}")
            logger.error("SELF-TEST FAILED: syn/non-syn regression raised %r", exc)

        try:
            ref_r = "AAACGTACGTAAACGT"
            qry_r = "AAACGTACGTGAACGT"
            nt_r = analyze_nucleotide(ref_r, qry_r)
            notations = [m["notation"] for m in nt_r["mutations"]]
            check("non-coding variants use nucleotide notation (e.g. A11G)",
                  any(re.fullmatch(r"[ACGT]\d+[ACGT]", n) for n in notations))
        except Exception as exc:
            failures.append(f"non-coding regression raised {exc!r}")
            logger.error("SELF-TEST FAILED: non-coding regression raised %r", exc)

        try:
            prot_stop = analyze_protein("selftest", "ATGAAACCCTAA", "ATGTAACCCTAA")
            check("stop-gain is present in the reported mutations",
                  any(m["type"] == "stop_gain" for m in prot_stop["mutations"]))
            check("silent_mutations stays empty (documented invariant)",
                  prot_stop["silent_mutations"] == [])
        except Exception as exc:
            failures.append(f"stop-gain regression raised {exc!r}")
            logger.error("SELF-TEST FAILED: stop-gain regression raised %r", exc)

    try:
        check(">=35% of 24 isolates requires 9 isolates, not 8",
              max(1, int(math.ceil(24 * 0.35))) == 9)
        check(">=25% of 24 isolates requires 6 isolates",
              max(1, int(math.ceil(24 * 0.25))) == 6)
    except Exception as exc:
        failures.append(f"threshold regression raised {exc!r}")
        logger.error("SELF-TEST FAILED: threshold regression raised %r", exc)

    try:
        import ast as _ast
        _lvl = {"debug", "info", "warning", "error", "exception", "critical"}
        _ph = re.compile(r"%(?:\([^)]*\))?[-#0 +]*\d*(?:\.\d+)*[hlL]?[diouxXeEfFgGcrsa%]")

        def _literal_text(node):
            if isinstance(node, _ast.Constant) and isinstance(node.value, str):
                return node.value
            if isinstance(node, _ast.JoinedStr):
                return "".join(v.value for v in node.values
                               if isinstance(v, _ast.Constant)
                               and isinstance(v.value, str))
            if isinstance(node, _ast.BinOp) and isinstance(node.op, _ast.Add):
                left, right = _literal_text(node.left), _literal_text(node.right)
                if left is None or right is None:
                    return None
                return left + right
            return None

        bad_calls = []
        with open(os.path.abspath(__file__), encoding="utf-8") as _fh:
            _tree = _ast.parse(_fh.read())
        for _n in _ast.walk(_tree):
            if not isinstance(_n, _ast.Call) or len(_n.args) < 2:
                continue
            _f = _n.func
            if not (isinstance(_f, _ast.Attribute) and _f.attr in _lvl):
                continue
            if not (isinstance(_f.value, _ast.Name) and _f.value.id == "logger"):
                continue
            _txt = _literal_text(_n.args[0])
            if _txt is None:
                continue
            _found = [x for x in _ph.findall(_txt) if x != "%%"]
            if len(_found) != len(_n.args) - 1:
                bad_calls.append(_n.lineno)
        check("every logger call's args match its %-placeholders",
              not bad_calls)
        if bad_calls:
            logger.error("  offending logger call(s) at line(s): %s", bad_calls)
    except Exception as exc:
        failures.append(f"logger-format audit raised {exc!r}")
        logger.error("SELF-TEST FAILED: logger-format audit raised %r", exc)

    logger.info("-" * 60)
    if failures:
        logger.error("SELF-TEST SUITE: %d of %d check(s) failed: %s",
                     len(failures), checks_run[0], failures)
        return False
    if Seq is None or PairwiseAligner is None:
        logger.warning("%d checks passed; 10 alignment-dependent checks unavailable (Biopython missing).", checks_run[0])
        return None
    logger.info("ALL SELF-TESTS PASSED (%d checks)", checks_run[0])
    return True

GENE_SHEET_ORDER = [
    "penA", "ponA1", "porB1b", "blaTEM",
    "parC", "parE", "gyrA",
    "mtrR", "mtrR_promoter", "mtrC", "mtrD",
    "macAB_promoter",
    "rpsJ", "folP", "rpsE", "16S_rDNA", "rpoB", "rpoD",
    "rplD", "rplV", "23S_rDNA",
    "ermA", "ermB", "ermC", "ermF", "ereA", "ereB", "mefA",
    "norM_promoter", "tetM",
]

def _save_fig(fig, path: str, dpi: int = 350) -> str:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path

def _make_image_pptx(image_paths: List[str], pptx_path: str,
                         title_prefix: str = "") -> Optional[str]:
    try:
        from pptx import Presentation
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN
    except Exception as exc:
        logger.warning("  [PPTX] python-pptx unavailable: %s", exc)
        return None
    if not image_paths:
        return None
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    blank = prs.slide_layouts[6]
    for idx, image_path in enumerate(image_paths, 1):
        slide = prs.slides.add_slide(blank)
        tb = slide.shapes.add_textbox(Inches(0.35), Inches(0.15),
                                      Inches(12.6), Inches(0.45))
        tf = tb.text_frame
        tf.clear()
        p = tf.paragraphs[0]
        p.text = f"{title_prefix} {idx:03d} — {Path(image_path).stem}"
        p.font.size = Pt(18)
        p.font.bold = True
        p.alignment = PP_ALIGN.LEFT
        slide.shapes.add_picture(
            image_path, Inches(0.35), Inches(0.72),
            width=Inches(12.6), height=Inches(6.35),
        )
    os.makedirs(os.path.dirname(pptx_path), exist_ok=True)
    prs.save(pptx_path)
    return pptx_path

def _gene_mutation_dataset(gene_results):
    isolates = set()
    genes = []
    curated = {}
    per_gene = {}
    for gr in gene_results or []:
        gene = str(gr.get("gene", "gene"))
        genes.append(gene)
        curated[gene] = set(gr.get("curated_set") or [])
        per_gene[gene] = {}
        for d in gr.get("all_data") or []:
            isolates.add(str(d.get("isolate", "")))
        for d in gr.get("all_mutation_data") or []:
            iso = str(d.get("isolate", ""))
            isolates.add(iso)
            per_gene[gene].setdefault(iso, set()).update(d.get("mutations") or [])
    genes = list(dict.fromkeys(genes))
    isolates = sorted([x for x in isolates if x], key=lambda s: (int(re.match(r"^\d+", s).group()) if re.match(r"^\d+", s) else 10**9, s))
    mutations = []
    for gene in genes:
        vals = sorted({m for ss in per_gene.get(gene, {}).values() for m in ss}, key=_mutation_sort_key)
        mutations.extend((gene, m) for m in vals)
    mat = np.zeros((len(isolates), len(mutations)), dtype=float)
    iso_idx = {x:i for i,x in enumerate(isolates)}
    for j,(gene,m) in enumerate(mutations):
        for iso, ms in per_gene.get(gene, {}).items():
            if m in ms and iso in iso_idx:
                mat[iso_idx[iso], j] = 1
    gene_burden = np.zeros((len(isolates), len(genes)), dtype=float)
    for gi,g in enumerate(genes):
        for iso,ms in per_gene.get(g, {}).items():
            if iso in iso_idx:
                gene_burden[iso_idx[iso],gi] = len(ms)
    labels = [f"{g}_{m}" for g,m in mutations]
    return isolates, genes, mutations, labels, mat, gene_burden, curated

def generate_exclusive_figures_package(gene_results: List[dict],
                                       output_root: str,
                                       curated_master: Optional[dict] = None,
                                       mic_workbook_path: Optional[str] = None) -> Optional[str]:
    if not MATPLOTLIB_OK:
        logger.warning("  [Exclusive_Figures] matplotlib unavailable — skipped")
        return None
    out_dir = os.path.join(output_root, "Exclusive_Figures")
    os.makedirs(out_dir, exist_ok=True)
    isolates, genes, mutations, labels, mat, gene_burden, curated = _gene_mutation_dataset(gene_results)
    n_iso, n_mut = mat.shape

    if n_iso == 0:
        isolates = ["No isolates"]
        mat = np.zeros((1, max(1, len(genes) or 1)))
        n_iso, n_mut = mat.shape
    if n_mut == 0:
        labels = list(genes) or ["No mutations"]
        mat = np.zeros((n_iso, len(labels)))
        n_mut = len(labels)

    short_iso = [short_isolate_label(x) for x in isolates]
    mut_freq = mat.sum(axis=0).astype(float)
    gene_counts = np.array(
        [gene_burden[:, i].sum() if gene_burden.shape[1] > i else 0 for i in range(len(genes))],
        dtype=float,
    )
    total_burden = mat.sum(axis=1).astype(float)

    TOP_MUT_HEAT = min(40, n_mut)
    TOP_MUT_CORR = min(25, n_mut)
    TOP_MUT_BAR = min(25, n_mut)
    TOP_MUT_ENT = min(30, n_mut)

    order_m_freq = np.argsort(-mut_freq)
    top_m = order_m_freq[:TOP_MUT_HEAT]
    top_c = order_m_freq[:TOP_MUT_CORR]
    top_b = order_m_freq[:TOP_MUT_BAR]
    top_e = order_m_freq[:TOP_MUT_ENT]
    order_i = np.argsort(-total_burden)

    try:
        from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
        from scipy.spatial.distance import pdist
        _SCIPY = True
    except Exception:
        _SCIPY = False

    leaf_iso = order_i
    Z_iso = None
    if _SCIPY and n_iso >= 2 and mat.sum() > 0:
        try:
            dist = np.nan_to_num(pdist(mat, metric="jaccard"))
            if np.all(dist == 0):
                dist = np.nan_to_num(pdist(mat, metric="euclidean"))
            Z_iso = linkage(dist, method="average")
            leaf_iso = leaves_list(Z_iso)
        except Exception:
            Z_iso = None
            leaf_iso = order_i

    def _short_mut(lab, n=22):
        s = str(lab)
        return s if len(s) <= n else (s[: n - 1] + "…")

    paths = []

    def save(fig, family, number):
        path = os.path.join(out_dir, f"{number:03d}_{family}.png")
        paths.append(_save_fig(fig, path, dpi=300))

    CMAPS_HEAT = ["Greys", "Blues", "viridis", "magma", "cividis"]
    CMAPS_CORR = ["RdBu_r", "magma", "coolwarm", "Greens", "Purples"]
    CMAPS_PRES = ["Blues", "Greens", "Purples", "Oranges", "Reds"]
    CMAPS_CLUS = ["viridis", "plasma", "inferno", "cividis", "magma"]
    SOLID = ["#1565c0", "#2E7D32", "#E65100", "#6A1B9A", "#C62828"]
    MARKERS = ["o", "s", "^", "D", "P"]

    for v in range(1, 6):
        sub = mat[order_i][:, top_m]
        fig_w = max(10, min(18, 0.32 * TOP_MUT_HEAT + 4))
        fig_h = max(5.5, min(12, 0.28 * n_iso + 2.5))
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        if v == 1:
            im = ax.imshow(sub, aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
        elif v == 2:
            im = ax.imshow(sub, aspect="auto", interpolation="nearest", cmap="Blues", vmin=0, vmax=1)
        elif v == 3:
            im = ax.imshow(sub, aspect="auto", interpolation="nearest", cmap="viridis", vmin=0, vmax=1)
        elif v == 4:
            im = ax.imshow(sub, aspect="auto", interpolation="nearest", cmap="magma", vmin=0, vmax=1)
        else:
            im = ax.imshow(sub, aspect="auto", interpolation="nearest", cmap="cividis", vmin=0, vmax=1)
        ax.set_yticks(range(n_iso))
        ax.set_yticklabels([short_iso[i] for i in order_i], fontsize=7)
        ax.set_xticks(range(len(top_m)))
        ax.set_xticklabels([_short_mut(labels[i]) for i in top_m], rotation=90, fontsize=7)
        ax.set_title(f"Mutation landscape — style {v} (top {len(top_m)} by frequency)", fontweight="bold")
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Presence")
        fig.tight_layout()
        save(fig, "Mutation_Landscape", 1 + (v - 1) * 25)

        k = TOP_MUT_BAR
        oo = order_m_freq[:k][::-1]
        pct = mut_freq[oo] / max(1, n_iso) * 100
        fig, ax = plt.subplots(figsize=(10, max(5.5, 0.28 * k + 2)))
        if v == 1:
            ax.barh(np.arange(k), pct, color=SOLID[0], edgecolor="white")
        elif v == 2:
            ax.barh(np.arange(k), pct, color=plt.cm.Blues(0.35 + 0.6 * pct / 100), edgecolor="white")
        elif v == 3:
            ax.hlines(np.arange(k), 0, pct, color="#90A4AE", linewidth=1.4)
            ax.scatter(pct, np.arange(k), s=55, color=SOLID[2], zorder=3, edgecolors="white")
        elif v == 4:
            ax.barh(np.arange(k), pct, color=SOLID[3], edgecolor="white", height=0.7)
            for yi, p in enumerate(pct):
                ax.text(p + 1.5, yi, f"{p:.0f}%", va="center", fontsize=7)
        else:
            ax.fill_betweenx(np.arange(k), 0, pct, alpha=0.35, color=SOLID[4])
            ax.plot(pct, np.arange(k), marker="o", color=SOLID[4], linewidth=1.6)
        ax.set_yticks(np.arange(k))
        ax.set_yticklabels([_short_mut(labels[i], 28) for i in oo], fontsize=8)
        ax.set_xlabel("Prevalence (%)")
        ax.set_xlim(0, 110)
        ax.set_title(f"Mutation prevalence — style {v} (top {k})", fontweight="bold")
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        save(fig, "Mutation_Prevalence", 2 + (v - 1) * 25)

        if genes:
            oo = np.argsort(-gene_counts)
            fig, ax = plt.subplots(figsize=(11, 6))
            vals = gene_counts[oo]
            if v == 1:
                ax.bar(np.arange(len(genes)), vals, color=SOLID[0], edgecolor="white")
            elif v == 2:
                ax.barh(np.arange(len(genes)), vals[::-1], color=SOLID[1], edgecolor="white")
                ax.set_yticks(np.arange(len(genes)))
                ax.set_yticklabels([genes[i] for i in oo[::-1]], fontsize=8)
                ax.set_xlabel("Mutation events")
                ax.set_title(f"Gene-level mutation burden — style {v}", fontweight="bold")
                ax.grid(axis="x", alpha=0.25)
                fig.tight_layout()
                save(fig, "Gene_Burden", 3 + (v - 1) * 25)
                genes_done = True
            elif v == 3:
                ax.plot(np.arange(len(genes)), vals, marker="o", linewidth=2, color=SOLID[2])
                ax.fill_between(np.arange(len(genes)), vals, alpha=0.2, color=SOLID[2])
            elif v == 4:
                ax.stem(np.arange(len(genes)), vals, linefmt="C3-", markerfmt="C3o", basefmt=" ")
            else:
                ax.bar(np.arange(len(genes)), vals, color=plt.cm.plasma(np.linspace(0.2, 0.85, len(genes))), edgecolor="white")
            if v != 2:
                ax.set_xticks(np.arange(len(genes)))
                ax.set_xticklabels([genes[i] for i in oo], rotation=70, fontsize=8)
                ax.set_ylabel("Mutation events")
                ax.set_title(f"Gene-level mutation burden — style {v}", fontweight="bold")
                ax.grid(axis="y", alpha=0.25)
                fig.tight_layout()
                save(fig, "Gene_Burden", 3 + (v - 1) * 25)

        fig, ax = plt.subplots(figsize=(9, 6))
        if v == 1:
            ax.hist(total_burden, bins=max(5, min(15, int(np.sqrt(max(1, n_iso))) + 2)),
                    color=SOLID[0], edgecolor="white")
            ax.set_xlabel("Mutations / isolate"); ax.set_ylabel("No. isolates")
        elif v == 2:
            ax.boxplot(total_burden, orientation="vertical", patch_artist=True,
                       boxprops=dict(facecolor="#90CAF9"))
            ax.set_ylabel("Mutations / isolate")
        elif v == 3:
            ax.violinplot(total_burden, showmeans=True, showmedians=True)
            ax.set_ylabel("Mutations / isolate")
        elif v == 4:
            vals = np.sort(total_burden)
            ax.plot(np.arange(1, len(vals) + 1), vals, marker="o", markersize=5, color=SOLID[1])
            ax.set_xlabel("Isolate rank"); ax.set_ylabel("Mutations / isolate")
        else:
            if np.all(total_burden >= 0) and len(total_burden):
                counts = np.bincount(total_burden.astype(int))
                ax.bar(np.arange(len(counts)), counts, color=SOLID[3], edgecolor="white")
            ax.set_xlabel("Mutations / isolate"); ax.set_ylabel("No. isolates")
        ax.set_title(f"Isolate burden — style {v}", fontweight="bold")
        fig.tight_layout()
        save(fig, "Isolate_Burden", 4 + (v - 1) * 25)

        sub_mat = mat[:, top_c]
        if sub_mat.shape[1] > 1:
            corr = np.corrcoef(sub_mat.T)
        else:
            corr = np.array([[1.0]])
        corr = np.nan_to_num(corr)
        if v == 1:
            arr = corr
        elif v == 2:
            arr = np.abs(corr)
        elif v == 3:
            arr = np.sign(corr)
        elif v == 4:
            arr = (corr > 0.25).astype(float)
        else:
            arr = (corr > 0.50).astype(float)
        kk = arr.shape[0]
        fig, ax = plt.subplots(figsize=(max(8, 0.35 * kk + 3), max(7, 0.35 * kk + 2.5)))
        im = ax.imshow(
            arr, cmap=CMAPS_CORR[v - 1],
            vmin=(-1 if v in (1, 3) else 0), vmax=1, aspect="auto",
        )
        tick_labs = [_short_mut(labels[i], 18) for i in top_c]
        ax.set_xticks(range(kk)); ax.set_xticklabels(tick_labs, rotation=90, fontsize=7)
        ax.set_yticks(range(kk)); ax.set_yticklabels(tick_labs, fontsize=7)
        ax.set_title(f"Mutation co-occurrence map — style {v} (top {kk})", fontweight="bold")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        save(fig, "Mutation_Cooccurrence", 5 + (v - 1) * 25)

    family_specs = [
        ("Gene_Mutation_Stacked", "stacked"),
        ("Mutation_Rank", "rank"),
        ("Gene_Presence_Heatmap", "presence"),
        ("Isolate_Gene_Bubble", "bubble"),
        ("Gene_Mutation_Scatter", "scatter"),
        ("Prevalence_Dotplot", "dot"),
        ("Mutation_Area", "area"),
        ("Gene_Radar", "radar"),
        ("Burden_ECDF", "ecdf"),
        ("Mutation_Box_By_Gene", "boxgene"),
        ("SNP_vs_Indel", "snpindel"),
        ("Curated_Noncurated", "curated"),
        ("Gene_Entropy", "entropy"),
        ("Mutation_Lollipop", "lollipop"),
        ("Pairwise_Burden", "pair"),
        ("Cumulative_Mutation", "cumulative"),
        ("Gene_Profile", "profile"),
        ("Mutation_Cluster", "cluster"),
        ("Gene_Mutation_DotMatrix", "dotmatrix"),
        ("Publication_Dashboard", "dashboard"),
    ]

    for fi, (family, kind) in enumerate(family_specs, 6):
        for v in range(1, 6):
            num = (fi - 6) * 5 + v + 25
            fig = plt.figure(figsize=(11, 7))
            ax = fig.add_subplot(111)

            if kind == "stacked":
                x = np.arange(len(genes))
                base = np.zeros(len(genes))
                n_show = min(n_iso, 8)
                for q in range(n_show):
                    vals = gene_burden[order_i[q], :] if gene_burden.shape[1] else np.zeros(len(genes))
                    ax.bar(x, vals, bottom=base, label=short_iso[order_i[q]], width=0.8)
                    base += vals
                ax.set_xticks(x)
                ax.set_xticklabels(genes, rotation=70, fontsize=7)
                ax.legend(fontsize=7, ncol=2, frameon=False)
                ax.set_ylabel("Mutation count")

            elif kind == "rank":
                vals = np.sort(mut_freq)[::-1]
                if v % 2 == 0:
                    vals = np.log1p(vals)
                ax.plot(
                    np.arange(1, len(vals) + 1), vals,
                    marker=MARKERS[v - 1], markersize=4, linewidth=1.5, color=SOLID[v - 1],
                )
                ax.set_xlabel("Mutation rank")
                ax.set_ylabel("Frequency" if v % 2 else "log1p(frequency)")
                ax.grid(alpha=0.25)

            elif kind == "presence":
                vals = (gene_burden > 0).astype(int) if gene_burden.size else np.zeros((n_iso, max(1, len(genes))))
                ord_p = np.argsort(-vals.sum(axis=1))
                ax.imshow(
                    vals[ord_p], aspect="auto",
                    cmap=CMAPS_PRES[v - 1],
                    interpolation="nearest", vmin=0, vmax=1,
                )
                ax.set_xticks(range(len(genes)))
                ax.set_xticklabels(genes, rotation=70, fontsize=7)
                ax.set_yticks(range(n_iso))
                ax.set_yticklabels([short_iso[i] for i in ord_p], fontsize=6)

            elif kind == "bubble":
                xg = np.repeat(np.arange(len(genes)), n_iso)
                yi = np.tile(np.arange(n_iso), len(genes))
                gb = gene_burden.T.flatten() if gene_burden.size else np.zeros_like(xg, dtype=float)
                if v == 1:
                    scale = 20
                    ax.scatter(xg, yi, s=np.clip(gb, 0, None) * scale + 10, alpha=0.5,
                               marker="o", c=gb, cmap="Blues", edgecolors="none")
                elif v == 2:
                    msk = gb > 0
                    ax.scatter(xg[msk], yi[msk], s=np.clip(gb[msk], 0, None) * 22 + 14, alpha=0.65,
                               marker="s", c=gb[msk], cmap="Greens", edgecolors="white", linewidths=0.3)
                elif v == 3:
                    msk = gb > 0
                    ax.scatter(xg[msk], yi[msk], s=np.clip(gb[msk], 0, None) * 18 + 16, alpha=0.7,
                               marker="^", c=gb[msk], cmap="Oranges", edgecolors="none")
                elif v == 4:
                    means = gene_burden.mean(axis=0) if gene_burden.size else np.zeros(len(genes))
                    ax.scatter(np.arange(len(genes)), means, s=80 + 40 * means, alpha=0.8,
                               marker="D", c=means, cmap="Purples", edgecolors="black", linewidths=0.4)
                    ax.set_xticks(range(len(genes)))
                    ax.set_xticklabels(genes, rotation=70, fontsize=7)
                    ax.set_ylabel("Mean mutations / isolate")
                else:
                    msk = gb > 0
                    ax.scatter(xg[msk], yi[msk], s=np.clip(gb[msk], 0, None) * 16 + 12, alpha=0.55,
                               marker="P", c=gb[msk], cmap="Reds", edgecolors="none")
                if v != 4:
                    ax.set_xticks(range(len(genes)))
                    ax.set_xticklabels(genes, rotation=70, fontsize=7)
                    ax.set_yticks(range(n_iso))
                    ax.set_yticklabels(short_iso, fontsize=6)
                    ax.set_xlim(-0.6, max(0, len(genes) - 0.4))
                    ax.set_ylim(-0.6, max(0, n_iso - 0.4))

            elif kind == "scatter":
                yv = mut_freq / max(1, n_iso)
                if v == 1:
                    xv = np.arange(n_mut)
                    yv_plot = yv
                    size_src = mut_freq
                elif v == 2:
                    ord_s = np.argsort(-yv)
                    xv = np.arange(n_mut)
                    yv_plot = yv[ord_s]
                    size_src = mut_freq[ord_s]
                elif v == 3:
                    ord_s = order_m_freq[:TOP_MUT_BAR]
                    xv = np.arange(len(ord_s))
                    yv_plot = yv[ord_s]
                    size_src = mut_freq[ord_s]
                elif v == 4:
                    ord_s = order_m_freq[:TOP_MUT_BAR]
                    xv = mut_freq[ord_s]
                    yv_plot = yv[ord_s]
                    size_src = mut_freq[ord_s]
                else:
                    ord_s = order_m_freq[:TOP_MUT_BAR]
                    xv = np.arange(len(ord_s))
                    yv_plot = np.log1p(mut_freq[ord_s])
                    size_src = mut_freq[ord_s]
                sizes = 20 + 8 * (size_src / max(1, size_src.max())) * (10 + 4 * v)
                sc = ax.scatter(
                    xv, yv_plot, s=sizes, c=yv_plot, cmap="viridis",
                    alpha=0.75, edgecolors="white", linewidths=0.3, marker=MARKERS[v - 1],
                )
                ax.set_xlabel({1: "Mutation index", 2: "Mutation rank", 3: "Top mutation rank",
                               4: "Isolate count", 5: "Top mutation rank"}[v])
                ax.set_ylabel("Prevalence" if v != 5 else "log1p(count)")
                if v != 5:
                    ax.set_ylim(-0.02, 1.08)
                top5 = np.argsort(-yv_plot)[:5]
                for ti in top5:
                    ax.annotate(f"{float(yv_plot[ti]):.0%}" if v != 5 else f"{yv_plot[ti]:.1f}",
                                (xv[ti], yv_plot[ti]),
                                textcoords="offset points", xytext=(0, 6),
                                ha="center", fontsize=7, color="#333333")
                fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)

            elif kind == "dot":
                oo = order_m_freq[: min(n_mut, 25)][::-1]
                ax.scatter(
                    mut_freq[oo] / max(1, n_iso) * 100, np.arange(len(oo)),
                    s=40 + mut_freq[oo] * 6, marker=MARKERS[v - 1],
                    color=SOLID[v - 1], alpha=0.85, edgecolors="white",
                )
                ax.set_yticks(np.arange(len(oo)))
                ax.set_yticklabels([_short_mut(labels[i], 28) for i in oo], fontsize=7)
                ax.set_xlabel("Prevalence (%)")
                ax.set_xlim(0, 105)
                ax.grid(axis="x", alpha=0.25)

            elif kind == "area":
                vals = np.sort(total_burden)
                if len(vals):
                    ax.fill_between(np.arange(len(vals)), vals, alpha=0.35, color=SOLID[v - 1])
                    ax.plot(vals, linewidth=1.8, color=SOLID[v - 1])
                ax.set_xlabel("Ordered isolates")
                ax.set_ylabel("Mutation burden")
                ax.grid(alpha=0.25)

            elif kind == "radar":
                if genes:
                    N = min(len(genes), 12)
                    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
                    vals = gene_counts[:N]
                    vals = vals / max(1.0, float(vals.max()) if len(vals) else 1.0)
                    angles += angles[:1]
                    vals = np.r_[vals, vals[:1]]
                    ax.remove()
                    ax = fig.add_subplot(111, polar=True)
                    ax.plot(angles, vals, linewidth=2, color=SOLID[v - 1])
                    ax.fill(angles, vals, alpha=0.25, color=SOLID[v - 1])
                    ax.set_xticks(angles[:-1])
                    ax.set_xticklabels(genes[:N], fontsize=7)
                    ax.set_ylim(0, 1)

            elif kind == "ecdf":
                vals = np.sort(total_burden)
                yy = np.arange(1, len(vals) + 1) / max(1, len(vals))
                ax.step(vals, yy, where="post", linewidth=2, color=SOLID[v - 1])
                ax.set_xlabel("Mutation burden")
                ax.set_ylabel("ECDF")
                ax.grid(alpha=0.25)

            elif kind == "boxgene":
                vals = [gene_burden[:, i] for i in range(min(len(genes), 15))]
                if vals:
                    ax.boxplot(vals, showfliers=False, patch_artist=True,
                               boxprops=dict(facecolor=["#BBDEFB", "#C8E6C9", "#FFE0B2", "#E1BEE7", "#FFCDD2"][v - 1]))
                    ax.set_xticks(range(1, len(vals) + 1))
                    ax.set_xticklabels(genes[: len(vals)], rotation=70, fontsize=7)
                    ax.set_ylabel("Mutations / isolate")

            elif kind == "snpindel":
                snp = np.array([
                    sum(1 for m in (gene_results[j].get("mutation_counts") or {}) if re.search(r"\d", str(m)))
                    for j in range(len(gene_results))
                ], dtype=float)
                indel = np.array([
                    sum(1 for m in (gene_results[j].get("mutation_counts") or {})
                        if any(x in str(m).lower() for x in ("ins", "del")))
                    for j in range(len(gene_results))
                ], dtype=float)
                x = np.arange(len(snp)); w = 0.38
                ax.bar(x - w / 2, snp, width=w, label="SNP-like", color=SOLID[0])
                ax.bar(x + w / 2, indel, width=w, label="Indel-like", color=SOLID[2])
                ax.set_xticks(x)
                ax.set_xticklabels(genes[: len(x)], rotation=70, fontsize=7)
                ax.legend(frameon=False)
                ax.set_ylabel("Count")

            elif kind == "curated":
                curvals = []; nonvals = []
                for g in genes:
                    vals = set()
                    for gr in gene_results:
                        if gr.get("gene") == g:
                            for d in gr.get("all_mutation_data") or []:
                                vals.update(d.get("mutations") or [])
                    curvals.append(sum(1 for m in vals if m in curated.get(g, set())))
                    nonvals.append(max(0, len(vals) - curvals[-1]))
                x = np.arange(len(genes))
                ax.bar(x, curvals, label="Curated", color="#2E7D32")
                ax.bar(x, nonvals, bottom=curvals, label="Non-curated", color="#E65100")
                ax.set_xticks(x)
                ax.set_xticklabels(genes, rotation=70, fontsize=7)
                ax.legend(frameon=False)
                ax.set_ylabel("Distinct mutations")

            elif kind == "entropy":
                oo = top_e[::-1]
                ent = []
                for j in oo:
                    p0 = mut_freq[j] / max(1, n_iso)
                    ent.append(-(p0 * np.log2(max(p0, 1e-12)) + (1 - p0) * np.log2(max(1 - p0, 1e-12))))
                ax.barh(np.arange(len(oo)), ent, color=SOLID[v - 1], edgecolor="white")
                ax.set_yticks(np.arange(len(oo)))
                ax.set_yticklabels([_short_mut(labels[i], 26) for i in oo], fontsize=7)
                ax.set_xlabel("Binary entropy")
                ax.grid(axis="x", alpha=0.25)

            elif kind == "lollipop":
                oo = order_m_freq[: min(n_mut, 25)][::-1]
                yy = np.arange(len(oo))
                ax.hlines(yy, 0, mut_freq[oo], linewidth=1.2, color="#90A4AE")
                ax.scatter(mut_freq[oo], yy, s=50 + 8 * v, color=SOLID[v - 1], zorder=3, marker=MARKERS[v - 1])
                ax.set_yticks(yy)
                ax.set_yticklabels([_short_mut(labels[i], 26) for i in oo], fontsize=7)
                ax.set_xlabel("Isolate count")
                ax.grid(axis="x", alpha=0.25)

            elif kind == "pair":
                if len(total_burden):
                    y = np.roll(total_burden, v)
                    ax.scatter(total_burden, y, alpha=0.75, s=40, color=SOLID[v - 1], edgecolors="white")
                    lo = min(total_burden.min(), y.min())
                    hi = max(total_burden.max(), y.max())
                    ax.plot([lo, hi], [lo, hi], linestyle="--", color="grey")
                    ax.set_xlabel("Burden")
                    ax.set_ylabel("Shifted burden")
                    ax.grid(alpha=0.25)

            elif kind == "cumulative":
                vals = np.sort(mut_freq)[::-1]
                cs = np.cumsum(vals)
                ax.plot(np.arange(1, len(cs) + 1), cs, linewidth=2, color=SOLID[v - 1])
                ax.fill_between(np.arange(1, len(cs) + 1), cs, alpha=0.2, color=SOLID[v - 1])
                ax.set_xlabel("Mutation rank")
                ax.set_ylabel("Cumulative frequency")
                ax.grid(alpha=0.25)

            elif kind == "profile":
                if gene_burden.size:
                    if v == 1:
                        profile = gene_burden.mean(axis=0); ylab = "Mean mutations"
                    elif v == 2:
                        profile = np.median(gene_burden, axis=0); ylab = "Median mutations"
                    elif v == 3:
                        profile = np.max(gene_burden, axis=0); ylab = "Max mutations"
                    elif v == 4:
                        profile = np.std(gene_burden, axis=0); ylab = "SD mutations"
                    else:
                        profile = (gene_burden > 0).mean(axis=0) * 100; ylab = "% isolates with ≥1 mut"
                    ax.plot(np.arange(len(profile)), profile, marker="o", linewidth=1.8, color=SOLID[v - 1])
                    ax.set_xticks(range(len(genes)))
                    ax.set_xticklabels(genes, rotation=70, fontsize=7)
                    ax.set_ylabel(ylab)
                    ax.grid(alpha=0.25)

            elif kind == "cluster":
                if v == 1 and Z_iso is not None:
                    ax.remove()
                    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 3.5], wspace=0.05)
                    axd = fig.add_subplot(gs[0, 0])
                    axh = fig.add_subplot(gs[0, 1])
                    dendrogram(Z_iso, ax=axd, orientation="left", labels=short_iso,
                               color_threshold=0, above_threshold_color="#555555", no_labels=False)
                    axd.tick_params(axis="y", labelsize=6)
                    axd.set_xticks([])
                    sub = mat[leaf_iso][:, top_m]
                    axh.imshow(sub, aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
                    axh.set_yticks(range(n_iso))
                    axh.set_yticklabels([short_iso[i] for i in leaf_iso], fontsize=6)
                    axh.set_xticks(range(len(top_m)))
                    axh.set_xticklabels([_short_mut(labels[i]) for i in top_m], rotation=90, fontsize=7)
                    axh.set_title("Clustered isolates × top mutations", fontsize=10)
                    ax = axh
                elif v == 2:
                    sub = mat[order_i][:, top_m]
                    rgb = np.zeros((n_iso, len(top_m), 3))
                    rgb[:, :] = (0.05, 0.05, 0.05)
                    rgb[sub > 0] = (0.95, 0.78, 0.10)
                    ax.imshow(rgb, aspect="auto", interpolation="nearest")
                    ax.set_yticks(range(n_iso))
                    ax.set_yticklabels([short_iso[i] for i in order_i], fontsize=6)
                    ax.set_xticks(range(len(top_m)))
                    ax.set_xticklabels([_short_mut(labels[i]) for i in top_m], rotation=90, fontsize=7)
                elif v == 3:
                    sub = mat[leaf_iso][:, top_m]
                    ax.imshow(sub, aspect="auto", interpolation="nearest", cmap="plasma", vmin=0, vmax=1)
                    ax.set_yticks(range(n_iso))
                    ax.set_yticklabels([short_iso[i] for i in leaf_iso], fontsize=6)
                    ax.set_xticks(range(len(top_m)))
                    ax.set_xticklabels([_short_mut(labels[i]) for i in top_m], rotation=90, fontsize=7)
                elif v == 4:
                    counts = mat[order_i][:, top_m].sum(axis=1)
                    ax.hlines(np.arange(n_iso), 0, counts, color="#90A4AE", linewidth=1)
                    ax.scatter(counts, np.arange(n_iso), s=36, color=SOLID[3], zorder=3)
                    ax.set_yticks(range(n_iso))
                    ax.set_yticklabels([short_iso[i] for i in order_i], fontsize=6)
                    ax.set_xlabel(f"# of top-{len(top_m)} mutations present")
                else:
                    sub = mat[order_i][:, top_m]
                    ax.imshow(sub, aspect="auto", interpolation="nearest", cmap="magma", vmin=0, vmax=1)
                    ax.set_yticks(range(n_iso))
                    ax.set_yticklabels([short_iso[i] for i in order_i], fontsize=6)
                    ax.set_xticks(range(len(top_m)))
                    ax.set_xticklabels([_short_mut(labels[i]) for i in top_m], rotation=90, fontsize=7)

            elif kind == "dotmatrix":
                sub = mat[:n_iso, top_m]
                rows, cols = np.where(sub > 0)
                ax.scatter(
                    cols, rows, s=28 + 6 * v,
                    marker=MARKERS[v - 1],
                    color=SOLID[v - 1], alpha=0.8, edgecolors="none",
                )
                ax.set_xticks(range(len(top_m)))
                ax.set_xticklabels([_short_mut(labels[i]) for i in top_m], rotation=90, fontsize=7)
                ax.set_yticks(range(n_iso))
                ax.set_yticklabels(short_iso[:n_iso], fontsize=6)
                ax.set_xlim(-0.5, len(top_m) - 0.5)
                ax.set_ylim(-0.5, n_iso - 0.5)
                ax.set_xlabel("Mutations"); ax.set_ylabel("Isolates")
                ax.invert_yaxis()

            elif kind == "dashboard":
                ax.axis("off")
                txt = (
                    f"EXCLUSIVE FIGURE {num:03d}\n\n"
                    f"Genes: {len(genes)}\nIsolates: {n_iso}\n"
                    f"Distinct mutation labels: {n_mut}\n"
                    f"Total mutation observations: {int(mat.sum())}\n"
                    f"Style variant: {v}\n\n"
                    f"Dense axes capped: heatmaps ≤{TOP_MUT_HEAT}, "
                    f"co-occurrence ≤{TOP_MUT_CORR}"
                )
                ax.text(
                    0.5, 0.55, txt, ha="center", va="center", fontsize=16,
                    fontweight="bold", transform=ax.transAxes,
                    family="sans-serif",
                )

            ax.set_title(f"{family.replace('_', ' ')} — variant {v}", fontweight="bold")
            fig.tight_layout()
            save(fig, family, num)

    next_num = 126

    for v in range(1, 6):
        fig = plt.figure(figsize=(12, max(6, 0.28 * n_iso + 2)))
        if v == 1 and Z_iso is not None:
            ax = fig.add_subplot(111)
            dendrogram(
                Z_iso, ax=ax, orientation="left",
                labels=short_iso,
                color_threshold=0.7 * max(Z_iso[:, 2]) if len(Z_iso) else 0,
                leaf_font_size=7,
            )
            ax.set_xlabel("Distance (Jaccard / Euclidean)")
            ax.set_title("Hierarchical clustering of mutation profiles (style 1)", fontweight="bold")
        elif v == 2 and Z_iso is not None:
            ax = fig.add_subplot(111, polar=True)
            n = n_iso
            angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
            radii = 0.4 + 0.6 * (total_burden[leaf_iso] / max(1.0, total_burden.max()))
            ax.scatter(angles, radii, c=total_burden[leaf_iso], cmap="viridis", s=50, edgecolors="white")
            for a, lab in zip(angles, [short_iso[i] for i in leaf_iso]):
                ax.text(a, 1.05, lab, fontsize=6, ha="center", va="center")
            ax.set_yticklabels([])
            ax.set_title("Mutation-profile radial map — burden-coloured tips (style 2)", fontweight="bold", pad=18)
        elif v == 3:
            ax = fig.add_subplot(111)
            vals = total_burden[leaf_iso]
            colors = plt.cm.viridis(vals / max(1.0, vals.max() if len(vals) else 1.0))
            ax.barh(np.arange(n_iso), vals, color=colors, edgecolor="white")
            ax.set_yticks(range(n_iso))
            ax.set_yticklabels([short_iso[i] for i in leaf_iso], fontsize=7)
            ax.set_xlabel("Mutation burden")
            ax.set_title("Mutation-profile cluster order × burden (style 3)", fontweight="bold")
        elif v == 4 and Z_iso is not None:
            ax = fig.add_subplot(111)
            dendrogram(
                Z_iso, ax=ax, orientation="top",
                labels=short_iso,
                leaf_rotation=90, leaf_font_size=7,
                color_threshold=0.5 * max(Z_iso[:, 2]) if len(Z_iso) else 0,
            )
            ax.set_ylabel("Distance")
            ax.set_title("Mutation-profile dendrogram — top orientation (style 4)", fontweight="bold")
        else:
            ax = fig.add_subplot(111)
            idx = leaf_iso[: min(n_iso, 30)]
            sub = mat[idx]
            if sub.shape[0] > 1:
                c = np.corrcoef(sub)
                c = np.nan_to_num(c)
            else:
                c = np.array([[1.0]])
            im = ax.imshow(c, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
            labs = [short_iso[i] for i in idx]
            ax.set_xticks(range(len(idx))); ax.set_xticklabels(labs, rotation=90, fontsize=6)
            ax.set_yticks(range(len(idx))); ax.set_yticklabels(labs, fontsize=6)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Profile correlation")
            ax.set_title("Mutation-profile heatmap — isolate profile correlation (style 5)", fontweight="bold")
        fig.tight_layout()
        save(fig, "Phylogenetic_Tree", next_num)
        next_num += 1

    for v in range(1, 6):
        fig = plt.figure(figsize=(11, 8))
        ax = fig.add_subplot(111)
        n_g = max(1, len(genes))
        gene_bin = (gene_burden > 0).astype(int) if gene_burden.size else np.zeros((n_iso, n_g))
        co = gene_bin.T @ gene_bin
        angles = np.linspace(0, 2 * np.pi, n_g, endpoint=False)
        xs = np.cos(angles)
        ys = np.sin(angles)
        max_co = float(co.max()) if co.size else 1.0
        for i in range(n_g):
            for j in range(i + 1, n_g):
                w = co[i, j]
                if w <= 0:
                    continue
                lw = 0.5 + 4.0 * (w / max(1.0, max_co))
                if v == 1:
                    ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color="#1565c0", alpha=0.35 + 0.5 * (w / max(1.0, max_co)), lw=lw, zorder=1)
                elif v == 2:
                    if w / max(1.0, max_co) >= 0.3:
                        ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color="#6A1B9A", alpha=0.7, lw=lw, zorder=1)
                elif v == 3:
                    ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color=plt.cm.YlOrRd(w / max(1.0, max_co)), alpha=0.8, lw=lw, zorder=1)
                elif v == 4:
                    if w >= np.percentile(co[np.triu_indices(n_g, 1)], 70) if n_g > 2 else True:
                        ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color="#C62828", alpha=0.75, lw=lw, zorder=1)
                else:
                    ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color="#2E7D32", alpha=0.4, lw=lw, linestyle="--", zorder=1)
        sizes = 200 + 800 * (gene_counts / max(1.0, gene_counts.max() if len(gene_counts) else 1.0))
        if v == 5:
            ax.clear()
            n_top = min(12, n_mut)
            top = order_m_freq[:n_top]
            gy = np.linspace(0.1, 0.9, n_g)
            my = np.linspace(0.1, 0.9, n_top)
            for gi, g in enumerate(genes):
                ax.scatter([0.15], [gy[gi]], s=180, color="#1565c0", zorder=3, edgecolors="white")
                ax.text(0.10, gy[gi], g, ha="right", va="center", fontsize=8, fontweight="bold")
            for mi, mj in enumerate(top):
                ax.scatter([0.85], [my[mi]], s=120, color="#E65100", zorder=3, edgecolors="white")
                ax.text(0.90, my[mi], _short_mut(labels[mj], 18), ha="left", va="center", fontsize=7)
            for gi, g in enumerate(genes):
                for mi, mj in enumerate(top):
                    if gene_burden.shape[1] > gi and mat[:, mj].sum() > 0:
                        strength = float(((gene_burden[:, gi] > 0) & (mat[:, mj] > 0)).sum())
                        if strength > 0:
                            ax.plot([0.15, 0.85], [gy[gi], my[mi]], color="#90A4AE",
                                    alpha=min(0.85, 0.15 + strength / max(1, n_iso)),
                                    lw=0.6 + 2.5 * strength / max(1, n_iso), zorder=1)
            ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
            ax.set_title("AMR relation network — gene↔mutation bipartite (style 5)", fontweight="bold")
        else:
            node_colors = {
                1: SOLID[0], 2: SOLID[3], 3: SOLID[2], 4: SOLID[4],
            }[v]
            ax.scatter(xs, ys, s=sizes, c=node_colors, alpha=0.9, edgecolors="white", linewidths=1.2, zorder=3)
            for i, g in enumerate(genes):
                ax.text(xs[i] * 1.12, ys[i] * 1.12, g, ha="center", va="center", fontsize=9, fontweight="bold")
            ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.4, 1.4)
            ax.set_aspect("equal"); ax.axis("off")
            titles = {
                1: "AMR relation network — gene co-occurrence (all edges)",
                2: "AMR relation network — strong edges only (≥30% max)",
                3: "AMR relation network — heat-coloured edges",
                4: "AMR relation network — top-quartile edges",
            }
            ax.set_title(titles[v] + f" (style {v})", fontweight="bold")
        fig.tight_layout()
        save(fig, "AMR_Relation_Network", next_num)
        next_num += 1

    _make_image_pptx(
        paths,
        os.path.join(out_dir, "Exclusive_Figures_Editable.pptx"),
        title_prefix="Exclusive Figure",
    )
    logger.info(
        "  [Exclusive_Figures] %d high-resolution PNG figures + PPTX written to %s",
        len(paths), out_dir,
    )
    return out_dir

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

def _safe_isolate_dirname(name: str, max_len: int = 80) -> str:
    s = re.sub(r'[<>:"/\\|?*\x00-\x1f]+', "_", str(name or "isolate"))
    s = re.sub(r"_+", "_", s).strip(" ._")
    return (s[:max_len] or "isolate")

def _write_fasta_record(path: str, header: str, sequence: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    seq = re.sub(r"\s+", "", str(sequence or ""))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f">{header}\n")
        if seq:
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")
        else:
            fh.write("\n")

def generate_extracted_sequences_package(
    gene_results: List[dict],
    output_root: str,
) -> Optional[str]:
    out_dir = os.path.join(output_root, EXTRACTED_SEQUENCES_FOLDER_NAME)
    os.makedirs(out_dir, exist_ok=True)
    logger.info("  [Extracted_Sequences] Writing per-isolate NT/AA under %s", out_dir)

    genes_order: List[str] = list(DEFAULT_GENES)
    by_iso: Dict[str, Dict[str, dict]] = {}
    for gr in gene_results or []:
        gene = str(gr.get("gene") or "gene")
        if gene not in genes_order:
            genes_order.append(gene)
        for d in gr.get("all_data") or []:
            iso = str(d.get("isolate") or "").strip()
            if not iso:
                continue
            by_iso.setdefault(iso, {})[gene] = {
                "nt": d.get("query_dna") or "",
                "aa": d.get("query_protein") or "",
                "status": d.get("status") or "",
                "gene_found": bool(d.get("gene_found")),
                "ref_name": d.get("ref_name") or gene,
            }

    for iso in list(by_iso.keys()):
        for gene in genes_order:
            by_iso[iso].setdefault(gene, {
                "nt": "", "aa": "", "status": "MISSING",
                "gene_found": False, "ref_name": gene,
            })

    n_iso = 0
    n_nt_ok = 0
    n_aa_ok = 0
    n_nt_files = 0
    n_aa_files = 0

    def _sort_iso(name: str):
        m = re.match(r"^(\d+)", name.strip())
        return (int(m.group(1)) if m else 10**9, name.lower())

    for iso, gene_map in sorted(by_iso.items(), key=lambda x: _sort_iso(x[0])):
        iso_dir = os.path.join(out_dir, _safe_isolate_dirname(iso))
        nt_dir = os.path.join(iso_dir, "NT")
        aa_dir = os.path.join(iso_dir, "AA")
        os.makedirs(nt_dir, exist_ok=True)
        os.makedirs(aa_dir, exist_ok=True)
        n_iso += 1

        manifest_rows = []
        for gene in genes_order:
            info = gene_map.get(gene) or {}
            nt_seq = info.get("nt") or ""
            aa_seq = info.get("aa") or ""
            found = bool(info.get("gene_found"))
            safe_gene = re.sub(r"[^A-Za-z0-9._-]+", "_", gene)
            gene_type = GENE_CONFIG.get(gene, {}).get("gene_type", "CDS")

            nt_path = os.path.join(nt_dir, f"{safe_gene}.fasta")
            if found and nt_seq:
                _write_fasta_record(nt_path, f"{iso}|{gene}|NT|len={len(nt_seq)}", nt_seq)
                n_nt_ok += 1
                nt_status = "OK"
            else:
                _write_fasta_record(
                    nt_path,
                    f"{iso}|{gene}|NT|MISSING|gene not found in assembly",
                    "",
                )
                nt_status = "MISSING"
            n_nt_files += 1

            aa_path = os.path.join(aa_dir, f"{safe_gene}.faa")
            if found and aa_seq:
                _write_fasta_record(aa_path, f"{iso}|{gene}|AA|len={len(aa_seq)}", aa_seq)
                n_aa_ok += 1
                aa_status = "OK"
            elif gene_type != "CDS":
                _write_fasta_record(
                    aa_path,
                    f"{iso}|{gene}|AA|N/A_noncoding|no protein product expected",
                    "",
                )
                aa_status = "N/A_noncoding"
            else:
                _write_fasta_record(
                    aa_path,
                    f"{iso}|{gene}|AA|MISSING|no translated AA available",
                    "",
                )
                aa_status = "MISSING"
            n_aa_files += 1

            manifest_rows.append({
                "Gene": gene,
                "Gene_type": gene_type,
                "NT_status": nt_status,
                "NT_length": len(nt_seq) if nt_seq else 0,
                "AA_status": aa_status,
                "AA_length": len(aa_seq) if aa_seq else 0,
            })

        safe_iso = _safe_isolate_dirname(iso)
        comb_nt = os.path.join(iso_dir, f"{safe_iso}_all_genes_NT.fasta")
        comb_aa = os.path.join(iso_dir, f"{safe_iso}_all_genes_AA.faa")
        with open(comb_nt, "w", encoding="utf-8") as fh_nt, \
             open(comb_aa, "w", encoding="utf-8") as fh_aa:
            for gene in genes_order:
                info = gene_map.get(gene) or {}
                nt_seq = re.sub(r"\s+", "", info.get("nt") or "")
                aa_seq = re.sub(r"\s+", "", info.get("aa") or "")
                if nt_seq:
                    fh_nt.write(f">{iso}|{gene}|NT|len={len(nt_seq)}\n")
                    for i in range(0, len(nt_seq), 70):
                        fh_nt.write(nt_seq[i:i + 70] + "\n")
                if aa_seq:
                    fh_aa.write(f">{iso}|{gene}|AA|len={len(aa_seq)}\n")
                    for i in range(0, len(aa_seq), 70):
                        fh_aa.write(aa_seq[i:i + 70] + "\n")

        man_path = os.path.join(iso_dir, "extraction_manifest.csv")
        try:
            if pd is not None:
                pd.DataFrame(manifest_rows).to_csv(man_path, index=False)
            else:
                with open(man_path, "w", newline="", encoding="utf-8") as fh:
                    w = csv.DictWriter(
                        fh,
                        fieldnames=["Gene", "Gene_type", "NT_status", "NT_length",
                                    "AA_status", "AA_length"],
                    )
                    w.writeheader()
                    w.writerows(manifest_rows)
        except Exception as e:
            logger.warning("  [Extracted_Sequences] manifest failed for %s: %s", iso, e)

    logger.info(
        "  [Extracted_Sequences] Done — %d isolate folder(s), "
        "%d NT files (%d OK), %d AA files (%d OK) → %s",
        n_iso, n_nt_files, n_nt_ok, n_aa_files, n_aa_ok, out_dir,
    )
    return out_dir

_AM_COLORS = {"substitution": "#1f77b4", "stop_gain": "#d62728", "insertion": "#2ca02c",
              "deletion": "#ff7f0e", "frameshift": "#9467bd", "other": "#7f7f7f"}

_AM_NAMES = {"substitution": "Substitution", "stop_gain": "Stop gain", "insertion": "Insertion",
             "deletion": "Deletion", "frameshift": "Frameshift", "other": "Other"}

def _am_group(t):
    t = str(t or "substitution").lower()
    if t.startswith("frameshift"):
        return "frameshift"
    if t in ("snp", "substitution", "missense"):
        return "substitution"
    if t in ("ins", "insertion"):
        return "insertion"
    if t in ("del", "deletion"):
        return "deletion"
    if t in _AM_COLORS:
        return t
    return "other"

def _am_position(name, pos):
    try:
        p = float(pos)
        if p > 0:
            return p
    except (TypeError, ValueError):
        pass
    m = re.search(r"(\d+)", str(name))
    return float(m.group(1)) if m else 0.0

def _am_build(all_data, all_mutation_data, mutation_counts, curated_set, total):
    info, ref_len = {}, 0
    for d in all_data or []:
        try:
            ref_len = max(ref_len, int(d.get("ref_length") or 0))
        except (TypeError, ValueError):
            pass
        for m in d.get("mutations") or []:
            nm = m.get("mutation")
            if nm and nm not in info:
                info[nm] = {"pos": _am_position(nm, m.get("position")), "grp": _am_group(m.get("type"))}
    muts = list(mutation_counts)
    for m in muts:
        if m not in info:
            info[m] = {"pos": _am_position(m, None), "grp": "substitution"}
    muts.sort(key=lambda k: (-mutation_counts[k], info[k]["pos"], k))
    rows = list(all_mutation_data or [])
    index = {m: i for i, m in enumerate(muts)}
    mat = np.zeros((len(rows), len(muts)), dtype=int)
    for r, d in enumerate(rows):
        for m in d.get("mutations", []):
            if m in index:
                mat[r, index[m]] = 1
    top = max([info[m]["pos"] for m in muts] + [1.0])
    length = ref_len if ref_len >= top else int(top * 1.05) + 1
    return {"muts": muts, "info": info, "iso": [d["isolate"] for d in rows], "mat": mat, "L": max(length, 1),
            "counts": dict(mutation_counts), "curated": set(curated_set or []), "total": max(int(total), 1)}

def _am_legend_handles(M, groups=None):
    from matplotlib.patches import Patch
    present = groups or [g for g in _AM_COLORS if any(M["info"][m]["grp"] == g for m in M["muts"])]
    return [Patch(facecolor=_AM_COLORS[g], label=_AM_NAMES[g]) for g in present]

def _am_landscape(M, gene, total):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    muts, mat, iso = M["muts"], M["mat"], M["iso"]
    if mat.size == 0:
        return None
    groups = list(_AM_COLORS)
    n_iso, n_mut = mat.shape
    order = sorted(range(n_iso), key=lambda r: tuple(-mat[r]))
    code = np.zeros_like(mat)
    for j, m in enumerate(muts):
        code[:, j] = mat[:, j] * (groups.index(M["info"][m]["grp"]) + 1)
    code = code[order]
    cmap = ListedColormap(["#eef1f4"] + [_AM_COLORS[g] for g in groups])
    fig = plt.figure(figsize=(max(8, 0.42 * n_mut + 4.5), max(5.5, 0.3 * n_iso + 3.2)))
    gs = fig.add_gridspec(2, 2, width_ratios=[max(8, n_mut * 0.6), 2.6],
                          height_ratios=[2.6, max(6, n_iso * 0.55)], hspace=0.04, wspace=0.04)
    ax = fig.add_subplot(gs[1, 0])
    axt = fig.add_subplot(gs[0, 0], sharex=ax)
    axr = fig.add_subplot(gs[1, 1], sharey=ax)
    ax.pcolormesh(np.arange(n_mut + 1), np.arange(n_iso + 1), code, cmap=cmap, vmin=0, vmax=len(groups),
                  edgecolors="white", linewidth=0.7)
    ax.set_ylim(n_iso, 0)
    ax.set_xticks(np.arange(n_mut) + 0.5)
    ax.set_xticklabels(muts, rotation=90, fontsize=7.5)
    ax.set_yticks(np.arange(n_iso) + 0.5)
    ax.set_yticklabels([iso[r][:28] for r in order], fontsize=7 if n_iso > 30 else 8.5)
    for lab, m in zip(ax.get_xticklabels(), muts):
        if m in M["curated"]:
            lab.set_color("#b71c1c")
            lab.set_fontweight("bold")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    pct = [M["counts"][m] / total * 100 for m in muts]
    axt.bar(np.arange(n_mut) + 0.5, pct, width=0.78, color=[_AM_COLORS[M["info"][m]["grp"]] for m in muts])
    axt.set_ylabel("% isolates", fontsize=8)
    axt.tick_params(labelbottom=False, length=0)
    burden = (code > 0).sum(axis=1)
    axr.barh(np.arange(n_iso) + 0.5, burden, height=0.78, color="#455a64")
    axr.set_xlabel("mutations", fontsize=8)
    axr.tick_params(labelleft=False, length=0)
    for a in (axt, axr):
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
    fig.legend(handles=_am_legend_handles(M), loc="lower center", ncol=6, frameon=False, fontsize=8,
               bbox_to_anchor=(0.45, -0.02))
    fig.suptitle(f"{gene}: mutation landscape (red labels = curated; n={total} isolates)", fontweight="bold", y=0.995)
    return fig

def _am_needle(M, gene, total):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    muts = M["muts"]
    xs = np.array([M["info"][m]["pos"] for m in muts], dtype=float)
    cnt = np.array([M["counts"][m] for m in muts], dtype=float)
    pct = cnt / total * 100
    L = float(M["L"])
    ymax = max(pct.max() * 1.45, 10)
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(13.5, 7), sharex=True,
                                  gridspec_kw={"height_ratios": [5, 1.3], "hspace": 0.06})
    ax.barh(0, L, left=0, height=ymax * 0.045, color="#cfd8dc", edgecolor="#90a4ae", zorder=1)
    cols = [_AM_COLORS[M["info"][m]["grp"]] for m in muts]
    ax.vlines(xs, 0, pct, colors=cols, linewidth=1.4, zorder=2)
    ax.scatter(xs, pct, s=45 + cnt * 22, c=cols, edgecolors=["black" if m in M["curated"] else "white" for m in muts],
               linewidths=1.6, zorder=3)
    placed = []
    for i in np.argsort(-cnt)[:min(len(muts), 22)]:
        y = pct[i] + ymax * 0.035
        while any(abs(xs[i] - px) < L * 0.022 and abs(y - py) < ymax * 0.15 for px, py in placed) and y < ymax:
            y += ymax * 0.15
        placed.append((xs[i], y))
        ax.text(xs[i], y, muts[i], rotation=90, ha="center", va="bottom", fontsize=7.5)
    ax.set_ylim(-ymax * 0.05, ymax)
    ax.set_ylabel("% of isolates")
    ax.set_title(f"{gene}: needle map along the sequence (n={total} isolates)", fontweight="bold")
    handles = _am_legend_handles(M) + [Line2D([0], [0], marker="o", color="w", markerfacecolor="#bdbdbd",
                                              markeredgecolor="black", markeredgewidth=1.6, markersize=9,
                                              label="Curated")]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=8, ncol=2)
    grid = np.linspace(0, L, 700)
    sigma = max(L * 0.012, 1.0)
    dens = np.zeros_like(grid)
    for x, c in zip(xs, cnt):
        dens += c * np.exp(-0.5 * ((grid - x) / sigma) ** 2)
    axd.fill_between(grid, dens, color="#78909c", alpha=0.85)
    axd.vlines(xs, 0, dens.max() * 0.18 if dens.max() > 0 else 1, colors="#263238", linewidth=0.6)
    axd.set_ylabel("density", fontsize=8)
    axd.set_yticks([])
    axd.set_xlim(0, L)
    axd.set_xlabel("Position")
    return fig

def _am_rainfall(M, gene, total):
    import matplotlib.pyplot as plt
    muts = M["muts"]
    if len(muts) < 3:
        return None
    order = sorted(muts, key=lambda m: (M["info"][m]["pos"], m))
    pos = np.array([M["info"][m]["pos"] for m in order], dtype=float)
    dist = np.diff(pos, prepend=0.0)
    y = np.log10(dist + 1)
    cnt = np.array([M["counts"][m] for m in order], dtype=float)
    fig, (ax, axh) = plt.subplots(1, 2, figsize=(12.5, 5.6), sharey=True,
                                  gridspec_kw={"width_ratios": [5, 1], "wspace": 0.03})
    cols = [_AM_COLORS[M["info"][m]["grp"]] for m in order]
    ax.scatter(pos, y, s=40 + cnt * 20, c=cols, alpha=0.85,
               edgecolors=["black" if m in M["curated"] else "white" for m in order], linewidths=1.3)
    ax.axhline(np.log10(11), color="#b71c1c", linestyle="--", linewidth=1)
    ax.text(0.995, np.log10(11), "  ≤10 residues apart", color="#b71c1c", fontsize=8, ha="right", va="bottom",
            transform=ax.get_yaxis_transform())
    for i in np.argsort(-cnt)[:min(len(order), 12)]:
        ax.annotate(order[i], (pos[i], y[i]), textcoords="offset points", xytext=(4, 6), fontsize=7.5)
    ax.set_xlim(0, M["L"])
    ax.set_xlabel("Position")
    ax.set_ylabel("log10(distance to previous mutation + 1)")
    ax.set_title(f"{gene}: rainfall plot, clustering of mutations (n={total} isolates)", fontweight="bold")
    ax.legend(handles=_am_legend_handles(M), frameon=False, fontsize=8, loc="upper left")
    axh.hist(y, bins=min(12, max(4, len(y) // 2)), orientation="horizontal", color="#90a4ae", edgecolor="white")
    axh.set_xlabel("count")
    axh.tick_params(labelleft=False)
    return fig

def _am_cluster(M, gene, total):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    try:
        from scipy.cluster import hierarchy as sch
        from scipy.spatial.distance import pdist
    except ImportError:
        return None
    mat, muts, iso = M["mat"], M["muts"], M["iso"]
    if mat.shape[0] < 2 or mat.shape[1] < 2:
        return None
    zr = sch.linkage(pdist(mat, "hamming"), "average")
    zc = sch.linkage(pdist(mat.T, "hamming"), "average")
    n_iso, n_mut = mat.shape
    fig = plt.figure(figsize=(max(8, 0.42 * n_mut + 4), max(6, 0.3 * n_iso + 3.5)))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.6, max(8, n_mut * 0.6)],
                          height_ratios=[1.8, 0.35, max(6, n_iso * 0.55)], hspace=0.03, wspace=0.03)
    axl = fig.add_subplot(gs[2, 0])
    axt = fig.add_subplot(gs[0, 1])
    axs = fig.add_subplot(gs[1, 1])
    axh = fig.add_subplot(gs[2, 1])
    dr = sch.dendrogram(zr, orientation="left", ax=axl, no_labels=True, color_threshold=0,
                        above_threshold_color="#546e7a")
    dc = sch.dendrogram(zc, orientation="top", ax=axt, no_labels=True, color_threshold=0,
                        above_threshold_color="#546e7a")
    ri, ci = dr["leaves"], dc["leaves"]
    data = mat[np.ix_(ri, ci)]
    axh.imshow(data, aspect="auto", origin="lower", cmap=ListedColormap(["#eceff1", "#c62828"]),
               interpolation="nearest")
    axh.set_xticks(np.arange(n_mut))
    axh.set_xticklabels([muts[c] for c in ci], rotation=90, fontsize=7.5)
    axh.set_yticks(np.arange(n_iso))
    axh.set_yticklabels([iso[r][:28] for r in ri], fontsize=7 if n_iso > 30 else 8.5)
    axh.yaxis.tick_right()
    axh.set_xticks(np.arange(-0.5, n_mut, 1), minor=True)
    axh.set_yticks(np.arange(-0.5, n_iso, 1), minor=True)
    axh.grid(which="minor", color="white", linewidth=0.7)
    axh.tick_params(which="both", length=0)
    for lab, c in zip(axh.get_xticklabels(), ci):
        if muts[c] in M["curated"]:
            lab.set_color("#b71c1c")
            lab.set_fontweight("bold")
    strip = np.array([[list(_AM_COLORS).index(M["info"][muts[c]]["grp"]) for c in ci]])
    axs.imshow(strip, aspect="auto", cmap=ListedColormap(list(_AM_COLORS.values())), vmin=0,
               vmax=len(_AM_COLORS) - 1, interpolation="nearest")
    axs.set_axis_off()
    for a in (axl, axt):
        a.set_axis_off()
    fig.legend(handles=_am_legend_handles(M), loc="lower left", ncol=3, frameon=False, fontsize=8,
               bbox_to_anchor=(0.02, 0.0))
    fig.suptitle(f"{gene}: clustered presence map (Hamming distance, average linkage; n={total} isolates)",
                 fontweight="bold", y=0.995)
    return fig

def _am_classes(M, gene, total):
    import matplotlib.pyplot as plt
    muts = M["muts"]
    groups = [g for g in _AM_COLORS if any(M["info"][m]["grp"] == g for m in muts)]
    sizes = [sum(1 for m in muts if M["info"][m]["grp"] == g) for g in groups]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.8), gridspec_kw={"width_ratios": [1, 2.2], "wspace": 0.18})
    wedges, _t = a1.pie(sizes, colors=[_AM_COLORS[g] for g in groups], startangle=90, counterclock=False,
                        wedgeprops={"width": 0.42, "edgecolor": "white", "linewidth": 2})
    for w, g, s in zip(wedges, groups, sizes):
        ang = np.deg2rad((w.theta1 + w.theta2) / 2)
        a1.text(0.79 * np.cos(ang), 0.79 * np.sin(ang), str(s), ha="center", va="center", color="white",
                fontweight="bold", fontsize=11)
    a1.text(0, 0, f"{len(muts)}\nunique", ha="center", va="center", fontsize=13, fontweight="bold")
    a1.legend(wedges, [f"{_AM_NAMES[g]} ({s})" for g, s in zip(groups, sizes)], loc="lower center",
              bbox_to_anchor=(0.5, -0.12), frameon=False, fontsize=8.5, ncol=2)
    a1.set_title("Mutation classes", fontweight="bold")
    nb = 14
    edges = np.linspace(0, M["L"], nb + 1)
    pos = np.array([M["info"][m]["pos"] for m in muts])
    bottom = np.zeros(nb)
    centers = (edges[:-1] + edges[1:]) / 2
    for g in groups:
        sel = np.array([M["info"][m]["grp"] == g for m in muts])
        h, _e = np.histogram(pos[sel], bins=edges)
        a2.bar(centers, h, width=(edges[1] - edges[0]) * 0.86, bottom=bottom, color=_AM_COLORS[g],
               edgecolor="white", label=_AM_NAMES[g])
        bottom += h
    ev = np.zeros(nb)
    for m in muts:
        b = min(int(np.searchsorted(edges, M["info"][m]["pos"], side="right") - 1), nb - 1)
        ev[max(b, 0)] += M["counts"][m]
    a3 = a2.twinx()
    a3.plot(centers, ev, color="#263238", marker="o", linewidth=1.8)
    a3.set_ylabel("isolate-level events", color="#263238")
    a2.set_ylabel("unique mutations")
    a2.set_xlabel("Position (binned)")
    a2.legend(frameon=False, fontsize=8, loc="upper left")
    a2.set_title("Where each class occurs along the sequence", fontweight="bold")
    a3.spines["top"].set_visible(False)
    fig.suptitle(f"{gene}: composition of mutation classes (n={total} isolates)", fontweight="bold", y=1.0)
    return fig

def _am_arcs(M, gene, total):
    import matplotlib.pyplot as plt
    muts = M["muts"]
    if len(muts) < 2:
        return None
    keep = muts[:30]
    keep = sorted(keep, key=lambda m: (M["info"][m]["pos"], m))
    idx = [M["muts"].index(m) for m in keep]
    sub = M["mat"][:, idx]
    co = sub.T @ sub
    n = len(keep)
    cmax = max(int(co[np.triu_indices(n, 1)].max()) if n > 1 else 1, 1)
    fig, ax = plt.subplots(figsize=(max(10, 0.55 * n + 3), 6.8))
    cmap = plt.get_cmap("viridis")
    t = np.linspace(0, np.pi, 60)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if co[i, j] > 0]
    for i, j in sorted(pairs, key=lambda p: co[p]):
        r = (j - i) / 2
        ax.plot((i + j) / 2 + r * np.cos(t), r * np.sin(t) * 0.7, color=cmap(co[i, j] / cmax),
                linewidth=0.8 + 4.2 * co[i, j] / cmax, alpha=0.6, solid_capstyle="round")
    cnt = np.array([M["counts"][m] for m in keep], dtype=float)
    ax.scatter(np.arange(n), np.zeros(n), s=70 + cnt * 28, c=[_AM_COLORS[M["info"][m]["grp"]] for m in keep],
               edgecolors=["black" if m in M["curated"] else "white" for m in keep], linewidths=1.6, zorder=3)
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels(keep, rotation=70, ha="right", fontsize=8)
    ax.set_yticks([])
    ax.set_xlim(-1, n)
    ax.set_ylim(-0.6, max(n * 0.37, 1.5))
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(1, cmax))
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("isolates sharing both")
    ax.legend(handles=_am_legend_handles(M), frameon=False, fontsize=8, loc="upper left")
    suffix = "" if len(M["muts"]) <= 30 else " (top 30 by frequency)"
    ax.set_title(f"{gene}: co-occurrence arcs in sequence order{suffix}", fontweight="bold")
    return fig

def _am_polar(M, gene, total):
    import matplotlib.pyplot as plt
    muts = M["muts"]
    pos = np.array([M["info"][m]["pos"] for m in muts], dtype=float)
    pct = np.array([M["counts"][m] for m in muts], dtype=float) / total * 100
    L = float(M["L"])
    base = max(pct.max(), 10) * 0.55
    top = base + max(pct.max(), 10) * 1.35
    theta = pos / L * 2 * np.pi
    fig = plt.figure(figsize=(9.5, 9.5))
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ring = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    ax.bar(ring, np.full_like(ring, base * 0.1), width=2 * np.pi / 360, bottom=base * 0.9, color="#cfd8dc",
           edgecolor="none")
    width = max(2 * np.pi / L * 2.2, 0.025)
    ax.bar(theta, pct, width=width, bottom=base, color=[_AM_COLORS[M["info"][m]["grp"]] for m in muts],
           edgecolor=["black" if m in M["curated"] else "white" for m in muts], linewidth=1.1, alpha=0.92)
    cnt = np.array([M["counts"][m] for m in muts])
    for i in np.argsort(-cnt)[:min(len(muts), 14)]:
        screen = 90 - np.degrees(theta[i])
        rot, ha = screen, "left"
        if 90 < screen % 360 < 270:
            rot, ha = screen + 180, "right"
        ax.text(theta[i], base + pct[i] + top * 0.025, muts[i], rotation=rot, rotation_mode="anchor", ha=ha,
                va="center", fontsize=8)
    ax.set_ylim(0, top)
    ax.set_yticks([])
    ticks = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(int(L * k / 8)) for k in range(8)], fontsize=8)
    ax.grid(color="#eceff1")
    ax.spines["polar"].set_visible(False)
    ax.text(0, 0, f"{gene}\nn={total}", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.legend(handles=_am_legend_handles(M), loc="lower center", bbox_to_anchor=(0.5, -0.07), ncol=3,
              frameon=False, fontsize=8.5)
    ax.set_title(f"{gene}: circular mutation map (bar height = % of isolates)", fontweight="bold", pad=18)
    return fig

def _am_combos(M, gene, total):
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    mat = M["mat"]
    profiles = {}
    for row in mat:
        key = tuple(np.nonzero(row)[0])
        if key:
            profiles[key] = profiles.get(key, 0) + 1
    if not profiles:
        return None
    ranked = sorted(profiles.items(), key=lambda kv: (-kv[1], len(kv[0])))[:20]
    used = sorted({j for k, _c in ranked for j in k})[:25]
    if not used:
        return None
    rows = {j: r for r, j in enumerate(used)}
    n_c, n_r = len(ranked), len(used)
    fig = plt.figure(figsize=(max(9, 0.6 * n_c + 5), max(5.5, 0.34 * n_r + 3.5)))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.2, max(6, n_c * 0.55)], height_ratios=[2.4, max(4, n_r * 0.5)],
                          hspace=0.04, wspace=0.32)
    axb = fig.add_subplot(gs[0, 1])
    axm = fig.add_subplot(gs[1, 1], sharex=axb)
    axl = fig.add_subplot(gs[1, 0], sharey=axm)
    axb.bar(np.arange(n_c), [c for _k, c in ranked], color="#37474f", width=0.7)
    for i, (_k, c) in enumerate(ranked):
        axb.text(i, c, str(c), ha="center", va="bottom", fontsize=8)
    axb.set_ylabel("isolates")
    axb.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    axb.tick_params(labelbottom=False, length=0)
    for i in range(n_c):
        for r in range(n_r):
            axm.scatter(i, r, s=75, color="#dfe3e6", zorder=2)
    for i, (k, _c) in enumerate(ranked):
        ys = [rows[j] for j in k if j in rows]
        if not ys:
            continue
        axm.plot([i, i], [min(ys), max(ys)], color="#37474f", linewidth=2.2, zorder=3)
        for y in ys:
            axm.scatter(i, y, s=90, color="#37474f", zorder=4)
    axm.set_ylim(n_r - 0.5, -0.5)
    axm.set_xlim(-0.6, n_c - 0.4)
    axm.set_xticks([])
    axm.set_yticks([])
    for a in (axb, axm):
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    axm.spines["left"].set_visible(False)
    axm.spines["bottom"].set_visible(False)
    axl.barh(np.arange(n_r), [M["counts"][M["muts"][j]] / total * 100 for j in used], height=0.65,
             color=[_AM_COLORS[M["info"][M["muts"][j]]["grp"]] for j in used])
    axl.set_yticks(np.arange(n_r))
    axl.set_yticklabels([M["muts"][j] for j in used], fontsize=8)
    for lab, j in zip(axl.get_yticklabels(), used):
        if M["muts"][j] in M["curated"]:
            lab.set_color("#b71c1c")
            lab.set_fontweight("bold")
    axl.invert_xaxis()
    axl.set_xlabel("% isolates", fontsize=8)
    axl.yaxis.tick_right()
    axl.yaxis.set_label_position("right")
    for sp in ("top", "left"):
        axl.spines[sp].set_visible(False)
    fig.suptitle(f"{gene}: mutation combinations per isolate (top {n_c} profiles, n={total} isolates)",
                 fontweight="bold", y=0.995)
    return fig

def _am_substitution(M, gene, total):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    pat = re.compile(r"^([A-Z\*])(\d+)([A-Z\*])$")
    cells = {}
    for m in M["muts"]:
        mm = pat.match(m)
        if mm:
            key = (mm.group(1), mm.group(3))
            cells.setdefault(key, []).append(m)
    if not cells:
        return None
    letters = {a for a, _b in cells} | {b for _a, b in cells}
    alphabet = list("ACGT") if letters <= set("ACGT") else list("ACDEFGHIKLMNPQRSTVWY*")
    alphabet += sorted(letters - set(alphabet))
    k = len(alphabet)
    grid = np.zeros((k, k))
    for (a, b), names in cells.items():
        grid[alphabet.index(a), alphabet.index(b)] = sum(M["counts"][n] for n in names)
    fig = plt.figure(figsize=(9.5 if k <= 4 else 11, 9.5 if k <= 4 else 11))
    gs = fig.add_gridspec(2, 2, width_ratios=[6, 1], height_ratios=[1, 6], hspace=0.04, wspace=0.04)
    ax = fig.add_subplot(gs[1, 0])
    axt = fig.add_subplot(gs[0, 0], sharex=ax)
    axr = fig.add_subplot(gs[1, 1], sharey=ax)
    masked = np.ma.masked_where(grid == 0, grid)
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad("#f4f6f8")
    ax.imshow(masked, cmap=cmap, aspect="equal", vmin=0, vmax=max(grid.max(), 1))
    for (a, b), names in cells.items():
        i, j = alphabet.index(a), alphabet.index(b)
        ax.text(j, i, str(int(grid[i, j])), ha="center", va="center", fontsize=8 if k > 4 else 12,
                color="white" if grid[i, j] > grid.max() * 0.6 else "#212121", fontweight="bold")
        if any(n in M["curated"] for n in names):
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="#1a237e", linewidth=2.2))
    ax.set_xticks(range(k))
    ax.set_xticklabels(alphabet)
    ax.set_yticks(range(k))
    ax.set_yticklabels(alphabet)
    ax.set_xlabel("to")
    ax.set_ylabel("from")
    ax.set_xticks(np.arange(-0.5, k, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, k, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="both", length=0)
    axt.bar(range(k), grid.sum(axis=0), color="#ef6c00", width=0.8)
    axt.tick_params(labelbottom=False, length=0)
    axr.barh(range(k), grid.sum(axis=1), color="#6a1b9a", height=0.8)
    axr.tick_params(labelleft=False, length=0)
    for a in (axt, axr):
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    fig.suptitle(f"{gene}: substitution matrix (boxed = curated; n={total} isolates)", fontweight="bold", y=0.96)
    return fig

def create_all_mutation_illustrations(all_data, all_mutation_data, mutation_counts, total_isolates, curated_set,
                                      output_dir, prefix, gene_name="gene"):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.info("  Skipping all-mutation illustrations: matplotlib not installed.")
        return []
    if not mutation_counts or total_isolates <= 0:
        return []
    os.makedirs(output_dir, exist_ok=True)
    count_all = len(all_mutation_data)
    if count_all > 60:
        all_mutation_data = sorted(all_mutation_data, key=lambda d:d["isolate"])[:60]
    M = _am_build(all_data, all_mutation_data, mutation_counts, curated_set, total_isolates)
    steps = [("1_landscape", _am_landscape), ("2_needle_map", _am_needle), ("3_rainfall", _am_rainfall),
             ("4_clustered_presence_map", _am_cluster), ("5_mutation_classes", _am_classes),
             ("6_cooccurrence_arcs", _am_arcs), ("7_circular_map", _am_polar),
             ("8_combination_profiles", _am_combos), ("9_substitution_matrix", _am_substitution)]
    written = []
    rc = {"font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False,
          "axes.grid": False, "figure.facecolor": "white"}
    for slug, fn in steps:
        fig = None
        try:
            with plt.rc_context(rc):
                fig = fn(M, gene_name, int(total_isolates))
                if fig is None:
                    continue
                if count_all > 60:
                    fig.text(.01,.003,f"Isolate displays: first 60 of {count_all} called clinical isolates; frequencies use the full cohort.", fontsize=8)
                path = os.path.join(output_dir, f"{prefix}_all_{slug}.png")
                fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
                written.append(path)
        except Exception as e:
            logger.warning(f"  all-mutation figure {slug} failed: {e}")
        finally:
            if fig is not None:
                plt.close(fig)
    return written

LIB_THRESHOLDS = {
    "16S_rDNA": (80, 95), "23S_rDNA": (80, 95), "blaTEM": (80, 95), "ereA": (80, 95), "ereB": (80, 95),
    "ermA": (80, 95), "ermB": (80, 95), "ermC": (80, 95), "ermF": (80, 95), "folP": (80, 95),
    "gyrA": (80, 95), "macAB_promoter": (80, 95), "mefA": (80, 95),
    "mtrR_promoter_mosaic_1": (99.0, 98), "mtrR_promoter_mosaic_2": (95.0, 98), "mtrR_promoter_mosaic_3": (98.5, 98),
    "mtrC": (80, 95), "mtrD_mosaic_1": (97.5, 98), "mtrD_mosaic_2": (95.0, 98), "mtrD_mosaic_3": (95.0, 98),
    "mtrR": (80, 95), "mtrR_promoter": (80, 95), "norM_promoter": (80, 95), "parC": (80, 95), "parE": (80, 95),
    "penA": (80, 95), "ponA1": (80, 95), "porB1b": (70, 95), "rpoB": (80, 95), "rpoD": (80, 95),
    "rplD": (80, 95), "rplV": (80, 95), "rpsE": (80, 95), "rpsJ": (80, 95), "tetM": (80, 95)}

RRNA = {"16S_rDNA", "23S_rDNA"}

PROMOTER_ANCHOR = {"mtrR_promoter": 251, "macAB_promoter": 226, "norM_promoter": 335}

MOSAIC_GROUPS = {"mtrD": ["mtrD_mosaic_1", "mtrD_mosaic_2", "mtrD_mosaic_3"],
                     "mtrR_promoter": ["mtrR_promoter_mosaic_1", "mtrR_promoter_mosaic_2", "mtrR_promoter_mosaic_3"]}

RRNA_COPY_MIN_PID = 98.65

FLANK = 60

_NT = _parasail.matrix_create("ACGTN", 2, -3) if _parasail is not None else None

_AA = _parasail.blosum62 if _parasail is not None else None

def _rc(s):
    return str(Seq(s).reverse_complement())

def _tr(s):
    s = s[: len(s) - len(s) % 3]
    return str(Seq(s).translate(table=11))

_SAFE = {"dir": None, "refs": {}, "n": 0}

def _safe_dir():
    if _SAFE["dir"] is None:
        _SAFE["dir"] = _tempfile.mkdtemp(prefix="ng_blast_")
        atexit.register(shutil.rmtree, _SAFE["dir"], True)
    return _SAFE["dir"]

def _run_command(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{os.path.basename(cmd[0])} exit {p.returncode}: {(p.stderr or p.stdout).strip()[:600]}")
    return p.stdout

def _blast(asm_path, refs_fasta, work=None):
    key = (os.path.abspath(asm_path), os.path.abspath(refs_fasta))
    if key in _CACHE["hits"]:
        return _CACHE["hits"][key]
    d = _safe_dir()
    if key[1] not in _SAFE["refs"]:
        q = os.path.join(d, f"refs_{len(_SAFE['refs'])}.fasta")
        shutil.copyfile(refs_fasta, q)
        _SAFE["refs"][key[1]] = q
    _SAFE["n"] += 1
    safe_asm = os.path.join(d, f"asm_{_SAFE['n']}.fna")
    shutil.copyfile(asm_path, safe_asm)
    db = os.path.join(d, f"asmdb_{_SAFE['n']}")
    _run_command(["makeblastdb", "-in", safe_asm, "-dbtype", "nucl", "-out", db])
    cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
    out = _run_command(["blastn", "-query", _SAFE["refs"][key[1]], "-db", db, "-outfmt", "6 " + cols,
                    "-evalue", "1e-10", "-dust", "no", "-max_target_seqs", "100", "-max_hsps", "100"])
    rows = []
    for line in out.strip().splitlines():
        f = line.split("\t")
        rows.append({"qseqid": f[0], "sseqid": [x for x in f[1].split("|") if x][-1], "pident": float(f[2]),
                     "length": int(f[3]), "qstart": int(f[6]), "qend": int(f[7]), "sstart": int(f[8]),
                     "send": int(f[9]), "bitscore": float(f[11]), "qlen": int(f[12]), "slen": int(f[13])})
    _CACHE["hits"][key] = rows
    return rows

def _contigs(asm_path):
    if asm_path not in _CACHE["contigs"]:
        seqs, covs = {}, {}
        for r in SeqIO.parse(asm_path, "fasta"):
            k = [p for p in r.id.split("|") if p][-1]
            seqs[k] = str(r.seq).upper()
            m = re.search(r"cov_([\d.]+)", r.description)
            covs[k] = float(m.group(1)) if m else None
        vals = [c for k, c in covs.items() if c is not None and len(seqs[k]) >= 5000]
        if len(_CACHE["contigs"]) >= 32:
            del _CACHE["contigs"][next(iter(_CACHE["contigs"]))]
        _CACHE["contigs"][asm_path] = (seqs, covs, _stats.median(vals) if vals else None)
    return _CACHE["contigs"][asm_path]

def _groups(hits, ref, gap=5000):
    h = [x for x in hits if x["qseqid"] == ref]
    for x in h:
        x["strand"] = "+" if x["send"] >= x["sstart"] else "-"
        x["s_lo"], x["s_hi"] = min(x["sstart"], x["send"]), max(x["sstart"], x["send"])
    buckets = {}
    for x in h:
        buckets.setdefault((x["sseqid"], x["strand"]), []).append(x)
    groups = []
    for (_c, _s), sub in buckets.items():
        sub.sort(key=lambda x: x["s_lo"])
        cur = []
        for x in sub:
            if cur and x["s_lo"] - max(y["s_hi"] for y in cur) > gap:
                groups.append(cur); cur = []
            cur.append(x)
        if cur:
            groups.append(cur)
    out = []
    for g in groups:
        qlen = g[0]["qlen"]; cov = set()
        for x in g:
            cov.update(range(x["qstart"], x["qend"] + 1))
        aln = sum(x["length"] for x in g)
        out.append({"contig": g[0]["sseqid"], "strand": g[0]["strand"], "s_lo": min(x["s_lo"] for x in g),
                    "s_hi": max(x["s_hi"] for x in g), "q_min": min(x["qstart"] for x in g),
                    "q_max": max(x["qend"] for x in g), "qcov": 100.0 * len(cov) / qlen,
                    "pid": sum(x["pident"] * x["length"] for x in g) / max(1, aln),
                    "bitscore": sum(x["bitscore"] for x in g), "qlen": qlen})
    return sorted(out, key=lambda x: -x["bitscore"])

def _region(seqs, g):
    seq = seqs[g["contig"]]; L = len(seq)
    if g["strand"] == "+":
        lo = g["s_lo"] - (g["q_min"] - 1) - FLANK; hi = g["s_hi"] + (g["qlen"] - g["q_max"]) + FLANK
    else:
        lo = g["s_lo"] - (g["qlen"] - g["q_max"]) - FLANK; hi = g["s_hi"] + (g["q_min"] - 1) + FLANK
    edge = lo < 1 or hi > L
    lo, hi = max(1, lo), min(L, hi)
    sub = seq[lo - 1:hi]
    return (_rc(sub) if g["strand"] == "-" else sub), edge

def _align(ref, region, edge=False):
    t = _parasail.sg_dx_trace_scan_32(ref, region, 5, 2, _NT).traceback
    a_r, a_q = t.query, t.ref
    i0 = next(i for i, c in enumerate(a_r) if c != "-")
    i1 = len(a_r) - next(i for i, c in enumerate(reversed(a_r)) if c != "-")
    a_r, a_q = a_r[i0:i1], a_q[i0:i1]
    off = 0
    if edge:
        j = 0
        while j < len(a_q) and a_q[j] == "-":
            j += 1
        k = len(a_q)
        while k > j and a_q[k - 1] == "-":
            k -= 1
        off = sum(1 for c in a_r[:j] if c != "-")
        a_r, a_q = a_r[j:k], a_q[j:k]
    return a_r, a_q, off

def _trim(a_r, a_q, off, lo, hi):
    kr, kq, pos, first = [], [], off, None
    for r, q in zip(a_r, a_q):
        if r != "-":
            pos += 1
        if lo <= pos <= hi and not (r == "-" and pos == hi):
            if first is None:
                first = pos - (1 if r != "-" else 0)
            kr.append(r); kq.append(q)
    o = first if first is not None else off
    while kq and kq[0] == "-":
        if kr[0] != "-":
            o += 1
        kr.pop(0); kq.pop(0)
    while kq and kq[-1] == "-":
        kr.pop(); kq.pop()
    return "".join(kr), "".join(kq), o

def _events(a_r, a_q, off=0, ambig="N"):
    ev, pos, i = [], off, 0
    while i < len(a_r):
        r, q = a_r[i], a_q[i]
        if r != "-" and q != "-":
            pos += 1
            if r != q and q != ambig and r != ambig:
                ev.append({"type": "SNP", "pos": pos, "ref": r, "alt": q})
            i += 1
        elif r == "-":
            j = i
            while j < len(a_r) and a_r[j] == "-":
                j += 1
            ev.append({"type": "INS", "pos": pos, "ref": "", "alt": a_q[i:j]}); i = j
        else:
            j = i; st = pos + 1
            while j < len(a_r) and a_q[j] == "-" and a_r[j] != "-":
                j += 1; pos += 1
            ev.append({"type": "DEL", "pos": st, "end": pos, "ref": a_r[i:j], "alt": ""}); i = j
    return ev

def _ins_eq(refseq, after, ins):
    target = refseq[:after] + ins + refseq[after:]
    return sorted({(a, target[a:a + len(ins)]) for a in range(len(refseq) + 1)
                   if refseq[:a] + target[a:a + len(ins)] + refseq[a:] == target})

def _del_eq(refseq, s, e):
    k = e - s + 1; target = refseq[:s - 1] + refseq[e:]
    return sorted({(a, a + k - 1) for a in range(1, len(refseq) - k + 2) if refseq[:a - 1] + refseq[a - 1 + k:] == target})

def _cds(ref_nt, a_r, a_q):
    q = a_q.replace("-", "")
    ev = _events(a_r, a_q)
    fs, shift, open_pos = [], 0, None
    for e in sorted([e for e in ev if e["type"] in ("INS", "DEL")], key=lambda x: x["pos"]):
        d = len(e["alt"]) if e["type"] == "INS" else -len(e["ref"])
        before = shift; shift += d
        if before % 3 == 0 and shift % 3 != 0:
            open_pos = e["pos"]
        elif before % 3 != 0 and shift % 3 == 0:
            if open_pos is not None and e["pos"] - open_pos > 30:
                fs.append((open_pos, e["pos"], "restored"))
            open_pos = None
    if shift % 3 != 0 and open_pos is not None:
        fs.append((open_pos, None, "unrestored"))
    ref_p = _tr(ref_nt).rstrip("*")
    q_full = _tr(q)
    start_ok = q[:3] in ("ATG", "GTG", "TTG")
    stop = q_full.find("*")
    q_p = q_full if stop < 0 else q_full[:stop]
    premature = stop >= 0 and stop < len(ref_p) and stop < len(q_full) - 1
    fs_aa = (min(f[0] for f in fs) - 1) // 3 + 1 if fs else None
    t = _parasail.nw_trace_scan_32(ref_p, q_p, 11, 1, _AA).traceback
    p_ev = _events(t.query, t.ref, ambig="X")
    if fs_aa is not None:
        p_ev = [e for e in p_ev if e["pos"] < fs_aa]
    if premature:
        p_ev = [e for e in p_ev if not (e["type"] == "DEL" and e.get("end", 0) >= len(ref_p) and e["pos"] > stop)]
    disrupted = any(f[2] == "unrestored" for f in fs) or premature or not start_ok
    return {"q_cds": q, "q_prot": q_p, "ref_prot": ref_p, "p_events": p_ev, "frameshifts": fs,
            "premature_stop_aa": (stop + 1) if premature else None, "start": q[:3], "start_ok": start_ok,
            "disrupted": disrupted}

def _tokens(e, refseq, kind, anchor=None):
    def c(p):
        if not anchor:
            return p
        x = p - anchor
        return x if x < 0 else x + 1
    if kind == "aa":
        if e["type"] == "SNP":
            t = f"{e['ref']}{e['pos']}{e['alt']}"; return t, {t}
        if e["type"] == "INS":
            eq = _ins_eq(refseq, e["pos"], e["alt"]); a, s = eq[-1]
            toks = set()
            for a2, s2 in eq:
                toks |= {f"ins{a2}{s2}", "-" * len(s2) + f"{a2}{s2}"}
            return f"ins{a}{s}", toks
        eq = _del_eq(refseq, e["pos"], e["end"]); s, en = eq[-1]
        canon = f"{refseq[s-1]}{s}-" if s == en else f"{refseq[s-1]}{s}_{refseq[en-1]}{en}del"
        return canon, [frozenset(f"{refseq[i-1]}{i}-" for i in range(a, b + 1)) for a, b in eq]
    if e["type"] == "SNP":
        t = f"{e['ref'].lower()}{c(e['pos'])}{e['alt'].lower()}"; return t, {t}
    if e["type"] == "INS":
        eq = _ins_eq(refseq, e["pos"], e["alt"]); a, s = eq[-1]
        return f"ins{c(a)}{s.lower()}", {f"ins{c(a2) if a2 else a2}{s2.lower()}" for a2, s2 in eq}
    eq = _del_eq(refseq, e["pos"], e["end"]); s, en = eq[-1]
    canon = (f"{refseq[s-1].lower()}{c(s)}del" if s == en else
             f"{refseq[s-1].lower()}{c(s)}_{refseq[en-1].lower()}{c(en)}del")
    return canon, [frozenset(f"{refseq[i-1].lower()}{c(i)}del" for i in range(a, b + 1)) for a, b in eq]

def _forms(t):
    f = {t}
    if t.endswith("-"):
        f.add(t[:-1] + "del")
    if t.endswith("del"):
        f.add(t[:-3] + "-")
    return f

def _locus_label(canon, toks, cur):
    if isinstance(toks, list):
        ok = [rep for rep in toks if rep and all(_forms(t) & cur for t in rep)]
        matched = sorted({x for rep in ok for t in rep for x in (_forms(t) & cur)})
    else:
        matched = sorted(t for t in toks if t in cur)
    return bool(matched), (matched[0] if matched else canon)

def _row(position, notation, typ, curated_flag, detail):
    return {"position": position, "mutation": notation, "type": typ, "curated": curated_flag, "detail": detail}

def _call_locus(gene, asm_path, refs_fasta, curated, work):
    if _CACHE["refs"] is None:
        _CACHE["refs"] = {r.id: str(r.seq).upper() for r in SeqIO.parse(refs_fasta, "fasta")}
    refs = _CACHE["refs"]
    cur = set(curated.get(gene, []))
    hits = _blast(asm_path, refs_fasta, work)
    seqs, covs, med = _contigs(asm_path)
    res = {"gene_found": False, "status": "absent", "pident": None, "qcov": None, "ref_length": None,
           "query_length": None, "ref_name": gene, "query_dna": "", "query_protein": "", "ref_protein": "",
           "mutations": [], "is_disrupted": False, "mosaic_typing": {}}
    rows = res["mutations"]

    if gene in MOSAIC_GROUPS:
        for m in MOSAIC_GROUPS[gene]:
            gs = _groups(hits, m)
            pt, ct = LIB_THRESHOLDS[m]
            b = gs[0] if gs else None
            res["mosaic_typing"][m] = {"pid": round(b["pid"], 3) if b else None, "qcov": round(b["qcov"], 2) if b else None,
                                       "present": bool(b and b["pid"] >= pt and b["qcov"] >= ct)}
    if gene == "mtrD":
        best, best_id = None, -1
        for m in MOSAIC_GROUPS["mtrD"]:
            gs = _groups(hits, m)
            if not gs:
                continue
            region, edge = _region(seqs, gs[0])
            a_r, a_q, _ = _align(refs[m], region, edge)
            ident = 100.0 * sum(1 for x, y in zip(a_r, a_q) if x == y) / len(a_r)
            if ident > best_id:
                best, best_id = (m, gs[0], a_r, a_q), ident
        if not best:
            return res
        m, g, a_r, a_q = best
        r = _cds(refs[m], a_r, a_q)
        res.update({"gene_found": True, "status": "disrupted" if r["disrupted"] else "intact", "pident": round(g["pid"], 2),
                    "qcov": round(g["qcov"], 2), "ref_length": len(r["ref_prot"]), "query_length": len(r["q_prot"]),
                    "ref_name": m, "query_dna": r["q_cds"], "query_protein": r["q_prot"], "ref_protein": r["ref_prot"]})
        for e in r["p_events"]:
            canon, _t = _tokens(e, r["ref_prot"], "aa")
            typ = {"SNP": "substitution", "INS": "insertion", "DEL": "deletion"}[e["type"]]
            rows.append(_row(e["pos"], canon, typ, False, f"vs {m} (closest mosaic reference; reference divergence)"))
        return res

    pt, ct = LIB_THRESHOLDS[gene]
    groups = _groups(hits, gene)
    if gene in RRNA:
        groups = [g for g in groups if g["pid"] >= RRNA_COPY_MIN_PID and g["qcov"] >= 10]
        cov_union = set()
        for g in groups:
            cov_union.update(range(g["q_min"], g["q_max"] + 1))
        present = bool(groups) and 100.0 * len(cov_union) / groups[0]["qlen"] >= ct
    else:
        present = bool(groups) and groups[0]["pid"] >= pt and groups[0]["qcov"] >= ct
    if not present:
        if groups:
            b = groups[0]
            res.update({"pident": round(b["pid"], 2), "qcov": round(b["qcov"], 2), "status": "below library threshold"})
            if b["pid"] >= pt and b["qcov"] < ct and _region(seqs, b)[1]:
                res["status"] = "partial (contig edge) - presence probable, below library coverage"
        return res
    b = groups[0]
    res.update({"gene_found": True, "pident": round(b["pid"], 2), "qcov": round(b["qcov"], 2), "ref_length": len(refs[gene])})
    anchor = PROMOTER_ANCHOR.get(gene)
    seen = set()
    for ci, g in enumerate(groups if gene in RRNA else [b], 1):
        region, edge = _region(seqs, g)
        a_r, a_q, off = _align(refs[gene], region, edge)
        if g["qcov"] < 99.0:
            a_r, a_q, off = _trim(a_r, a_q, off, g["q_min"], g["q_max"])
        covered = sum(1 for x in a_r if x != "-")
        if gene in CODING:
            if not (off == 0 and covered == len(refs[gene])):
                res.update({"status": "partial (contig edge) - variants not called"}); return res
            r = _cds(refs[gene], a_r, a_q)
            res.update({"status": "disrupted" if r["disrupted"] else "intact", "is_disrupted": r["disrupted"],
                        "query_dna": r["q_cds"], "query_protein": r["q_prot"], "ref_protein": r["ref_prot"],
                        "ref_length": len(r["ref_prot"]), "query_length": len(r["q_prot"])})
            for e in r["p_events"]:
                canon, toks = _tokens(e, r["ref_prot"], "aa")
                is_cur, disp = _locus_label(canon, toks, cur)
                typ = {"SNP": "substitution", "INS": "insertion", "DEL": "deletion"}[e["type"]]
                if e["type"] == "SNP" and e["alt"] == "*":
                    typ = "stop_gain"
                det = f"canonical {canon}"
                if gene == "blaTEM":
                    p = e["pos"]; amb = p + 2 if p <= 236 else (p + 3 if p <= 249 else p + 4)
                    det += f"; Ambler {amb}"
                rows.append(_row(e["pos"], disp, typ, is_cur, det))
            if r["disrupted"]:
                why = []
                if r["frameshifts"]:
                    why.append(f"frameshift(s) {r['frameshifts']}")
                if r["premature_stop_aa"]:
                    why.append(f"premature stop at aa {r['premature_stop_aa']}")
                if not r["start_ok"]:
                    why.append(f"start codon {r['start']}")
                rows.append(_row(0, "disrupted", "disruption", "disrupted" in cur, "; ".join(why)))
        else:
            res.update({"status": "present", "query_dna": a_q.replace("-", ""), "query_length": len(a_q.replace("-", ""))})
            for e in _events(a_r, a_q, off):
                canon, toks = _tokens(e, refs[gene], "nt", anchor)
                is_cur, disp = _locus_label(canon, toks, cur)
                if disp in seen:
                    continue
                seen.add(disp)
                typ = {"SNP": "SNP", "INS": "insertion", "DEL": "deletion"}[e["type"]]
                det = f"ref index {e['pos']}"
                if anchor:
                    det += f"; PAARSNP-formula coord {e['pos'] - (len(refs[gene]) - 49)}"
                if gene in RRNA:
                    cc = covs.get(g["contig"])
                    det += f"; copy {ci}; est. collapsed copies {round(cc / med, 2) if cc and med else 'n/a'}"
                rows.append(_row(e["pos"], disp, typ, is_cur, det))
    return res

def process_gene(
    gene_name: str,
    refs_fasta: str,
    curated: dict,
    merged: dict,
    output_root: str,
    thresholds: List[int],
) -> Optional[dict]:
    logger.info("\n" + "=" * 70)
    logger.info(f"GENE: {gene_name}")
    PREFIX = gene_name
    logger.info("=" * 70)
    gene_out = os.path.join(output_root, "Mutation_Tables")
    os.makedirs(gene_out, exist_ok=True)
    work = os.path.join(output_root, "_blastdb")


    def _k(key):
        disp = merged[key].get("display") or key
        m = re.match(r"^(\d+)", disp.strip())
        return (int(m.group(1)) if m else 10 ** 9, disp)

    curated_set = set(curated.get(gene_name, []))
    all_data, all_mutation_data = [], []
    processed = skipped = 0
    for key in sorted(merged, key=_k):
        entry = merged[key]
        disp = entry.get("display") or key
        asm = entry.get("plane")
        if not asm:
            logger.warning(f"  {disp}: no assembly (plane) file - calls only from assemblies; skipped")
            all_data.append({"isolate": disp, "gene_found": False, "status": "no assembly", "mutations": []})
            skipped += 1
            continue
        try:
            r = _call_locus(gene_name, asm, refs_fasta, curated, work)
        except Exception as exc:
            logger.error("%s %s: %s", gene_name, disp, exc)
            r = {"gene_found": False, "status": f"assembly failed: {exc}", "mutations": []}
        r["isolate"] = disp
        all_data.append(r)
        if r["gene_found"]:
            processed += 1
            all_mutation_data.append({"isolate": disp, "mutations": [m["mutation"] for m in r["mutations"]]})
            curated_set |= {m["mutation"] for m in r["mutations"] if m["curated"]}
        else:
            skipped += 1
            logger.info(f"  {disp}: {r['status']}")

    csv_out = os.path.join(gene_out, f"{PREFIX}_protein_summary.csv")
    with open(csv_out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Isolate", "Pident", "Qcov", "Ref_Length", "Query_Length", "Status", "Mutation", "Type",
                    "Position", "Curated", "Detail", "Mosaic_typing"])
        for s in all_data:
            mos = "; ".join(f"{k}:{'present' if v['present'] else 'absent'}({v['pid']})" for k, v in s.get("mosaic_typing", {}).items())
            if not s.get("mutations"):
                w.writerow([s["isolate"], s.get("pident"), s.get("qcov"), s.get("ref_length"), s.get("query_length"),
                            s.get("status"), "" if s.get("gene_found") else "GENE NOT FOUND", "", "", "",
                            "No mutations" if s.get("gene_found") else s.get("status"), mos])
            for m in s.get("mutations", []):
                w.writerow([s["isolate"], s.get("pident"), s.get("qcov"), s.get("ref_length"), s.get("query_length"),
                            s.get("status"), m["mutation"], m["type"], m["position"],
                            "Curated" if m["curated"] else "Non-Curated", m["detail"], mos])
    all_unique = sorted({m for d in all_mutation_data for m in d["mutations"]})
    matrix_csv = os.path.join(gene_out, f"{PREFIX}_mutation_matrix.csv")
    with open(matrix_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Isolate"] + all_unique)
        for d in all_mutation_data:
            s = set(d["mutations"])
            w.writerow([d["isolate"]] + ["YES" if m in s else "NO" for m in all_unique])
    clinical_data = [d for d in all_data if d["isolate"] != QC]
    clinical_calls = [d for d in all_mutation_data if d["isolate"] != QC]
    total_isolates = len(clinical_data)
    mutation_counts: Dict[str, int] = {}
    for d in clinical_calls:
        for m in d["mutations"]:
            mutation_counts[m] = mutation_counts.get(m, 0) + 1
    if mutation_counts:
        overall_dir = os.path.join(output_root, "Mutation_Figures", "overall_mutation_charts")
        os.makedirs(overall_dir, exist_ok=True)
        try:
            create_overall_mutation_charts(mutation_counts=mutation_counts, total_isolates=total_isolates,
                                           curated_set=curated_set, output_dir=overall_dir, prefix=PREFIX,
                                           gene_name=gene_name)
        except Exception as e:
            logger.error(f"  overall charts failed: {e}")
        all_mut_dir = os.path.join(output_root, "Mutation_Figures", "all_mutation_charts")
        try:
            create_all_mutation_illustrations(all_data=clinical_data, all_mutation_data=clinical_calls,
                                              mutation_counts=mutation_counts, total_isolates=total_isolates,
                                              curated_set=curated_set, output_dir=all_mut_dir, prefix=PREFIX,
                                              gene_name=gene_name)
        except Exception as e:
            logger.error(f"  all-mutation illustrations failed: {e}")
    return {"gene": gene_name, "processed": processed, "skipped": skipped, "mutation_counts": mutation_counts,
            "total_isolates": total_isolates, "csv": csv_out, "matrix_csv": matrix_csv, "excel": str(Path(output_root) / "Mutation_Reports" / "By_Gene" / (safe_name(gene_name) + ".xlsx")),
            "all_mutation_data": all_mutation_data, "all_data": all_data, "curated_set": curated_set,
            "n_found": processed, "n_collection": total_isolates}

def style_workbook(wb):
    from openpyxl.styles import Alignment, Font, PatternFill
    for sheet in wb:
        sheet.freeze_panes="B2"
        sheet.auto_filter.ref=sheet.dimensions
        sheet.row_dimensions[1].height=38
        for cell in sheet[1]:
            cell.font=Font(color="FFFFFF",bold=True)
            cell.fill=PatternFill("solid",fgColor="234D62")
            cell.alignment=Alignment(wrap_text=True,vertical="center")
            sheet.column_dimensions[cell.column_letter].width=min(52,max(17,len(str(cell.value))*.75))

def category(value):
    return {"s":"S","sensitive":"S","susceptible":"S","i":"I","intermediate":"I","r":"R","resistant":"R"}.get(str(value).strip().lower())

def method_profiles(book, root, field, audit):
    """Paginated clinical MIC profiles, with bounds shown as entered."""
    keys=sorted(k for k in book.ids if k!=QC and any(book.values.get((k,d),{}).get(field) is not None for d in ANTIBIOTIC_ORDER))
    if not keys:return
    target=root/("EST_rounded" if field==ROUNDED else "AD")
    target.mkdir()
    for offset in range(0,len(keys),32):
        selected=keys[offset:offset+32]
        vals=np.full((len(selected),len(ANTIBIOTIC_ORDER)),np.nan)
        labels={}
        for i,key in enumerate(selected):
            for j,drug in enumerate(ANTIBIOTIC_ORDER):
                value=book.values.get((key,drug),{}).get(field)
                op,amount=measurement(value)
                if amount:
                    labels[i,j]=str(value)
                    if op in ("","="):vals[i,j]=math.log2(amount)
        fig,ax=plt.subplots(figsize=(11,max(3,len(selected)*.34+2)),layout="constrained")
        cmap=plt.get_cmap("Blues").copy();cmap.set_bad("#f0f0f0")
        lo,hi=(np.nanmin(vals),np.nanmax(vals)) if np.isfinite(vals).any() else (0,1)
        im=ax.imshow(np.ma.masked_invalid(vals),aspect="auto",cmap=cmap,vmin=lo,vmax=max(hi,lo+.01))
        for (i,j),v in labels.items():
            ax.text(j,i,v,ha="center",va="center",fontsize=8,color="white" if vals[i,j]>lo+.65*(hi-lo) else "#152a36")
        ax.set_xticks(range(len(ANTIBIOTIC_ORDER)),[ANTIBIOTIC_SHORT[d] for d in ANTIBIOTIC_ORDER])
        ax.set_yticks(range(len(selected)),selected,fontsize=8)
        ax.set_title(f"{field} (µg/mL) • isolates {offset+1}–{offset+len(selected)} of {len(keys)}")
        fig.colorbar(im,label="log₂ exact MIC",shrink=.75)
        fig.savefig(target/f"Clinical_MIC_{offset+1}.png",dpi=150)
        plt.close(fig)

def agreement_tables(wb,book):
    """Recorded categories with explicit method pairs and denominators."""
    from itertools import combinations
    from scipy.stats import binomtest, spearmanr
    ws=wb.create_sheet("Method pairs")
    ws.append(["Antibiotic","Method 1","Method 2","Paired S/I/R","Same S/I/R","Agreement (%)",
               "Paired R or NS binary","Only method 1 R or NS","Only method 2 R or NS","Exact McNemar p"])
    ref=wb.create_sheet("AD comparator errors")
    ref.append(["Antibiotic","S/I/R pairs","Major errors (EST R, AD S)","AD S denominator","Major errors (%)",
                "Very major errors (EST S, AD R)","AD R denominator","Very major errors (%)","Minor errors (I versus S or R)","Minor errors (%)"])
    for drug in ANTIBIOTIC_ORDER:
        records=[r for (k,d),r in book.values.items() if d==drug and k!=QC]
        for left,right in combinations(CATEGORIES,2):
            binary=[];pairs=[]
            for row in records:
                a,b=row.get(left),row.get(right)
                x,y=category(a),category(b)
                if x and y:pairs.append((x,y))
                def resistant(value):
                    k=re.sub(r"[^a-z]","",str(value).lower())
                    return True if k in {"ns","nonsusceptible"} else None if not category(value) else category(value)=="R"
                x,y=resistant(a),resistant(b)
                if x is not None and y is not None:binary.append((x,y))
            n=len(pairs);matches=sum(x==y for x,y in pairs)
            a=sum(x and not y for x,y in binary);b=sum(y and not x for x,y in binary)
            ws.append([drug,left,right,n,matches,100*matches/n if n else None,len(binary),a,b,
                       binomtest(a,a+b,.5).pvalue if a+b else 1 if binary else None])
            if left==EST_CAT and right==AD_CAT:
                major=sum(x=="R" and y=="S" for x,y in pairs);vmajor=sum(x=="S" and y=="R" for x,y in pairs)
                sd=sum(y=="S" for _,y in pairs);rd=sum(y=="R" for _,y in pairs)
                minor=sum((x=="I")!=(y=="I") for x,y in pairs)
                ref.append([drug,n,major,sd,100*major/sd if sd else None,vmajor,rd,100*vmajor/rd if rd else None,minor,100*minor/n if n else None])

def clinical_gene_profiles(results,book,root,audit):
    # Retrospective field-specific plots require a called locus and an exact
    # measured MIC for the same clinical isolate. Failed calls remain missing.
    for method in (ROUNDED,AD_MIC):
        folder=root/("EST_rounded_genotype_phenotype" if method==ROUNDED else "AD_genotype_phenotype")
        found=False
        for gr in results:
            rows=[r for r in gr["all_data"] if r.get("gene_found") and r["isolate"]!=QC]
            table=[]
            for drug in ANTIBIOTIC_ORDER:
                groups=[[],[]]
                for r in rows:
                    op,val=measurement(book.values.get((r["isolate"],drug),{}).get(method))
                    if val and op in ("","="):groups[bool(r["mutations"])].append(math.log2(val))
                table.append((drug,groups))
            if not any(a or b for _,(a,b) in table):continue
            found=True;folder.mkdir(exist_ok=True)
            fig,axes=plt.subplots(2,4,figsize=(14,7),layout="constrained")
            for ax,(drug,groups) in zip(axes.flat,table):
                for j,vals in enumerate(groups):
                    if vals:ax.boxplot(vals,positions=[j],widths=.5)
                ax.set_xticks([0,1],[f"No call of a variant\nn={len(groups[0])}",f"Variant observed\nn={len(groups[1])}"],fontsize=8)
                ax.set(title=drug,ylabel="log₂ MIC (µg/mL)")
            fig.suptitle(f"{gr['gene']} • {method} • clinical isolates with a called locus")
            fig.savefig(folder/(safe_name(gr["gene"])+".png"),dpi=150)
            plt.close(fig)
        if not found:audit.add(method+" genomic profiles","skipped","no matched clinical locus calls and exact MICs")

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

def build_mic_illustration_pptx(gallery_dir: str, manifest: List[str], output_root: str) -> Optional[str]:
    try:
        from pptx import Presentation
        from pptx.util import Inches, Pt
        from pptx.dml.color import RGBColor
        from pptx.enum.text import PP_ALIGN
    except ImportError as e:
        logger.error("  [PPTX] python-pptx not installed - skip MIC illustration gallery PPTX")
        logger.info(f"         Detail: {e}")
        return None
    if not manifest:
        logger.info("  [PPTX] No illustration PNGs to embed")
        return None

    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    blank = prs.slide_layouts[6]

    slide = prs.slides.add_slide(blank)
    box = slide.shapes.add_textbox(Inches(0.8), Inches(2.2), Inches(11.7), Inches(2.5))
    tf = box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = "MIC and mutation reports"
    p.font.size = Pt(28); p.font.bold = True
    p.font.color.rgb = RGBColor(0x15, 0x65, 0xC0)
    p.alignment = PP_ALIGN.CENTER
    p2 = tf.add_paragraph()
    p2.text = f"{len(manifest)} figures | {Path(gallery_dir).name}"
    p2.font.size = Pt(14); p2.font.color.rgb = RGBColor(0x42, 0x42, 0x42)
    p2.alignment = PP_ALIGN.CENTER

    for idx, abs_p in enumerate(manifest, 1):
        slide = prs.slides.add_slide(blank)
        cap = slide.shapes.add_textbox(Inches(0.3), Inches(0.15), Inches(12.7), Inches(0.45))
        cp = cap.text_frame.paragraphs[0]
        cp.text = f"[{idx}/{len(manifest)}]  {os.path.basename(abs_p)}"
        cp.font.size = Pt(10); cp.font.color.rgb = RGBColor(0x37, 0x47, 0x4F)
        try:
            slide.shapes.add_picture(abs_p, Inches(0.3), Inches(0.7), width=Inches(12.7))
        except Exception as e:
            note = slide.shapes.add_textbox(Inches(1), Inches(3), Inches(11), Inches(1))
            note.text_frame.paragraphs[0].text = f"Could not embed: {abs_p}\n{e}"

    out_path = os.path.join(gallery_dir, "NG_AMR_MIC_Comparison_Illustration.pptx")
    try:
        prs.save(out_path)
        logger.info(f"  [PPTX] MIC illustration gallery -> {out_path}")
        return out_path
    except Exception as e:
        logger.error(f"  [PPTX] Failed to save illustration gallery: {e}")
        return None

def build_run_summary_pptx(output_root: str, summary: dict) -> Optional[str]:
    try:
        from pptx import Presentation
        from pptx.util import Inches, Pt
        from pptx.dml.color import RGBColor
        from pptx.enum.text import PP_ALIGN
        from pptx.chart.data import CategoryChartData
        from pptx.enum.chart import XL_CHART_TYPE
    except ImportError:
        logger.info("  [Run-Summary] python-pptx not available — skip summary PPTX")
        return None

    BLUE = RGBColor(0x15, 0x65, 0xC0)
    DARK = RGBColor(0x21, 0x21, 0x21)
    GREY = RGBColor(0x61, 0x61, 0x61)
    GREEN = RGBColor(0x2E, 0x7D, 0x32)
    ORANGE = RGBColor(0xEF, 0x6C, 0x00)
    RED = RGBColor(0xC6, 0x28, 0x28)
    WHITE = RGBColor(0xFF, 0xFF, 0xFF)
    HEADER = BLUE
    ALT = RGBColor(0xE3, 0xF2, 0xFD)

    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    blank = prs.slide_layouts[6]

    def _p(para, text, size=14, bold=False, color=DARK):
        para.text = str(text)
        para.font.size = Pt(size)
        para.font.bold = bold
        para.font.color.rgb = color
        para.font.name = "Calibri"

    def _title(slide, text, size=20):
        box = slide.shapes.add_textbox(Inches(0.4), Inches(0.2), Inches(12.5), Inches(0.55))
        _p(box.text_frame.paragraphs[0], text, size=size, bold=True, color=BLUE)

    def _cell(cell, text, bold=False, size=11, fill=None, font=None):
        cell.text = str(text)
        for para in cell.text_frame.paragraphs:
            para.font.size = Pt(size)
            para.font.bold = bold
            para.font.name = "Calibri"
            para.alignment = PP_ALIGN.CENTER
            if font is not None:
                para.font.color.rgb = font
        if fill is not None:
            try:
                cell.fill.solid()
                cell.fill.fore_color.rgb = fill
            except Exception as _exc:
                logger.debug("Ignored non-critical exception: %s", _exc)
                pass

    slide = prs.slides.add_slide(blank)
    _title(slide, "NG AMR Master Mutation Pipeline — Run Summary", size=24)
    box = slide.shapes.add_textbox(Inches(0.6), Inches(1.2), Inches(12), Inches(5.5))
    tf = box.text_frame
    tf.word_wrap = True
    lines = [
        f"Status: {summary.get('status', 'COMPLETE')}",
        f"Started:  {summary.get('start_time', '—')}",
        f"Finished: {summary.get('end_time', '—')}",
        f"Elapsed:  {summary.get('elapsed_human', '—')}  ({summary.get('elapsed_sec', 0):.1f} seconds)",
        "",
        f"Output folder: {summary.get('output_root', '—')}",
        "Script: run_amr.py",
    ]
    for i, line in enumerate(lines):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        _p(para, line, size=16, bold=(i == 0), color=DARK)

    slide = prs.slides.add_slide(blank)
    _title(slide, "Timing and processing counts")
    rows_data = [
        ["Metric", "Value"],
        ["Wall-clock time", summary.get("elapsed_human", "—")],
        ["Elapsed (seconds)", f"{summary.get('elapsed_sec', 0):.2f}"],
        ["Genes requested", str(summary.get("genes_requested", 0))],
        ["Genes processed OK", str(summary.get("genes_ok", 0))],
        ["Genes failed / skipped", str(summary.get("genes_failed", 0))],
        ["Total distinct mutations (sum)", str(summary.get("total_distinct_muts", 0))],
        ["Total mutation events (sum)", str(summary.get("total_events", 0))],
        ["HTML alignment export", "Yes" if summary.get("write_html") else "No"],
        ["MIC workbook used", "Yes" if summary.get("mic_used") else "No"],
        ["Thresholds", summary.get("thresholds_str", "—")],
    ]
    tbl = slide.shapes.add_table(len(rows_data), 2, Inches(1.5), Inches(1.0),
                                 Inches(10), Inches(min(5.8, 0.38 * len(rows_data)))).table
    for ri, row in enumerate(rows_data):
        for ci, val in enumerate(row):
            fill = HEADER if ri == 0 else (ALT if ri % 2 == 0 else WHITE)
            font = WHITE if ri == 0 else DARK
            _cell(tbl.cell(ri, ci), val, bold=(ri == 0 or ci == 0), size=12, fill=fill, font=font)

    slide = prs.slides.add_slide(blank)
    _title(slide, "Pipeline stages executed")
    stages = summary.get("stages") or []
    headers = ["#", "Stage", "Status", "Note"]
    n = len(stages) + 1
    tbl = slide.shapes.add_table(n, 4, Inches(0.4), Inches(0.95), Inches(12.5),
                                 Inches(min(6.0, 0.32 * n + 0.2))).table
    for c, h in enumerate(headers):
        _cell(tbl.cell(0, c), h, bold=True, size=11, fill=HEADER, font=WHITE)
    for i, st in enumerate(stages):
        status = st.get("status", "?")
        if status == "OK":
            sc = GREEN
        elif status == "FAILED":
            sc = RED
        elif status == "SKIPPED":
            sc = ORANGE
        else:
            sc = GREY
        fill = ALT if i % 2 == 0 else WHITE
        vals = [i + 1, st.get("name", ""), status, st.get("note", "")]
        for c, v in enumerate(vals):
            _cell(tbl.cell(i + 1, c), v, bold=(c == 2), size=10, fill=fill,
                  font=(sc if c == 2 else DARK))

    gene_rows = summary.get("gene_rows") or []
    if gene_rows:
        slide = prs.slides.add_slide(blank)
        _title(slide, "Per-gene processing summary")
        headers = ["Gene", "Processed", "Missing", "Distinct muts"]
        show = gene_rows[:28]
        tbl = slide.shapes.add_table(len(show) + 1, 4, Inches(0.5), Inches(0.95),
                                     Inches(12.3), Inches(min(6.0, 0.22 * len(show) + 0.3))).table
        for c, h in enumerate(headers):
            _cell(tbl.cell(0, c), h, bold=True, size=10, fill=HEADER, font=WHITE)
        for i, gr in enumerate(show):
            fill = ALT if i % 2 == 0 else WHITE
            for c, v in enumerate([gr.get("gene"), gr.get("processed"), gr.get("missing"), gr.get("muts")]):
                _cell(tbl.cell(i + 1, c), v, bold=(c == 0), size=9, fill=fill, font=DARK)

    if stages:
        slide = prs.slides.add_slide(blank)
        _title(slide, "Stage outcome counts")
        from collections import Counter
        c = Counter(st.get("status", "?") for st in stages)
        cats = list(c.keys())
        vals = [c[k] for k in cats]
        try:
            cd = CategoryChartData()
            cd.categories = cats
            cd.add_series("Stages", vals)
            slide.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(2.0), Inches(1.1),
                                   Inches(9.5), Inches(5.5), cd)
        except Exception as _exc:
            logger.debug("Ignored non-critical exception: %s", _exc)
            pass

    slide = prs.slides.add_slide(blank)
    _title(slide, "Notes")
    box = slide.shapes.add_textbox(Inches(0.6), Inches(1.1), Inches(12), Inches(5.5))
    tf = box.text_frame
    tf.word_wrap = True
    note_lines = [
        "This summary is generated automatically at the end of the pipeline run.",
        "Stage status OK = completed without exception; FAILED = exception caught; SKIPPED = not run.",
        "Elapsed time is wall-clock (includes all I/O, figure generation, and PPTX building).",
        "Major output packages (when OK): master matrix, analysis figures, thesis package,",
        "extended galleries + editable PPTX, Final_for_Thesis (+ editable thesis deck).",
        f"High-level stages recorded: {len(stages)}",
        f"Genes with results: {summary.get('genes_ok', 0)} / {summary.get('genes_requested', 0)}",
    ]
    for i, line in enumerate(note_lines):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        _p(para, line, size=14, bold=False, color=DARK)

    out_path = os.path.join(output_root, "NG_AMR_Pipeline_Run_Summary.pptx")
    prs.save(out_path)
    logger.info(f"\n  [Run-Summary] → {out_path}  ({len(prs.slides)} slides)")
    return out_path
