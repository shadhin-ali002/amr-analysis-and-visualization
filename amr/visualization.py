"""Visualization routines."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .inputs import (
    logger,
    re,
)
from .mutation_reports import (
    Dict,
    List,
    Optional,
    os,
)

MATPLOTLIB_OK = True

FIG_TOP_MUT_HEAT = 40

FIG_TOP_MUT_BAR = 25

FIG_TOP_MUT_ANNOT = 12

def _fig_short_label(lab, n: int = 22) -> str:
    s = str(lab or "")
    return s if len(s) <= n else (s[: max(1, n - 1)] + "…")

def _fig_top_indices(counts, k: int):
    import numpy as _np
    arr = _np.asarray(counts, dtype=float)
    if arr.size == 0:
        return _np.array([], dtype=int)
    k = int(min(max(1, k), arr.size))
    return _np.argsort(-arr)[:k]

def _mutation_sort_key(notation: str):
    m = re.search(r"(\d+)", notation or "")
    pos = int(m.group(1)) if m else 0
    return (pos, notation or "")

def create_master_analysis_figures(gene_results: List[dict], output_root: str) -> Optional[str]:
    if not MATPLOTLIB_OK:
        logger.info("  matplotlib missing – skip master analysis figures")
        return None
    if not gene_results:
        return None

    import numpy as np
    from matplotlib.patches import Patch, FancyBboxPatch, Rectangle
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.sankey import Sankey

    try:
        from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
        from scipy.spatial.distance import pdist
        SCIPY_OK = True
    except ImportError:
        SCIPY_OK = False

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
    if not isolates:
        logger.info("  No isolates – skip master analysis figures")
        return None

    genes = [gr["gene"] for gr in gene_results]
    n_iso = len(isolates)
    iso_idx = {iso: i for i, iso in enumerate(isolates)}

    gene_present = np.zeros((n_iso, len(genes)), dtype=float)
    mut_labels: List[str] = []
    mut_gene: List[str] = []
    mut_curated: List[bool] = []
    mut_counts: List[int] = []
    gene_mut_isos: Dict[str, Dict[str, set]] = {}
    curated_by_gene: Dict[str, set] = {}

    for gi, gr in enumerate(gene_results):
        gene = gr["gene"]
        curated_set = set(gr.get("curated_set") or set())
        curated_by_gene[gene] = curated_set
        data_by_iso = {d["isolate"]: d for d in (gr.get("all_data") or [])}
        mutdata_by_iso = {
            d["isolate"]: set(d.get("mutations") or [])
            for d in (gr.get("all_mutation_data") or [])
        }
        gene_mut_isos[gene] = {}
        for iso in isolates:
            ad = data_by_iso.get(iso)
            found = bool(ad and ad.get("gene_found", ad.get("status") == "OK"))
            gene_present[iso_idx[iso], gi] = 1.0 if found else 0.0
            if found:
                for m in mutdata_by_iso.get(iso, set()):
                    gene_mut_isos[gene].setdefault(m, set()).add(iso)

    for gene in genes:
        for m in sorted(gene_mut_isos.get(gene, {}), key=_mutation_sort_key):
            mut_labels.append(f"{gene}_{m}")
            mut_gene.append(gene)
            mut_curated.append(m in curated_by_gene.get(gene, set()))
            mut_counts.append(len(gene_mut_isos[gene][m]))

    n_mut = len(mut_labels)
    mut_matrix = np.zeros((n_iso, n_mut), dtype=float)
    for j, (gene, label) in enumerate(zip(mut_gene, mut_labels)):
        m = label.split("_", 1)[1] if "_" in label else label
        m = label[len(gene) + 1:]
        for iso in gene_mut_isos[gene].get(m, set()):
            mut_matrix[iso_idx[iso], j] = 1.0

    def _short(iso: str) -> str:
        m = re.search(r"IsolateID[_-]?(.+?)(?:_GC[AF]\d|$)", iso, re.I)
        if m:
            return m.group(1).strip("_")[:18]
        m = re.match(r"^(\d+)", iso)
        return m.group(1) if m else iso[:18]

    short_isos = [_short(i) for i in isolates]
    out_dir = os.path.join(output_root, "master_analysis_figures")
    os.makedirs(out_dir, exist_ok=True)

    def _save(fig, name):
        path = os.path.join(out_dir, name)
        fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        return path

    curated_color = "#2e7d32"
    noncurated_color = "#e65100"
    present_color = "#1565c0"
    absent_color = "#eceff1"

    try:
        if n_mut > 0:
            mut_order = _fig_top_indices(mut_counts, FIG_TOP_MUT_HEAT)
            n_show = len(mut_order)
            burden = mut_matrix.sum(axis=1)
            iso_order = np.argsort(-burden)
            mat = mut_matrix[iso_order][:, mut_order]
            labels_o = [mut_labels[j] for j in mut_order]
            isos_o = [short_isos[i] for i in iso_order]
            color_mat = np.zeros((n_iso, n_show, 3))
            color_mat[:] = 1.0
            for j_new, j_old in enumerate(mut_order):
                rgb = (0.18, 0.49, 0.20) if mut_curated[j_old] else (0.90, 0.32, 0.0)
                for i_new in range(n_iso):
                    if mat[i_new, j_new] > 0:
                        color_mat[i_new, j_new] = rgb
            fig_h = max(4, n_iso * 0.28)
            fig_w = max(8, min(18, n_show * 0.32 + 3))
            fig, ax = plt.subplots(figsize=(fig_w, fig_h))
            ax.imshow(color_mat, aspect="auto", interpolation="nearest")
            ax.set_yticks(range(n_iso))
            ax.set_yticklabels(isos_o, fontsize=7)
            ax.set_xticks(range(n_show))
            ax.set_xticklabels([_fig_short_label(x) for x in labels_o], rotation=90, fontsize=7)
            ax.set_title(f"Oncoprint: Mutation presence (top {n_show} by frequency)", fontweight="bold")
            legend_el = [
                Patch(facecolor=curated_color, label="Curated"),
                Patch(facecolor=noncurated_color, label="Non-curated"),
            ]
            ax.legend(handles=legend_el, loc="upper left", bbox_to_anchor=(1.01, 1))
            _save(fig, "01_oncoprint.png")
    except Exception as e:
        logger.error(f"  Oncoprint failed: {e}")

    try:
        if n_mut > 0:
            top_h = _fig_top_indices(mut_counts, FIG_TOP_MUT_HEAT)
            sub = mut_matrix[:, top_h]
            labels_h = [mut_labels[j] for j in top_h]
            fig_h = max(4, n_iso * 0.28)
            fig_w = max(8, min(18, len(top_h) * 0.32 + 3))
            fig, ax = plt.subplots(figsize=(fig_w, fig_h))
            cmap = ListedColormap([absent_color, present_color])
            ax.imshow(sub, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                      interpolation="nearest")
            ax.set_yticks(range(n_iso))
            ax.set_yticklabels(short_isos, fontsize=7)
            ax.set_xticks(range(len(top_h)))
            ax.set_xticklabels([_fig_short_label(x) for x in labels_h], rotation=90, fontsize=7)
            ax.set_title(f"Heatmap: Isolate × Mutation (top {len(top_h)}; blue = present)", fontweight="bold")
            _save(fig, "02_heatmap.png")
    except Exception as e:
        logger.error(f"  Heatmap failed: {e}")

    try:
        if n_mut > 0:
            order = _fig_top_indices(mut_counts, FIG_TOP_MUT_BAR)
            labels_f = [mut_labels[j] for j in order]
            counts_f = [mut_counts[j] for j in order]
            pcts = [100.0 * c / n_iso for c in counts_f]
            colors = [curated_color if mut_curated[j] else noncurated_color for j in order]
            fig, ax = plt.subplots(figsize=(max(8, min(16, len(order) * 0.4 + 2)), 6))
            ax.bar(range(len(order)), pcts, color=colors, edgecolor="white", width=0.8)
            ax.set_xticks(range(len(order)))
            ax.set_xticklabels([_fig_short_label(x) for x in labels_f], rotation=90, fontsize=7)
            ax.set_ylabel("% of isolates")
            ax.set_title(f"Frequency Bar: Top {len(order)} mutations (n={n_iso})", fontweight="bold")
            ax.legend(handles=[
                Patch(facecolor=curated_color, label="Curated"),
                Patch(facecolor=noncurated_color, label="Non-curated"),
            ], loc="upper right")
            ax.grid(axis="y", linestyle="--", alpha=0.3)
            _save(fig, "03_frequency_bar.png")
    except Exception as e:
        logger.error(f"  Frequency bar failed: {e}")

    try:
        burden = mut_matrix.sum(axis=1) if n_mut else np.zeros(n_iso)
        curated_burden = np.zeros(n_iso)
        noncur_burden = np.zeros(n_iso)
        for j in range(n_mut):
            if mut_curated[j]:
                curated_burden += mut_matrix[:, j]
            else:
                noncur_burden += mut_matrix[:, j]
        order = np.argsort(-burden)
        fig, ax = plt.subplots(figsize=(max(8, n_iso * 0.4), 5))
        x = range(n_iso)
        ax.bar(x, curated_burden[order], color=curated_color, label="Curated")
        ax.bar(x, noncur_burden[order], bottom=curated_burden[order],
               color=noncurated_color, label="Non-curated")
        ax.set_xticks(list(x))
        ax.set_xticklabels([short_isos[i] for i in order], rotation=90, fontsize=7)
        ax.set_ylabel("Mutation count")
        ax.set_title("Burden Plot: Mutations per isolate", fontweight="bold")
        ax.legend()
        ax.grid(axis="y", linestyle="--", alpha=0.3)
        _save(fig, "04_burden_plot.png")
    except Exception as e:
        logger.error(f"  Burden plot failed: {e}")

    try:
        gene_cur = {g: 0 for g in genes}
        gene_non = {g: 0 for g in genes}
        for j, g in enumerate(mut_gene):
            if mut_curated[j]:
                gene_cur[g] += mut_counts[j]
            else:
                gene_non[g] += mut_counts[j]
        high = mid = low = 0
        for c in mut_counts:
            pct = 100.0 * c / n_iso if n_iso else 0
            if pct >= 50:
                high += c
            elif pct >= 25:
                mid += c
            else:
                low += c
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[])
        ax.axis("off")
        ax.set_title("Sankey-style flow: Gene → Curated status → Prevalence class",
                     fontweight="bold", pad=20)
        col1_x, col2_x, col3_x = 0.08, 0.42, 0.78
        g_tot = {g: gene_cur[g] + gene_non[g] for g in genes}
        active_genes = [g for g in genes if g_tot[g] > 0] or list(genes)
        max_flow = max(max(g_tot.values()) if g_tot else 1, high, mid, low, 1)
        y_g = np.linspace(0.85, 0.15, len(active_genes))
        for yi, g in zip(y_g, active_genes):
            h = 0.04 + 0.08 * (g_tot[g] / max_flow)
            ax.add_patch(FancyBboxPatch(
                (col1_x - 0.06, yi - h / 2), 0.12, h,
                boxstyle="round,pad=0.01", facecolor="#90caf9", edgecolor="#1565c0"))
            ax.text(col1_x, yi, g, ha="center", va="center", fontsize=8, fontweight="bold")
        for label, y, color, val in [
            ("Curated", 0.65, curated_color, sum(gene_cur.values())),
            ("Non-curated", 0.35, noncurated_color, sum(gene_non.values())),
        ]:
            h = 0.05 + 0.1 * (val / max_flow)
            ax.add_patch(FancyBboxPatch(
                (col2_x - 0.08, y - h / 2), 0.16, h,
                boxstyle="round,pad=0.01", facecolor=color, edgecolor="black", alpha=0.85))
            ax.text(col2_x, y, f"{label}\n(n={val})", ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")
        for label, y, val in [
            ("High (≥50%)", 0.75, high),
            ("Mid (25–50%)", 0.50, mid),
            ("Low (<25%)", 0.25, low),
        ]:
            h = 0.04 + 0.08 * (val / max_flow)
            ax.add_patch(FancyBboxPatch(
                (col3_x - 0.08, y - h / 2), 0.16, h,
                boxstyle="round,pad=0.01", facecolor="#7e57c2", edgecolor="black", alpha=0.85))
            ax.text(col3_x, y, f"{label}\n(n={val})", ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")
        for yi, g in zip(y_g, active_genes):
            if gene_cur[g]:
                ax.annotate("", xy=(col2_x - 0.08, 0.65), xytext=(col1_x + 0.06, yi),
                            arrowprops=dict(arrowstyle="->", color=curated_color, alpha=0.4,
                                            lw=1 + 2 * gene_cur[g] / max_flow))
            if gene_non[g]:
                ax.annotate("", xy=(col2_x - 0.08, 0.35), xytext=(col1_x + 0.06, yi),
                            arrowprops=dict(arrowstyle="->", color=noncurated_color, alpha=0.4,
                                            lw=1 + 2 * gene_non[g] / max_flow))
        for y_src, val in [(0.65, sum(gene_cur.values())), (0.35, sum(gene_non.values()))]:
            for y_dst, dest_val in [(0.75, high), (0.50, mid), (0.25, low)]:
                if dest_val > 0 and val > 0:
                    ax.annotate("", xy=(col3_x - 0.08, y_dst), xytext=(col2_x + 0.08, y_src),
                                arrowprops=dict(arrowstyle="->", color="#7e57c2", alpha=0.25, lw=1.5))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        _save(fig, "05_sankey.png")
    except Exception as e:
        logger.error(f"  Sankey failed: {e}")

    try:
        if n_mut > 0:
            order = np.argsort(-np.array(mut_counts))[:min(40, n_mut)]
            theta = np.linspace(0, 2 * np.pi, len(order), endpoint=False)
            radii = [100.0 * mut_counts[j] / n_iso for j in order]
            colors = [curated_color if mut_curated[j] else noncurated_color for j in order]
            width = 2 * np.pi / len(order) * 0.85
            fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
            ax.bar(theta, radii, width=width, color=colors, edgecolor="white", alpha=0.9)
            ax.set_xticks(theta)
            ax.set_xticklabels([_fig_short_label(mut_labels[j], 16) for j in order], fontsize=7)
            ax.set_title("Circular: Mutation prevalence (top mutations)", fontweight="bold", pad=20)
            _save(fig, "06_circular.png")
    except Exception as e:
        logger.error(f"  Circular failed: {e}")

    try:
        gene_cur_c = []
        gene_non_c = []
        for g in genes:
            gc = sum(mut_counts[j] for j, gg in enumerate(mut_gene) if gg == g and mut_curated[j])
            gn = sum(mut_counts[j] for j, gg in enumerate(mut_gene) if gg == g and not mut_curated[j])
            gene_cur_c.append(gc)
            gene_non_c.append(gn)
        fig, ax = plt.subplots(figsize=(max(8, len(genes) * 0.7), 5))
        x = range(len(genes))
        ax.bar(x, gene_cur_c, color=curated_color, label="Curated (isolate hits)")
        ax.bar(x, gene_non_c, bottom=gene_cur_c, color=noncurated_color,
               label="Non-curated (isolate hits)")
        ax.set_xticks(list(x))
        ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=9)
        ax.set_ylabel("Sum of isolate×mutation detections")
        ax.set_title("Stacked Bar: Curated vs non-curated burden by gene", fontweight="bold")
        ax.legend()
        ax.grid(axis="y", linestyle="--", alpha=0.3)
        _save(fig, "07_stacked_bar.png")
    except Exception as e:
        logger.error(f"  Stacked bar failed: {e}")

    try:
        if n_mut > 0:
            fig, ax = plt.subplots(figsize=(12, 6))
            xs, ys, ss, cs, labs = [], [], [], [], []
            for j in range(n_mut):
                xs.append(j)
                ys.append(100.0 * mut_counts[j] / n_iso)
                ss.append(50 + 800 * (mut_counts[j] / n_iso))
                cs.append(curated_color if mut_curated[j] else noncurated_color)
                labs.append(mut_labels[j])
            ax.scatter(xs, ys, s=ss, c=cs, alpha=0.7, edgecolors="black", linewidths=0.4)
            ranked = sorted(zip(xs, ys, labs), key=lambda t: -t[1])[:FIG_TOP_MUT_ANNOT]
            for x, y, lab in ranked:
                if y >= 25:
                    ax.annotate(_fig_short_label(lab, 18), (x, y), textcoords="offset points",
                                xytext=(0, 6), ha="center", fontsize=7)
            ax.set_xlabel("Mutation index")
            ax.set_ylabel("% of isolates")
            ax.set_title("Bubble: Mutation prevalence (size ∝ frequency)", fontweight="bold")
            ax.legend(handles=[
                Patch(facecolor=curated_color, label="Curated"),
                Patch(facecolor=noncurated_color, label="Non-curated"),
            ])
            ax.grid(True, linestyle="--", alpha=0.3)
            _save(fig, "08_bubble.png")
    except Exception as e:
        logger.error(f"  Bubble failed: {e}")

    try:
        if n_mut > 0 and n_iso >= 2:
            fig, ax = plt.subplots(figsize=(max(8, n_iso * 0.35), 6))
            if SCIPY_OK and mut_matrix.sum() > 0:
                try:
                    dist = pdist(mut_matrix, metric="jaccard")
                    if np.all(np.isnan(dist)) or np.all(dist == 0):
                        dist = pdist(mut_matrix, metric="euclidean")
                    dist = np.nan_to_num(dist, nan=0.0)
                    Z = linkage(dist, method="average")
                    dendrogram(Z, labels=short_isos, leaf_rotation=90, leaf_font_size=8, ax=ax)
                    ax.set_title("Hierarchical clustering of "
                                 "mutation profiles", fontweight="bold")
                    ax.set_ylabel("Distance")
                except Exception:
                    order = np.argsort(-mut_matrix.sum(axis=1))
                    ax.bar(range(n_iso), mut_matrix.sum(axis=1)[order], color=present_color)
                    ax.set_xticks(range(n_iso))
                    ax.set_xticklabels([short_isos[i] for i in order], rotation=90, fontsize=7)
                    ax.set_title("Isolates ordered by observed mutation count",
                                 fontweight="bold")
            else:
                order = np.argsort(-mut_matrix.sum(axis=1))
                ax.bar(range(n_iso), mut_matrix.sum(axis=1)[order], color=present_color)
                ax.set_xticks(range(n_iso))
                ax.set_xticklabels([short_isos[i] for i in order], rotation=90, fontsize=7)
                ax.set_title("Observed mutation counts by isolate",
                             fontweight="bold")
            _save(fig, "09_phylogenetic_tree.png")
    except Exception as e:
        logger.error(f"  Phylogenetic tree failed: {e}")

    try:
        if len(genes) >= 2:
            gene_burden = np.zeros((n_iso, len(genes)))
            for j, g in enumerate(mut_gene):
                gi = genes.index(g)
                gene_burden[:, gi] += mut_matrix[:, j]
            with np.errstate(invalid="ignore"):
                if n_iso >= 3:
                    corr = np.corrcoef(gene_burden.T)
                    corr = np.nan_to_num(corr, nan=0.0)
                else:
                    corr = np.eye(len(genes))
            fig, ax = plt.subplots(figsize=(max(6, len(genes) * 0.7), max(5, len(genes) * 0.6)))
            im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
            ax.set_xticks(range(len(genes)))
            ax.set_yticks(range(len(genes)))
            ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
            ax.set_yticklabels(genes, fontsize=8)
            for i in range(len(genes)):
                for j in range(len(genes)):
                    ax.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center", fontsize=7)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            ax.set_title("Correlation Matrix: Gene mutation-burden co-occurrence",
                         fontweight="bold")
            _save(fig, "10_correlation_matrix.png")
    except Exception as e:
        logger.error(f"  Correlation matrix failed: {e}")

    try:
        burden = mut_matrix.sum(axis=1) if n_mut else np.zeros(n_iso)
        order = np.argsort(-burden)
        vals = burden[order]
        colors_w = [present_color if v > 0 else absent_color for v in vals]
        fig, ax = plt.subplots(figsize=(max(8, n_iso * 0.4), 5))
        ax.bar(range(n_iso), vals, color=colors_w, edgecolor="white")
        ax.set_xticks(range(n_iso))
        ax.set_xticklabels([short_isos[i] for i in order], rotation=90, fontsize=7)
        ax.set_ylabel("Total mutations")
        ax.set_title("Waterfall: Isolates ranked by mutation burden", fontweight="bold")
        ax.axhline(y=float(np.mean(vals)) if len(vals) else 0, color="red",
                   linestyle="--", alpha=0.6, label=f"Mean={np.mean(vals):.1f}")
        ax.legend()
        ax.grid(axis="y", linestyle="--", alpha=0.3)
        _save(fig, "11_waterfall.png")
    except Exception as e:
        logger.error(f"  Waterfall failed: {e}")

    try:
        if n_mut > 0:
            order = _fig_top_indices(mut_counts, FIG_TOP_MUT_BAR)[::-1]
            labels_l = [_fig_short_label(mut_labels[j], 28) for j in order]
            pcts = [100.0 * mut_counts[j] / n_iso for j in order]
            colors = [curated_color if mut_curated[j] else noncurated_color for j in order]
            fig, ax = plt.subplots(figsize=(10, max(4, len(order) * 0.32 + 1)))
            y = range(len(order))
            ax.hlines(y, 0, pcts, color="gray", linewidth=1, alpha=0.6)
            ax.scatter(pcts, y, s=60, c=colors, zorder=3, edgecolors="black", linewidths=0.4)
            for yi, p, c in zip(y, pcts, [mut_counts[j] for j in order]):
                ax.text(p + 1, yi, f"{p:.0f}% ({c}/{n_iso})", va="center", fontsize=7)
            ax.set_yticks(list(y))
            ax.set_yticklabels(labels_l, fontsize=7)
            ax.set_xlabel("% of isolates")
            ax.set_title(f"Lollipop: Top {len(order)} mutations by prevalence", fontweight="bold")
            ax.legend(handles=[
                Patch(facecolor=curated_color, label="Curated"),
                Patch(facecolor=noncurated_color, label="Non-curated"),
            ], loc="lower right")
            ax.set_xlim(0, max(pcts) * 1.25 if pcts else 1)
            ax.grid(axis="x", linestyle="--", alpha=0.3)
            _save(fig, "12_lollipop.png")
    except Exception as e:
        logger.error(f"  Lollipop failed: {e}")

    try:
        fig, ax = plt.subplots(figsize=(max(6, len(genes) * 0.7), max(4, n_iso * 0.28)))
        cmap = ListedColormap(["#ef9a9a", "#a5d6a7"])
        ax.imshow(gene_present, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                  interpolation="nearest")
        ax.set_yticks(range(n_iso))
        ax.set_yticklabels(short_isos, fontsize=7)
        ax.set_xticks(range(len(genes)))
        ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
        ax.set_title("Gene Presence Heatmap (green=found, red=missing)", fontweight="bold")
        legend_el = [
            Patch(facecolor="#a5d6a7", label="Gene found"),
            Patch(facecolor="#ef9a9a", label="Gene missing"),
        ]
        ax.legend(handles=legend_el, loc="upper left", bbox_to_anchor=(1.01, 1))
        _save(fig, "13_gene_presence_heatmap.png")
    except Exception as e:
        logger.error(f"  Gene presence heatmap failed: {e}")

    logger.info(f"\nMaster analysis figures → {out_dir}")
    return out_dir
