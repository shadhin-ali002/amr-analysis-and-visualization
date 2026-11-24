"""Cross-gene summary tables: mutation catalogue, burden, thresholds, drug classes and key findings."""
from __future__ import annotations

import csv
import os
import re
from typing import Dict, List, Optional

import pandas as pd

from ng_amr.core.audit import (
    logger,
)

DRUG_CLASS = {
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





def _summary_copy_master_matrix(output_root: str, tables_root: str) -> Optional[str]:
    src = os.path.join(output_root, "NG_AMR_Mutation_MIC_Matrix.xlsx")
    if not os.path.isfile(src):return None
    dst = os.path.join(tables_root, "00_NG_AMR_Master_Mutation_Matrix.xlsx")
    try:
        import shutil
        shutil.copy2(src, dst)
        logger.info(f"  [Summary_Tables] Copied master matrix → {os.path.basename(dst)}")
        return dst
    except Exception as e:
        logger.error(f"  [Summary_Tables] Copy master matrix failed: {e}")
        return None


def _summary_aggregate_rows(gene_results, curated_master=None):
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
            "drug_class": DRUG_CLASS.get(gene, "Other / context"),
        })
    return rows


def _summary_write_mutation_catalogue(rows, tables_root, thresholds):
    out_dir = os.path.join(tables_root, "01_Mutation_Catalogue")
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
    logger.info(f"  [Summary_Tables] Mutation catalogue → {out_dir}")
    return out_dir


def _summary_burden_and_carriage(rows, tables_root):
    out_dir = os.path.join(tables_root, "02_Mutation_Burden")
    os.makedirs(out_dir, exist_ok=True)
    summ_path = os.path.join(out_dir, "Burden_summary.csv")
    with open(summ_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Gene", "Isolates", "Mean_burden", "Pct_with_ge1", "Distinct", "Events"])
        for r in rows:
            w.writerow([r["gene"], r["total"], f"{r['mean_burden']:.3f}",
                        f"{r['pct_with']:.2f}", r["distinct"], r["events"]])
    logger.info(f"  [Summary_Tables] Burden analysis → {out_dir}")
    return out_dir


def _summary_thresholds(rows, tables_root, thresholds):
    out_dir = os.path.join(tables_root, "03_Threshold_Stratification")
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
    logger.info(f"  [Summary_Tables] Threshold stratification → {out_dir}")
    return out_dir


def _summary_curated_vs_non(rows, tables_root):
    out_dir = os.path.join(tables_root, "04_Curated_vs_NonCurated")
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
    logger.info(f"  [Summary_Tables] Curated vs non-curated → {out_dir}")
    return out_dir


def _summary_drug_class_summary(rows, tables_root):
    out_dir = os.path.join(tables_root, "05_DrugClass_Summary")
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
    logger.info(f"  [Summary_Tables] Drug-class summary → {out_dir}")
    return out_dir


def _summary_key_findings(rows, tables_root, thresholds):
    out_dir = os.path.join(tables_root, "06_Key_Findings_and_Index")
    os.makedirs(out_dir, exist_ok=True)
    lines = []
    lines.append("# Summary tables — key findings")
    lines.append("")
    lines.append("Generated automatically by ng-amr.")
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
    lines.append("See `03_Threshold_Stratification/` for mutations meeting each threshold.")
    lines.append("")
    lines.append("## Folder index")
    lines.append("| Folder | Content |")
    lines.append("|--------|---------|")
    lines.append("| `00_NG_AMR_Master_Mutation_Matrix.xlsx` | Copy of master matrix |")
    lines.append("| `01_Mutation_Catalogue` | Full mutation catalogue (CSV/XLSX) |")
    lines.append("| `02_Mutation_Burden` | Burden summary per gene |")
    lines.append("| `03_Threshold_Stratification` | ≥25–50% prevalence strata |")
    lines.append("| `04_Curated_vs_NonCurated` | AMR-relevant vs other mutations |")
    lines.append("| `05_DrugClass_Summary` | Grouped by antibiotic class |")
    lines.append("| `06_Key_Findings_and_Index` | This narrative + index |")
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
    logger.info(f"  [Summary_Tables] Key findings → {path}")
    return out_dir


def _summary_advanced_genomics_modules(rows, tables_root, gene_results, output_root):
    d07 = os.path.join(tables_root, "07_Typing_MLST_NGSTAR")
    os.makedirs(d07, exist_ok=True)
    readme10 = os.path.join(d07, "README_Typing_requirements.md")
    with open(readme10, "w", encoding="utf-8") as fh:
        fh.write("""# Molecular typing (MLST / NG-MAST / NG-STAR)

## Why this folder exists
WHO EGASP and published NG WGS studies report MLST, NG-MAST and NG-STAR alongside AMR determinants.

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
    with open(os.path.join(d07, "Typing_template.csv"), "w", newline="", encoding="utf-8") as fh:
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
    logger.info(f"  [Summary_Tables] Typing module (templates) → {d07}")

    d08 = os.path.join(tables_root, "08_Phylogeny")
    os.makedirs(d08, exist_ok=True)
    with open(os.path.join(d08, "README_Phylogeny_requirements.md"), "w", encoding="utf-8") as fh:
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
    logger.info(f"  [Summary_Tables] Phylogeny module (checklist) → {d08}")

    d09 = os.path.join(tables_root, "09_Plasmid_MGE")
    os.makedirs(d09, exist_ok=True)
    with open(os.path.join(d09, "README_Plasmid_requirements.md"), "w", encoding="utf-8") as fh:
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
        with open(os.path.join(d09, "blaTEM_tetM_detection_from_pipeline.csv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(plasmid_rows[0].keys()))
            w.writeheader()
            w.writerows(plasmid_rows)
    logger.info(f"  [Summary_Tables] Plasmid module → {d09}")

    d10 = os.path.join(tables_root, "10_Structural_penA_PBP2")
    os.makedirs(d10, exist_ok=True)
    PBP2_MOTIFS = {
        "SxxK_active_site": list(range(310, 320)),
        "SxN_motif": list(range(360, 366)),
        "KTG_motif": list(range(497, 505)),
        "beta3_beta4_loop": list(range(500, 520)),
        "alpha2_helix_cluster": list(range(480, 500)),
        "Ala501_cluster": [501],
        "Gly545_region": list(range(542, 548)),
    }
    with open(os.path.join(d10, "README_Structural_mapping.md"), "w", encoding="utf-8") as fh:
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
        with open(os.path.join(d10, "penA_mutations_structural_motif_map.csv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(mapped[0].keys()))
            w.writeheader()
            w.writerows(mapped)
    logger.info(f"  [Summary_Tables] Structural penA mapping → {d10}")
    return d07, d08, d09, d10


def generate_summary_tables(
    gene_results: List[dict],
    output_root: str,
    curated_master: Optional[dict] = None,
    thresholds: Optional[List[int]] = None,
    mic_excel: Optional[str] = None,
) -> str:
    logger.info("\n" + "#" * 70)
    logger.info("# SUMMARY TABLES")
    logger.info("#" * 70)
    tables_root = os.path.join(output_root, "Summary_Tables")
    os.makedirs(tables_root, exist_ok=True)

    rows = _summary_aggregate_rows(gene_results, curated_master)
    if not rows:
        logger.info("  [Summary_Tables] No gene_results — skip")
        return tables_root

    _summary_copy_master_matrix(output_root, tables_root)
    _summary_write_mutation_catalogue(rows, tables_root, thresholds)
    _summary_burden_and_carriage(rows, tables_root)
    _summary_thresholds(rows, tables_root, thresholds)
    _summary_curated_vs_non(rows, tables_root)
    _summary_drug_class_summary(rows, tables_root)
    _summary_key_findings(rows, tables_root, thresholds)

    try:
        _summary_advanced_genomics_modules(rows, tables_root, gene_results, output_root)
    except Exception as e:
        logger.error(f"  [Summary_Tables] Advanced genomics modules failed: {e}")

    try:
        idx_path = os.path.join(tables_root, "06_Key_Findings_and_Index", "KEY_FINDINGS.md")
        extra = (
            "\n\n## Also in this package\n"
            "- `07_Typing_MLST_NGSTAR/` — templates + checklist for external typing\n"
            "- `08_Phylogeny/` — phylogeny checklist\n"
            "- `09_Plasmid_MGE/` — plasmid module + blaTEM/tetM detection if analysed\n"
            "- `10_Structural_penA_PBP2/` — penA mutations mapped to PBP2 motifs\n"
        )
        if os.path.isfile(idx_path):
            with open(idx_path, "a", encoding="utf-8") as fh:
                fh.write(extra)
    except Exception as _exc:
        logger.debug("Ignored non-critical exception: %s", _exc)
        pass

    logger.info("\nSummary tables complete → %s", tables_root)
    return tables_root


