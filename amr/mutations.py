"""Mutations orchestration."""
from __future__ import annotations

from .engine import (
    DEFAULT_GENES,
    QC,
    SeqIO,
    _CACHE,
    _parasail,
    clinical_gene_profiles,
    clinical_view,
    combined_master,
    create_gene_mutation_cooccurrence_figure,
    create_master_analysis_figures,
    create_multigene_mutation_cooccurrence_count_heatmap,
    csv_write,
    generate_aa_comparison_package,
    generate_exclusive_figures_package,
    generate_extended_thesis_publication_package,
    generate_extracted_sequences_package,
    generate_final_for_thesis_package,
    generate_thesis_publication_package,
    grouped_mutation_reports,
    load_curated,
    process_gene,
    shutil,
    stage,
    write_master_matrix,
)


def genomes(base,assemblies,book,root,args,audit):
    refs=base/"Refs"/"Genes_RefSeq.fasta"
    catalog=args.catalog or base/"Refs"/"Curated_Variants.json"
    reasons=[]
    if not any(a.get("plane") for a in assemblies.values()):reasons.append("no unambiguous multi-contig assemblies")
    for path in (refs,catalog):
        if not path.is_file():reasons.append("missing "+str(path))
    for name,installed in (("biopython",SeqIO),("parasail",_parasail),("blastn",shutil.which("blastn")),("makeblastdb",shutil.which("makeblastdb"))):
        if not installed:reasons.append("missing "+name)
    if reasons:
        audit.add("WGS","skipped","; ".join(reasons))
        return
    _CACHE.update(hits={},contigs={},refs=None)
    curated=load_curated(str(catalog))
    genes=args.genes.split(",") if args.genes else DEFAULT_GENES
    results=[]
    for gene in dict.fromkeys(g.strip() for g in genes):
        if gene not in DEFAULT_GENES:
            audit.add("WGS gene","failed",f"unknown gene {gene}");continue
        gr=stage(audit,"WGS "+gene,lambda:process_gene(gene,str(refs),curated,assemblies,str(root),args.thresholds))
        if gr:
            results.append(gr)
            for r in gr["all_data"]:
                if r.get("status","").startswith("assembly failed"):audit.add("WGS isolate","failed",f"{r['isolate']} {gene}: {r['status']}")
    if not results:return
    stage(audit, "grouped mutation reports", lambda: grouped_mutation_reports(results, book, root, not args.no_html))
    cases=clinical_view(results)
    csv_write(root/"WGS_calls.csv",["Canonical ID","Cohort","Gene","Assembly status","Reference identity (%)","Query coverage (%)","Variant count when called"],
              ([r["isolate"],"QC" if r["isolate"]==QC else "clinical",g["gene"],r.get("status"),r.get("pident"),r.get("qcov"),len(r["mutations"]) if r.get("gene_found") else None] for g in results for r in g["all_data"]))
    stage(audit,"clinical mutation and MIC matrix",lambda:write_master_matrix(cases,str(root),args.thresholds,book))
    if cases and any(gr["all_data"] for gr in cases):
        stage(audit,"combined master matrix",lambda:combined_master(root))
    stage(audit,"extracted sequences",lambda:generate_extracted_sequences_package(results,str(root)))
    aa=base/"Assembly"/"AA"
    if aa.is_dir():
        stage(audit,"AA comparison",lambda:generate_aa_comparison_package(results,str(root),str(aa)))
    else:audit.add("AA comparison","skipped","original AA directory missing")
    stage(audit,"clinical genomic and MIC profiles",lambda:clinical_gene_profiles(cases,book,root,audit))
    all_ids=sorted({d["isolate"] for g in cases for d in g["all_data"]})
    if not all_ids:return cases
    # Full data drive the tables. Bound the source's isolate-axis illustration
    # galleries, and make the displayed cohort explicit in the folder and manifest.
    plotroot=root
    visual=cases
    if len(all_ids)>60:
        selected=all_ids[:60];visual=clinical_view(results,set(selected))
        plotroot=root/f"Clinical_figure_subset_60_of_{len(all_ids)}"
        plotroot.mkdir(exist_ok=True)
        csv_write(plotroot/"displayed_isolates.csv",["Canonical ID"],([k] for k in selected))
        audit.add("WGS figures","subset",f"first 60 of {len(all_ids)} clinical IDs in lexical order in {plotroot.name}")
    for name,fn in (("master genomic figures",create_master_analysis_figures),
                    ("gene co-occurrence",create_gene_mutation_cooccurrence_figure),
                    ("multi-gene co-occurrence",create_multigene_mutation_cooccurrence_count_heatmap)):
        stage(audit,name,lambda:fn(visual,str(plotroot)))
    stage(audit,"exclusive genomic figures",lambda:generate_exclusive_figures_package(visual,str(plotroot),curated))
    for name,fn in (("publication figures",generate_thesis_publication_package),
                    ("extended publication figures",generate_extended_thesis_publication_package)):
        stage(audit,name,lambda:fn(visual,str(plotroot),curated,args.thresholds))
    stage(audit,"thesis tables and reports",lambda:generate_final_for_thesis_package(cases,str(root),curated,args.thresholds))
    return cases
