"""Genomic workflow: per-gene calling, then reports, matrices, sequences and summary tables."""
from __future__ import annotations

import shutil

try:
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
    from Bio.Seq import Seq
except ImportError:
    SeqIO = Seq = PairwiseAligner = None

try:
    import parasail as _parasail
except ImportError:
    _parasail = None

from ng_amr.core.audit import (
    csv_write,
    stage,
)
from ng_amr.core.constants import (
    DEFAULT_GENES,
    QC,
)
from ng_amr.genomics.aa_comparison import (
    generate_aa_comparison_package,
)
from ng_amr.genomics.alignment import (
    load_curated,
)
from ng_amr.genomics.locus import (
    _CACHE,
    process_gene,
)
from ng_amr.genomics.sequences import (
    generate_extracted_sequences_package,
)
from ng_amr.reports.grouped import (
    grouped_mutation_reports,
)
from ng_amr.reports.master_matrix import (
    clinical_view,
    combined_master,
    write_master_matrix,
)
from ng_amr.reports.summary_tables import (
    generate_summary_tables,
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
    all_ids=sorted({d["isolate"] for g in cases for d in g["all_data"]})
    if not all_ids:return cases
    stage(audit,"summary tables",lambda:generate_summary_tables(cases,str(root),curated,args.thresholds))
    return cases
