"""Run summary orchestration."""
from __future__ import annotations

from .engine import (
    DEFAULT_GENES,
    Path,
    build_run_summary_pptx,
    hashlib,
    json,
    sys,
)


def report_run(root, args, results, audit, started):
    from datetime import datetime
    now = datetime.now().astimezone()
    elapsed = (now - started).total_seconds()
    selected = args.genes.split(",") if args.genes else DEFAULT_GENES
    gene_rows = [{"gene": g["gene"], "processed": g["processed"], "missing": g["skipped"],
                  "muts": len(g["mutation_counts"])} for g in results]
    summary = dict(start_time=started.isoformat(), end_time=now.isoformat(), elapsed_sec=elapsed,
                   elapsed_human=str(now - started), output_root=str(root), genes_requested=len(selected),
                   genes_ok=len(results), genes_failed=len(selected)-len(results), gene_rows=gene_rows,
                   total_distinct_muts=sum(r["muts"] for r in gene_rows),
                   total_events=sum(sum(g["mutation_counts"].values()) for g in results),
                   thresholds_str=", ".join(map(str, args.thresholds)), write_html=not args.no_html,
                   mic_used=any(name == "master workbook" and status == "ok" for name,status,_ in audit.rows),
                   status="See stage outcomes", stages=[{"name":name,"status":status.upper(),"note":detail}
                            for name,status,detail in audit.rows if name in
                            {"master workbook","genomic workflow","genome phylogeny","virulence reporting",
                             "E-strip photo workflow","AD photo workflow","MIC illustration galleries"}])
    manifest = dict(arguments=vars(args), started=started.isoformat(), finished=now.isoformat(),
                    python=sys.version, source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in [Path(__file__).resolve().parents[1] / "run_amr.py", *sorted(Path(__file__).parent.glob("*.py"))]},
                    workbook_sha256=hashlib.sha256(args.workbook.read_bytes()).hexdigest() if args.workbook is not None else None,
                    stages=audit.rows)
    other = root / "05_Other_Outputs"
    (other / "NG_AMR_Run_Manifest.json").write_text(json.dumps(manifest, default=str, indent=2), encoding="utf-8")
    if not build_run_summary_pptx(str(other), summary):
        raise RuntimeError("Could not create run summary presentation")
