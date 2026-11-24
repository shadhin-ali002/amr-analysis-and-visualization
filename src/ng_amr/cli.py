"""Command-line interface: MIC summaries and, optionally, genomic mutation reports."""
from __future__ import annotations

import argparse
import csv
import logging
from datetime import datetime
from pathlib import Path

from ng_amr import __version__


def assemblies(args, audit):
    from ng_amr.core.ids import canonical
    from ng_amr.io.assemblies import discover_assemblies

    if args.wgs_base is None:
        return {}
    mapping = {}
    if args.id_map:
        with args.id_map.open(encoding="utf-8-sig", newline="") as handle:
            rows = csv.DictReader(handle)
            if not {"filename", "isolate_id"} <= set(rows.fieldnames or []):
                raise ValueError("Expected filename and isolate_id in the ID map")
            for row in rows:
                name = row["filename"]
                sample = canonical(row["isolate_id"])
                if not name or not sample or name in mapping:
                    raise ValueError("Empty or duplicate filename or empty sample ID")
                mapping[name] = sample
    return discover_assemblies(args.wgs_base, audit, mapping)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="ng-amr",
        description="MIC workbook summaries and resistance-gene mutation reports for Neisseria gonorrhoeae.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--workbook", type=Path, help="MIC workbook")
    parser.add_argument("--wgs-base", type=Path, help="Parent directory of Assembly and Refs")
    parser.add_argument("--catalog", type=Path, help="Override the curated JSON filename")
    parser.add_argument("--id-map", type=Path, help="Optional filename,isolate_id CSV")
    parser.add_argument("--aliases", type=Path, help="Private alias,canonical_id CSV")
    parser.add_argument("--output", type=Path, default=Path("outputs"), help="Parent of timestamped run folders")
    parser.add_argument("--genes", help="Comma-separated gene subset")
    parser.add_argument("--thresholds", default="25,30,35,40,50", help="Comma-separated reporting percentages")
    parser.add_argument("--no-html", action="store_true", help="Skip HTML alignment reports")
    parser.add_argument("--self-test", action="store_true", help="Run the built-in regression checks and exit")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.self_test:
        from ng_amr.selftest import run_self_tests

        logging.basicConfig(level=logging.INFO, format="%(message)s")
        verdict = run_self_tests()
        return 2 if verdict is None else int(not verdict)
    try:
        args.thresholds = sorted({int(value) for value in args.thresholds.split(",")})
        if not all(1 <= value <= 100 for value in args.thresholds):
            raise ValueError
    except ValueError:
        parser.error("Thresholds must be whole percentages from 1 to 100")
    if not (args.workbook or args.wgs_base):
        parser.error("Supply --workbook or --wgs-base")
    for name, value in vars(args).items():
        if name == "output" or not isinstance(value, Path):
            continue
        if name in {"wgs_base", "est_raw"}:
            if not value.is_dir() or args.output.resolve().is_relative_to(value.resolve()):
                parser.error(f"Supply a valid input folder distinct from the output tree: {value}")
        elif not value.is_file():
            parser.error(f"Missing input file: {value}")

    from ng_amr.core.audit import Audit, csv_write
    from ng_amr.core.ids import ALIASES, load_aliases

    try:
        aliases = load_aliases(args.aliases)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    ALIASES.clear()
    ALIASES.update(aliases)
    started = datetime.now().astimezone()
    root = args.output / ("Run_" + started.strftime("%Y%m%d_%H%M%S_%f"))
    root.mkdir(parents=True)
    other = root / "05_Other_Outputs"
    other.mkdir()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(other / "analysis.log", encoding="utf-8")])
    audit = Audit()
    try:
        from ng_amr.genomics.pipeline import genomes
        from ng_amr.io.workbook import Phenotypes, read_master
        from ng_amr.phenotype.mic import phenotype_outputs
        from ng_amr.reports.run_summary import report_run

        book = read_master(args.workbook, audit) if args.workbook else Phenotypes()
        if args.workbook:
            audit.add("master workbook", "ok", str(args.workbook))
        samples = assemblies(args, audit)
        folder = root / "01_Mutation_Analysis"
        folder.mkdir()
        if args.workbook or samples:
            phenotype_outputs(book, samples, folder, audit)
        results = genomes(args.wgs_base, samples, book, folder, args, audit) if args.wgs_base else []
        report_run(root, args, results or [], audit, started)
        if args.wgs_base and not results:
            audit.add("genomic workflow", "failed", "No genomic reports produced; inspect WGS status")
    except Exception as exc:
        logging.exception("Workflow error")
        audit.add("run", "failed", str(exc))
    finally:
        csv_write(other / "run_status.csv", ["Stage", "Status", "Detail"], audit.rows)
    print(root)
    return int(any(status == "failed" for _, status, _ in audit.rows))


if __name__ == "__main__":
    raise SystemExit(main())
