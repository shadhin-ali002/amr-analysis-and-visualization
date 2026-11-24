"""Grouped mutation reports by gene and by isolate."""
from __future__ import annotations

import json
import os
import re
import tempfile as _tempfile
from collections import defaultdict
from pathlib import Path

import openpyxl

from ng_amr.core.constants import (
    AD_MIC,
    ANTIBIOTIC_ORDER,
    CODING,
    ORIGINAL,
    QC,
    ROUNDED,
)
from ng_amr.genomics.alignment import (
    A4,
    REPORTLAB_OK,
    analyze_nucleotide,
    analyze_protein,
    rl_canvas,
    safe_name,
)
from ng_amr.genomics.locus import (
    _CACHE,
)
from ng_amr.reports.html import (
    _inject_disruption_banner,
    create_combined_html_alignment,
)
from ng_amr.reports.pdf import (
    _draw_combined_pdf,
    _PDFPage,
)


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
    from textwrap import wrap

    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter
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
