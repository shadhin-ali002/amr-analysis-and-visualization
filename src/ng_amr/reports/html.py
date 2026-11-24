"""HTML mutation and alignment reports."""
from __future__ import annotations

import os

from ng_amr.genomics.alignment import (
    categorize_nt_mutations,
)


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
