"""Sequence translation, pairwise alignment and mutation calling."""
from __future__ import annotations

import json
import re
from typing import Dict, List, Optional, Tuple

try:
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
    from Bio.Seq import Seq
except ImportError:
    SeqIO = Seq = PairwiseAligner = None

from ng_amr.core.constants import (
    DISRUPTION_MAX_LENGTH_CHANGE,
    DISRUPTION_MIN_DELETION_AA,
    DISRUPTION_MIN_DELETION_EVENTS,
    GENE_CONFIG,
)


def safe_name(value):
    return re.sub(r'[^\w. -]', '_', str(value)).strip('. ') or 'unnamed'


try:
    from reportlab.lib.colors import Color, HexColor, black, white
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import inch, mm
    from reportlab.pdfgen import canvas as rl_canvas
    REPORTLAB_OK = True
except ImportError:
    REPORTLAB_OK = False
    A4 = (595.2755905511812, 841.8897637795277)
    inch = 72.0
    mm = 2.834645669291339
    def HexColor(_x):
        return None
    black = white = Color = None
    rl_canvas = None


def translate_dna_to_protein(dna_seq: str) -> str:
    return str(Seq(dna_seq).translate(to_stop=False))


def _pairwise_global(a: str, b: str, open_gap: float = -10) -> Tuple[str, str]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(a, b)
    if not alignments:
        return a, b
    aln = alignments[0]
    a_out, b_out = [], []
    a_prev = b_prev = 0
    for (a_s, a_e), (b_s, b_e) in zip(aln.aligned[0], aln.aligned[1]):
        while a_prev < a_s:
            a_out.append(a[a_prev]); b_out.append("-"); a_prev += 1
        while b_prev < b_s:
            a_out.append("-"); b_out.append(b[b_prev]); b_prev += 1
        for x, y in zip(a[a_s:a_e], b[b_s:b_e]):
            a_out.append(x); b_out.append(y)
        a_prev, b_prev = a_e, b_e
    while a_prev < len(a):
        a_out.append(a[a_prev]); b_out.append("-"); a_prev += 1
    while b_prev < len(b):
        a_out.append("-"); b_out.append(b[b_prev]); b_prev += 1
    return "".join(a_out), "".join(b_out)


def align_protein_sequences(ref_prot: str, qry_prot: str,
                            open_gap: float = -5) -> Tuple[str, str]:
    return _pairwise_global(ref_prot, qry_prot, open_gap=open_gap)


def call_protein_mutations(ref_aln: str, qry_aln: str) -> List[dict]:
    mutations = []
    ref_pos = 0
    i = 0
    while i < len(ref_aln):
        r, q = ref_aln[i], qry_aln[i]
        if r != "-" and q != "-":
            ref_pos += 1
            if r != q:
                if q == "*":
                    mutations.append({
                        "type": "stop_gain", "position": ref_pos,
                        "ref_aa": r, "alt_aa": q, "notation": f"{r}{ref_pos}Ter",
                    })
                else:
                    mutations.append({
                        "type": "substitution", "position": ref_pos,
                        "ref_aa": r, "alt_aa": q, "notation": f"{r}{ref_pos}{q}",
                    })
            i += 1
        elif r == "-" and q != "-":
            inserted = []
            start_pos = ref_pos
            while i < len(ref_aln) and ref_aln[i] == "-" and qry_aln[i] != "-":
                inserted.append(qry_aln[i]); i += 1
            aa_seq = "".join(inserted)
            ins_pos = start_pos + 1
            notation = (f"ins{ins_pos}{aa_seq}" if len(aa_seq) == 1
                        else f"ins{ins_pos}_{ins_pos + len(aa_seq) - 1}{aa_seq}")
            mutations.append({
                "type": "insertion", "position": ins_pos, "ref_position": start_pos,
                "aa_seq": aa_seq, "notation": notation, "length": len(aa_seq),
            })
        elif r != "-" and q == "-":
            deleted = []
            start_pos = ref_pos
            while i < len(ref_aln) and ref_aln[i] != "-" and qry_aln[i] == "-":
                deleted.append(ref_aln[i]); ref_pos += 1; i += 1
            end_pos = ref_pos - 1
            aa_seq = "".join(deleted)
            notation = (f"del{start_pos}" if len(aa_seq) == 1
                        else f"del{start_pos}_{end_pos}")
            mutations.append({
                "type": "deletion", "position": start_pos, "end_position": end_pos,
                "aa_seq": aa_seq, "notation": notation, "length": len(aa_seq),
            })
        else:
            i += 1
    return mutations


def detect_frameshifts(protein_mutations: List[dict], nt_mutations: List[dict],
                       window: int = 12) -> List[dict]:
    nt_indels = [m for m in (nt_mutations or [])
                 if m.get("type") in ("insertion", "deletion")]
    if not nt_indels:
        return []
    nt_indels = sorted(nt_indels, key=lambda m: (m.get("position") or 0))

    stop_positions = sorted(
        m["position"] for m in (protein_mutations or [])
        if m.get("type") == "stop_gain" and m.get("position") is not None
    )

    frameshifts: List[dict] = []
    net = 0
    open_fs: Optional[dict] = None

    for nm in nt_indels:
        length = int(nm.get("length") or 0)
        if length <= 0:
            continue
        nt_pos = int(nm.get("position") or 1)
        net += length if nm["type"] == "insertion" else -length

        if net % 3 != 0:
            if open_fs is None:
                aa_pos = max(1, (nt_pos + 2) // 3)
                stop_pos = next((s for s in stop_positions if s > aa_pos), None)
                open_fs = {
                    "type": f"frameshift_{nm['type']}",
                    "position": aa_pos,
                    "nt_position": nt_pos,
                    "length": abs(net),
                    "aa_length": max(1, abs(net) // 3),
                    "net_shift": net,
                    "notation": f"fs{aa_pos}Ter{stop_pos}" if stop_pos else f"fs{aa_pos}",
                }
                if stop_pos:
                    open_fs["stop_position"] = stop_pos
                if nm["type"] == "deletion":
                    end_nt = int(nm.get("end_position") or nt_pos)
                    open_fs["end_position"] = max(open_fs["position"], (end_nt + 2) // 3)
                frameshifts.append(open_fs)
            else:
                open_fs["length"] = abs(net)
                open_fs["aa_length"] = max(1, abs(net) // 3)
                open_fs["net_shift"] = net
        else:
            if open_fs is not None:
                nt_open = int(open_fs.get("nt_position") or 0)
                if abs(nt_pos - nt_open) <= window and open_fs in frameshifts:
                    frameshifts.remove(open_fs)
            open_fs = None

    return frameshifts


def analyze_protein(gene_name: str, ref_dna: str, query_dna: str,
                    nt_mutations: Optional[List[dict]] = None,
                    open_gap: float = -5) -> dict:
    ref_prot = translate_dna_to_protein(ref_dna)
    qry_prot = translate_dna_to_protein(query_dna)
    ref_aln, qry_aln = align_protein_sequences(ref_prot, qry_prot, open_gap=open_gap)
    mutations = call_protein_mutations(ref_aln, qry_aln)
    frameshifts = detect_frameshifts(mutations, nt_mutations or [])
    silent: List[dict] = []
    all_muts = mutations + frameshifts + silent
    seen = set()
    unique = []
    for m in all_muts:
        if m["notation"] not in seen:
            seen.add(m["notation"])
            unique.append(m)
    return {
        "gene": gene_name,
        "ref_protein": ref_prot,
        "query_protein": qry_prot,
        "ref_alignment": ref_aln,
        "query_alignment": qry_aln,
        "mutations": unique,
        "substitutions": [m for m in unique if m["type"] == "substitution"],
        "insertions": [m for m in unique if m["type"] == "insertion"],
        "deletions": [m for m in unique if m["type"] == "deletion"],
        "stop_gains": [m for m in unique if m["type"] == "stop_gain"],
        "frameshifts": frameshifts,
        "silent_mutations": silent,
    }


def analyze_nucleotide(ref_dna, qry_dna):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(ref_dna, qry_dna)
    if not alignments:
        ref_aln = ref_dna
        qry_aln = qry_dna
    else:
        alignment = alignments[0]
        ref_aln_list = []
        qry_aln_list = []
        r_blocks = alignment.aligned[0]
        q_blocks = alignment.aligned[1]
        r_prev = 0
        q_prev = 0
        for (r_start, r_end), (q_start, q_end) in zip(r_blocks, q_blocks):
            while r_prev < r_start:
                ref_aln_list.append(ref_dna[r_prev])
                qry_aln_list.append("-")
                r_prev += 1
            while q_prev < q_start:
                ref_aln_list.append("-")
                qry_aln_list.append(qry_dna[q_prev])
                q_prev += 1
            for r, q in zip(ref_dna[r_start:r_end], qry_dna[q_start:q_end]):
                ref_aln_list.append(r)
                qry_aln_list.append(q)
            r_prev = r_end
            q_prev = q_end
        while r_prev < len(ref_dna):
            ref_aln_list.append(ref_dna[r_prev])
            qry_aln_list.append("-")
            r_prev += 1
        while q_prev < len(qry_dna):
            ref_aln_list.append("-")
            qry_aln_list.append(qry_dna[q_prev])
            q_prev += 1
        ref_aln = "".join(ref_aln_list)
        qry_aln = "".join(qry_aln_list)
    ref_pos_map = []
    ref_counter = 0
    for nt in ref_aln:
        if nt != '-':
            ref_counter += 1
            ref_pos_map.append(ref_counter)
        else:
            ref_pos_map.append(None)
    mutations = []
    i = 0
    while i < len(ref_aln):
        r = ref_aln[i]
        q = qry_aln[i]
        if r != '-' and q != '-':
            if r != q:
                pos = ref_pos_map[i]
                mutations.append({
                    'type': 'SNP',
                    'position': pos,
                    'ref_nt': r,
                    'alt_nt': q,
                    'notation': f"{r}{pos}{q}"
                })
            i += 1
        elif r == '-' and q != '-':
            inserted = []
            pos_before = None
            for j in range(i-1, -1, -1):
                if ref_pos_map[j] is not None:
                    pos_before = ref_pos_map[j]
                    break
            if pos_before is None:
                pos_before = 0
            while i < len(ref_aln) and ref_aln[i] == '-' and qry_aln[i] != '-':
                inserted.append(qry_aln[i])
                i += 1
            alt_seq = ''.join(inserted)
            ins_pos = pos_before + 1
            if len(alt_seq) == 1:
                notation = f"ins{ins_pos}{alt_seq}"
            else:
                notation = f"ins{ins_pos}_{ins_pos+len(alt_seq)-1}{alt_seq}"
            mutations.append({
                'type': 'insertion',
                'position': ins_pos,
                'ref_nt': '-',
                'alt_nt': alt_seq,
                'notation': notation,
                'length': len(alt_seq)
            })
        elif r != '-' and q == '-':
            deleted = []
            start_pos = ref_pos_map[i]
            while i < len(ref_aln) and ref_aln[i] != '-' and qry_aln[i] == '-':
                deleted.append(ref_aln[i])
                i += 1
            end_pos = ref_pos_map[i-1] if i-1 >= 0 else start_pos
            ref_seq = ''.join(deleted)
            if len(deleted) == 1:
                notation = f"del{start_pos}"
            else:
                notation = f"del{start_pos}_{end_pos}"
            mutations.append({
                'type': 'deletion',
                'position': start_pos,
                'end_position': end_pos,
                'ref_nt': ref_seq,
                'alt_nt': '-',
                'notation': notation,
                'length': len(deleted)
            })
        else:
            i += 1
    return {
        'ref_sequence': ref_dna,
        'query_sequence': qry_dna,
        'ref_alignment': ref_aln,
        'query_alignment': qry_aln,
        'mutations': mutations,
        'snp_count': sum(1 for m in mutations if m['type'] == 'SNP'),
        'insertion_count': sum(1 for m in mutations if m['type'] == 'insertion'),
        'deletion_count': sum(1 for m in mutations if m['type'] == 'deletion')
    }


def categorize_nt_mutations(nt_mutations, ref_dna, query_dna, is_coding=True):
    snp_list = []
    ins_list = []
    del_list = []
    ref_dna = (ref_dna or "").upper()
    query_dna = (query_dna or "").upper()

    subst_at: Dict[int, str] = {}
    deleted_ref_positions = set()
    shift_events: List[Tuple[int, int]] = []

    for m in nt_mutations or []:
        mtype = m.get("type")
        if mtype == "SNP":
            if m.get("position") is not None:
                subst_at[int(m["position"])] = str(m.get("alt_nt") or "")
        elif mtype == "insertion":
            shift_events.append((int(m.get("position") or 1),
                                 int(m.get("length") or 0)))
        elif mtype == "deletion":
            start = int(m.get("position") or 1)
            end = int(m.get("end_position") or start)
            for p in range(start, end + 1):
                deleted_ref_positions.add(p)
            shift_events.append((start, -int(m.get("length") or (end - start + 1))))

    shift_events.sort(key=lambda t: t[0])

    def _net_shift_before(ref_position: int) -> int:
        total = 0
        for pos, delta in shift_events:
            if pos <= ref_position:
                total += delta
            else:
                break
        return total

    def _query_base_for(ref_position: int) -> Optional[str]:
        if ref_position in deleted_ref_positions:
            return None
        if ref_position in subst_at:
            return subst_at[ref_position]
        if 1 <= ref_position <= len(ref_dna):
            return ref_dna[ref_position - 1]
        return None

    missense_count = 0
    nonsense_count = 0
    syn_count = 0
    undetermined_count = 0

    for m in nt_mutations or []:
        if m["type"] == "SNP":
            pos = m.get("position")
            if not is_coding:
                m["effect"] = "SNP (noncoding)"
                m["effect_class"] = "noncoding"
                undetermined_count += 1
                snp_list.append(m)
                continue
            if pos is None:
                m["effect"] = "SNP"
                m["effect_class"] = "undetermined"
                undetermined_count += 1
                snp_list.append(m)
                continue

            pos = int(pos)
            codon_index = (pos - 1) // 3
            codon_start = codon_index * 3 + 1
            codon_positions = [codon_start, codon_start + 1, codon_start + 2]

            if codon_positions[-1] > len(ref_dna):
                m["effect"] = "SNP (incomplete reference codon)"
                m["effect_class"] = "undetermined"
                undetermined_count += 1
                snp_list.append(m)
                continue

            if _net_shift_before(codon_start) % 3 != 0:
                m["effect"] = "SNP (downstream of frameshift)"
                m["effect_class"] = "undetermined"
                undetermined_count += 1
                snp_list.append(m)
                continue

            qry_bases = [_query_base_for(p) for p in codon_positions]
            if any(b is None or len(b) != 1 for b in qry_bases):
                m["effect"] = "SNP (indel-affected codon)"
                m["effect_class"] = "undetermined"
                undetermined_count += 1
                snp_list.append(m)
                continue

            ref_codon = ref_dna[codon_start - 1: codon_start + 2]
            qry_codon = "".join(qry_bases)
            ref_aa = translate_dna_to_protein(ref_codon)
            qry_aa = translate_dna_to_protein(qry_codon)
            m["ref_codon"] = ref_codon
            m["query_codon"] = qry_codon
            m["codon_number"] = codon_index + 1
            m["ref_aa"] = ref_aa
            m["query_aa"] = qry_aa
            if ref_aa == qry_aa:
                m["effect"] = "Synonymous (silent)"
                m["effect_class"] = "synonymous"
                syn_count += 1
            elif qry_aa == "*":
                m["effect"] = "Nonsense"
                m["effect_class"] = "nonsense"
                nonsense_count += 1
            else:
                m["effect"] = "Non-synonymous (missense)"
                m["effect_class"] = "missense"
                missense_count += 1
            snp_list.append(m)
        elif m["type"] == "insertion":
            ins_list.append(m)
        elif m["type"] == "deletion":
            del_list.append(m)

    syn_snps = [m for m in snp_list if m.get("effect_class") == "synonymous"]
    non_syn_snps = [m for m in snp_list if m.get("effect_class") != "synonymous"]
    return {
        "syn_snps": syn_snps,
        "non_syn_snps": non_syn_snps,
        "ins_list": ins_list,
        "del_list": del_list,
        "syn_count": syn_count,
        "non_syn_count": len(non_syn_snps),
        "missense_count": missense_count,
        "nonsense_count": nonsense_count,
        "undetermined_count": undetermined_count,
        "snp_count": len(snp_list),
        "ins_count": len(ins_list),
        "del_count": len(del_list),
    }


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


def load_curated(path: str) -> dict:
    with open(path) as f:
        return json.load(f)
