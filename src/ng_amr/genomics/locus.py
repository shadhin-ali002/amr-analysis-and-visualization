"""BLAST+ search and per-locus variant calling against the reference panel."""
from __future__ import annotations

import atexit
import csv
import os
import re
import shutil
import statistics as _stats
import subprocess
import tempfile as _tempfile
from pathlib import Path
from typing import Dict, List, Optional

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
    logger,
)
from ng_amr.core.constants import (
    CODING,
    QC,
)
from ng_amr.genomics.alignment import (
    safe_name,
)

_CACHE = {"hits": {}, "contigs": {}, "refs": None}
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
    return {"gene": gene_name, "processed": processed, "skipped": skipped, "mutation_counts": mutation_counts,
            "total_isolates": total_isolates, "csv": csv_out, "matrix_csv": matrix_csv, "excel": str(Path(output_root) / "Mutation_Reports" / "By_Gene" / (safe_name(gene_name) + ".xlsx")),
            "all_mutation_data": all_mutation_data, "all_data": all_data, "curated_set": curated_set,
            "n_found": processed, "n_collection": total_isolates}
