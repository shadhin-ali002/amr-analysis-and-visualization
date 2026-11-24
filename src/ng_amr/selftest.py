"""Built-in regression checks (`ng-amr --self-test`)."""
from __future__ import annotations

import math
import os
import re
from typing import Optional

try:
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
    from Bio.Seq import Seq
except ImportError:
    SeqIO = Seq = PairwiseAligner = None

from ng_amr.core.audit import (
    logger,
)
from ng_amr.core.ids import (
    normalize_sample_id,
    short_isolate_label,
)
from ng_amr.genomics.alignment import (
    analyze_nucleotide,
    analyze_protein,
    assess_disruption,
    call_protein_mutations,
    categorize_nt_mutations,
    get_gene_type,
)


def run_self_tests() -> Optional[bool]:
    failures = []
    checks_run = [0]

    def check(label, condition):
        checks_run[0] += 1
        if not condition:
            failures.append(label)
            logger.error("SELF-TEST FAILED: %s", label)
        else:
            logger.info("SELF-TEST OK: %s", label)

    try:
        check("short_isolate_label truncates long names", len(short_isolate_label("A" * 50, max_len=18)) <= 18)
        check("normalize_sample_id is deterministic", normalize_sample_id("Sample_1") == normalize_sample_id("Sample_1"))
        check("normalize_sample_id lower-cases output", normalize_sample_id("SAMPLE_ABC") == normalize_sample_id("SAMPLE_ABC").lower())
    except Exception as exc:
        failures.append(f"ID normalisation raised {exc!r}")
        logger.error("SELF-TEST FAILED: ID normalisation raised %r", exc)

    try:
        muts = call_protein_mutations("MSTV", "MSTV")
        check("call_protein_mutations finds no changes for identical sequences", muts == [])

        muts = call_protein_mutations("MSTV", "MSTF")
        check("call_protein_mutations detects a simple substitution",
              any(m["type"] == "substitution" and m["notation"] == "V4F" for m in muts))

        muts = call_protein_mutations("MST-V", "MSTAV")
        check("call_protein_mutations detects an insertion", any(m["type"] == "insertion" for m in muts))

        muts = call_protein_mutations("MSTAV", "MST-V")
        check("call_protein_mutations detects a deletion", any(m["type"] == "deletion" for m in muts))

        muts = call_protein_mutations("MSTV", "MST*")
        check("call_protein_mutations detects a stop-gain", any(m["type"] == "stop_gain" for m in muts))
    except Exception as exc:
        failures.append(f"call_protein_mutations raised {exc!r}")
        logger.error("SELF-TEST FAILED: call_protein_mutations raised %r", exc)

    check("16S_rDNA and 23S_rDNA are configured as non-CDS",
          get_gene_type("16S_rDNA") != "CDS" and get_gene_type("23S_rDNA") != "CDS")
    check("promoter regions are configured as non-CDS",
          get_gene_type("mtrR_promoter") != "CDS"
          and get_gene_type("macAB_promoter") != "CDS")

    if Seq is not None and PairwiseAligner is not None:
        try:
            ref = "ATGAAACCCGGGTTTAAACCCGGGTAA"
            qry_fs = "ATGAACCCGGGTTTAAACCCGGGTAA"
            nt_fs = analyze_nucleotide(ref, qry_fs)
            prot_fs = analyze_protein("selftest", ref, qry_fs,
                                      nt_mutations=nt_fs["mutations"])
            check("1 bp deletion is called as a frameshift",
                  len(prot_fs["frameshifts"]) == 1
                  and prot_fs["frameshifts"][0]["type"] == "frameshift_deletion")
            check("frameshifted gene is assessed as DISRUPTED",
                  assess_disruption(prot_fs) is True)

            qry_inframe = "ATG" + "CCCGGGTTTAAACCCGGGTAA"
            nt_if = analyze_nucleotide(ref, qry_inframe)
            prot_if = analyze_protein("selftest", ref, qry_inframe,
                                     nt_mutations=nt_if["mutations"])
            check("in-frame 3 bp deletion is NOT called a frameshift",
                  prot_if["frameshifts"] == [])
        except Exception as exc:
            failures.append(f"frameshift regression raised {exc!r}")
            logger.error("SELF-TEST FAILED: frameshift regression raised %r", exc)

        try:
            ref_syn = "ATGAAACCCGGGTTTTAA"
            qry_syn = "ATGAAACCCGGGTTCTAA"
            nt_syn = analyze_nucleotide(ref_syn, qry_syn)
            cat_syn = categorize_nt_mutations(nt_syn["mutations"], ref_syn, qry_syn, True)
            check("synonymous SNP is called synonymous (no indel present)",
                  cat_syn["syn_count"] == 1 and cat_syn["missense_count"] == 0)

            ref_i = "ATGAAACCCGGGTTTTAA"
            qry_i = "ATGAAAGCCCGGGTTCTAA"
            nt_i = analyze_nucleotide(ref_i, qry_i)
            cat_i = categorize_nt_mutations(nt_i["mutations"], ref_i, qry_i, True)
            snp_effects = [m.get("effect") for m in nt_i["mutations"] if m["type"] == "SNP"]
            check("SNP downstream of an indel is not mis-called missense",
                  all(e != "Non-synonymous (missense)" for e in snp_effects))
            check("frame-affected SNPs are labelled, not silently guessed",
                  all((e or "").startswith("SNP (") for e in snp_effects))
            check("categorize_nt_mutations keeps syn + non_syn == snp_count",
                  cat_i["syn_count"] + cat_i["non_syn_count"] == cat_i["snp_count"])
        except Exception as exc:
            failures.append(f"syn/non-syn regression raised {exc!r}")
            logger.error("SELF-TEST FAILED: syn/non-syn regression raised %r", exc)

        try:
            ref_r = "AAACGTACGTAAACGT"
            qry_r = "AAACGTACGTGAACGT"
            nt_r = analyze_nucleotide(ref_r, qry_r)
            notations = [m["notation"] for m in nt_r["mutations"]]
            check("non-coding variants use nucleotide notation (e.g. A11G)",
                  any(re.fullmatch(r"[ACGT]\d+[ACGT]", n) for n in notations))
        except Exception as exc:
            failures.append(f"non-coding regression raised {exc!r}")
            logger.error("SELF-TEST FAILED: non-coding regression raised %r", exc)

        try:
            prot_stop = analyze_protein("selftest", "ATGAAACCCTAA", "ATGTAACCCTAA")
            check("stop-gain is present in the reported mutations",
                  any(m["type"] == "stop_gain" for m in prot_stop["mutations"]))
            check("silent_mutations stays empty (documented invariant)",
                  prot_stop["silent_mutations"] == [])
        except Exception as exc:
            failures.append(f"stop-gain regression raised {exc!r}")
            logger.error("SELF-TEST FAILED: stop-gain regression raised %r", exc)

    try:
        check(">=35% of 24 isolates requires 9 isolates, not 8",
              max(1, int(math.ceil(24 * 0.35))) == 9)
        check(">=25% of 24 isolates requires 6 isolates",
              max(1, int(math.ceil(24 * 0.25))) == 6)
    except Exception as exc:
        failures.append(f"threshold regression raised {exc!r}")
        logger.error("SELF-TEST FAILED: threshold regression raised %r", exc)

    try:
        import ast as _ast
        _lvl = {"debug", "info", "warning", "error", "exception", "critical"}
        _ph = re.compile(r"%(?:\([^)]*\))?[-#0 +]*\d*(?:\.\d+)*[hlL]?[diouxXeEfFgGcrsa%]")

        def _literal_text(node):
            if isinstance(node, _ast.Constant) and isinstance(node.value, str):
                return node.value
            if isinstance(node, _ast.JoinedStr):
                return "".join(v.value for v in node.values
                               if isinstance(v, _ast.Constant)
                               and isinstance(v.value, str))
            if isinstance(node, _ast.BinOp) and isinstance(node.op, _ast.Add):
                left, right = _literal_text(node.left), _literal_text(node.right)
                if left is None or right is None:
                    return None
                return left + right
            return None

        bad_calls = []
        _pkg = os.path.dirname(os.path.abspath(__file__))
        _nodes = []
        for _dir, _subdirs, _files in os.walk(_pkg):
            for _name in sorted(_files):
                if _name.endswith(".py"):
                    _path = os.path.join(_dir, _name)
                    with open(_path, encoding="utf-8") as _fh:
                        _rel = os.path.relpath(_path, _pkg)
                        _nodes += [(_rel, _n) for _n in _ast.walk(_ast.parse(_fh.read()))]
        for _file, _n in _nodes:
            if not isinstance(_n, _ast.Call) or len(_n.args) < 2:
                continue
            _f = _n.func
            if not (isinstance(_f, _ast.Attribute) and _f.attr in _lvl):
                continue
            if not (isinstance(_f.value, _ast.Name) and _f.value.id == "logger"):
                continue
            _txt = _literal_text(_n.args[0])
            if _txt is None:
                continue
            _found = [x for x in _ph.findall(_txt) if x != "%%"]
            if len(_found) != len(_n.args) - 1:
                bad_calls.append(f"{_file}:{_n.lineno}")
        check("every logger call's args match its %-placeholders",
              not bad_calls)
        if bad_calls:
            logger.error("  offending logger call(s) at line(s): %s", bad_calls)
    except Exception as exc:
        failures.append(f"logger-format audit raised {exc!r}")
        logger.error("SELF-TEST FAILED: logger-format audit raised %r", exc)

    logger.info("-" * 60)
    if failures:
        logger.error("SELF-TEST SUITE: %d of %d check(s) failed: %s",
                     len(failures), checks_run[0], failures)
        return False
    if Seq is None or PairwiseAligner is None:
        logger.warning("%d checks passed; 10 alignment-dependent checks unavailable (Biopython missing).", checks_run[0])
        return None
    logger.info("ALL SELF-TESTS PASSED (%d checks)", checks_run[0])
    return True
