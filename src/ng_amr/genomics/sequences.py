"""Export of extracted gene and protein sequences."""
from __future__ import annotations

import csv
import os
import re
from typing import Dict, List, Optional

import pandas as pd

from ng_amr.core.audit import (
    logger,
)
from ng_amr.core.constants import (
    DEFAULT_GENES,
    EXTRACTED_SEQUENCES_FOLDER_NAME,
    GENE_CONFIG,
)


def _safe_isolate_dirname(name: str, max_len: int = 80) -> str:
    s = re.sub(r'[<>:"/\\|?*\x00-\x1f]+', "_", str(name or "isolate"))
    s = re.sub(r"_+", "_", s).strip(" ._")
    return (s[:max_len] or "isolate")


def _write_fasta_record(path: str, header: str, sequence: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    seq = re.sub(r"\s+", "", str(sequence or ""))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f">{header}\n")
        if seq:
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")
        else:
            fh.write("\n")


def generate_extracted_sequences_package(
    gene_results: List[dict],
    output_root: str,
) -> Optional[str]:
    out_dir = os.path.join(output_root, EXTRACTED_SEQUENCES_FOLDER_NAME)
    os.makedirs(out_dir, exist_ok=True)
    logger.info("  [Extracted_Sequences] Writing per-isolate NT/AA under %s", out_dir)

    genes_order: List[str] = list(DEFAULT_GENES)
    by_iso: Dict[str, Dict[str, dict]] = {}
    for gr in gene_results or []:
        gene = str(gr.get("gene") or "gene")
        if gene not in genes_order:
            genes_order.append(gene)
        for d in gr.get("all_data") or []:
            iso = str(d.get("isolate") or "").strip()
            if not iso:
                continue
            by_iso.setdefault(iso, {})[gene] = {
                "nt": d.get("query_dna") or "",
                "aa": d.get("query_protein") or "",
                "status": d.get("status") or "",
                "gene_found": bool(d.get("gene_found")),
                "ref_name": d.get("ref_name") or gene,
            }

    for iso in list(by_iso.keys()):
        for gene in genes_order:
            by_iso[iso].setdefault(gene, {
                "nt": "", "aa": "", "status": "MISSING",
                "gene_found": False, "ref_name": gene,
            })

    n_iso = 0
    n_nt_ok = 0
    n_aa_ok = 0
    n_nt_files = 0
    n_aa_files = 0

    def _sort_iso(name: str):
        m = re.match(r"^(\d+)", name.strip())
        return (int(m.group(1)) if m else 10**9, name.lower())

    for iso, gene_map in sorted(by_iso.items(), key=lambda x: _sort_iso(x[0])):
        iso_dir = os.path.join(out_dir, _safe_isolate_dirname(iso))
        nt_dir = os.path.join(iso_dir, "NT")
        aa_dir = os.path.join(iso_dir, "AA")
        os.makedirs(nt_dir, exist_ok=True)
        os.makedirs(aa_dir, exist_ok=True)
        n_iso += 1

        manifest_rows = []
        for gene in genes_order:
            info = gene_map.get(gene) or {}
            nt_seq = info.get("nt") or ""
            aa_seq = info.get("aa") or ""
            found = bool(info.get("gene_found"))
            safe_gene = re.sub(r"[^A-Za-z0-9._-]+", "_", gene)
            gene_type = GENE_CONFIG.get(gene, {}).get("gene_type", "CDS")

            nt_path = os.path.join(nt_dir, f"{safe_gene}.fasta")
            if found and nt_seq:
                _write_fasta_record(nt_path, f"{iso}|{gene}|NT|len={len(nt_seq)}", nt_seq)
                n_nt_ok += 1
                nt_status = "OK"
            else:
                _write_fasta_record(
                    nt_path,
                    f"{iso}|{gene}|NT|MISSING|gene not found in assembly",
                    "",
                )
                nt_status = "MISSING"
            n_nt_files += 1

            aa_path = os.path.join(aa_dir, f"{safe_gene}.faa")
            if found and aa_seq:
                _write_fasta_record(aa_path, f"{iso}|{gene}|AA|len={len(aa_seq)}", aa_seq)
                n_aa_ok += 1
                aa_status = "OK"
            elif gene_type != "CDS":
                _write_fasta_record(
                    aa_path,
                    f"{iso}|{gene}|AA|N/A_noncoding|no protein product expected",
                    "",
                )
                aa_status = "N/A_noncoding"
            else:
                _write_fasta_record(
                    aa_path,
                    f"{iso}|{gene}|AA|MISSING|no translated AA available",
                    "",
                )
                aa_status = "MISSING"
            n_aa_files += 1

            manifest_rows.append({
                "Gene": gene,
                "Gene_type": gene_type,
                "NT_status": nt_status,
                "NT_length": len(nt_seq) if nt_seq else 0,
                "AA_status": aa_status,
                "AA_length": len(aa_seq) if aa_seq else 0,
            })

        safe_iso = _safe_isolate_dirname(iso)
        comb_nt = os.path.join(iso_dir, f"{safe_iso}_all_genes_NT.fasta")
        comb_aa = os.path.join(iso_dir, f"{safe_iso}_all_genes_AA.faa")
        with open(comb_nt, "w", encoding="utf-8") as fh_nt, \
             open(comb_aa, "w", encoding="utf-8") as fh_aa:
            for gene in genes_order:
                info = gene_map.get(gene) or {}
                nt_seq = re.sub(r"\s+", "", info.get("nt") or "")
                aa_seq = re.sub(r"\s+", "", info.get("aa") or "")
                if nt_seq:
                    fh_nt.write(f">{iso}|{gene}|NT|len={len(nt_seq)}\n")
                    for i in range(0, len(nt_seq), 70):
                        fh_nt.write(nt_seq[i:i + 70] + "\n")
                if aa_seq:
                    fh_aa.write(f">{iso}|{gene}|AA|len={len(aa_seq)}\n")
                    for i in range(0, len(aa_seq), 70):
                        fh_aa.write(aa_seq[i:i + 70] + "\n")

        man_path = os.path.join(iso_dir, "extraction_manifest.csv")
        try:
            if pd is not None:
                pd.DataFrame(manifest_rows).to_csv(man_path, index=False)
            else:
                with open(man_path, "w", newline="", encoding="utf-8") as fh:
                    w = csv.DictWriter(
                        fh,
                        fieldnames=["Gene", "Gene_type", "NT_status", "NT_length",
                                    "AA_status", "AA_length"],
                    )
                    w.writeheader()
                    w.writerows(manifest_rows)
        except Exception as e:
            logger.warning("  [Extracted_Sequences] manifest failed for %s: %s", iso, e)

    logger.info(
        "  [Extracted_Sequences] Done — %d isolate folder(s), "
        "%d NT files (%d OK), %d AA files (%d OK) → %s",
        n_iso, n_nt_files, n_nt_ok, n_aa_files, n_aa_ok, out_dir,
    )
    return out_dir
