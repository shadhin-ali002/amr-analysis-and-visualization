"""Generate the synthetic genomic demo dataset used by the BLAST+ test.

Every sequence is random (fixed seed). Nothing here is real N. gonorrhoeae sequence; the
gene names are only labels the pipeline recognises. Re-running this script reproduces the
committed files byte for byte.

    python tests/fixtures/make_genomic_demo.py
"""
from __future__ import annotations

import json
import random
from pathlib import Path

OUT = Path(__file__).resolve().parent / "genomic_demo"
STOPS = {"TAA", "TAG", "TGA"}
rng = random.Random(49226)


def dna(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def cds(codons: int) -> str:
    body = []
    while len(body) < codons:
        c = dna(3)
        if c not in STOPS:
            body.append(c)
    return "ATG" + "".join(body) + "TAA"


def set_codon(seq: str, aa_pos: int, codon: str) -> str:
    """Replace codon number aa_pos (1-based, ATG = 1)."""
    i = (aa_pos - 1) * 3
    return seq[:i] + codon + seq[i + 3:]


def flip(seq: str, pos: int) -> str:
    """Change the base at 1-based position pos."""
    base = seq[pos - 1]
    return seq[: pos - 1] + {"A": "G", "G": "A", "C": "T", "T": "C"}[base] + seq[pos:]


def wrap(seq: str) -> str:
    return "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70))


REFS = {
    "penA": cds(200),
    "gyrA": cds(150),
    "norM_promoter": dna(400),
    "23S_rDNA": dna(900),
}
# Make known codons so the demo mutations have readable, fixed names.
REFS["penA"] = set_codon(REFS["penA"], 10, "GCT")   # Ala10
REFS["penA"] = set_codon(REFS["penA"], 120, "GGT")  # Gly120
REFS["gyrA"] = set_codon(REFS["gyrA"], 91, "TCT")   # Ser91

BACKGROUND = [dna(6000) for _ in range(4)]


def genome(genes: dict[str, str], rrna_copies: list[str]) -> str:
    """Four contigs: genes embedded in fixed random background, rRNA in two copies."""
    contigs = []
    names = list(genes)
    for i, bg in enumerate(BACKGROUND):
        seq = bg
        if i < len(names):
            seq = bg[:2500] + genes[names[i]] + bg[2500:]
        if i == 3:
            seq = bg[:1500] + rrna_copies[0] + bg[1500:3500] + rrna_copies[1] + bg[3500:]
        contigs.append(f">NODE_{i + 1}_length_{len(seq)}_cov_50.0\n{wrap(seq)}")
    return "\n".join(contigs) + "\n"


def isolate(name: str) -> str:
    g = {k: v for k, v in REFS.items() if k != "23S_rDNA"}
    r23 = [REFS["23S_rDNA"], REFS["23S_rDNA"]]
    if name == "DEMO-B":
        g["penA"] = set_codon(g["penA"], 10, "GTT")        # A10V
        g["gyrA"] = set_codon(g["gyrA"], 91, "TTT")        # S91F
        g["norM_promoter"] = flip(g["norM_promoter"], 300)
        r23 = [flip(c, 450) for c in r23]                  # same SNP in both copies
    if name == "DEMO-C":
        g["penA"] = set_codon(g["penA"], 120, "TGA")       # premature stop
        del g["gyrA"]                                      # locus absent
    return genome(g, r23)


def main() -> None:
    refs = OUT / "Refs"
    plane = OUT / "Assembly" / "NT" / "Plane_Fasta"
    refs.mkdir(parents=True, exist_ok=True)
    plane.mkdir(parents=True, exist_ok=True)
    (refs / "Genes_RefSeq.fasta").write_text(
        "".join(f">{k}\n{wrap(v)}\n" for k, v in REFS.items()), encoding="utf-8", newline="\n")
    (refs / "Curated_Variants.json").write_text(
        json.dumps({"gyrA": ["S91F"], "penA": ["A10V"]}, indent=2) + "\n", encoding="utf-8", newline="\n")
    for name in ("DEMO-A", "DEMO-B", "DEMO-C", "ATCC49226"):
        (plane / f"{name}.fasta").write_text(isolate(name), encoding="utf-8", newline="\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
