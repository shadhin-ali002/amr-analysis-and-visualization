"""Fixed configuration: antibiotics, gene panel, gene types and output folder names."""
from __future__ import annotations

import re
from typing import Dict, List

QC = "ATCC49226"
SEQUENCE_EXT = {".fa", ".fasta", ".fna"}
ANTIBIOTIC_ORDER = [
    "Spectinomycin", "Penicillin G", "Ciprofloxacin", "Azithromycin",
    "Tetracycline", "Ceftriaxone", "Cefixime", "Cefepime",
]


ANTIBIOTIC_SHORT = {
    "Spectinomycin": "SPC",
    "Penicillin G": "PG",
    "Ciprofloxacin": "CIP",
    "Azithromycin": "AZM",
    "Tetracycline": "TET",
    "Ceftriaxone": "CRO",
    "Cefixime": "CFM",
    "Cefepime": "FEP",
}


ORIGINAL, ROUNDED, AD_MIC = MIC_FIELDS = ("E-strip original MIC", "E-strip rounded MIC", "AD MIC")
EST_CAT, AD_CAT, DD_CAT, WGS_CAT = CATEGORIES = (
    "E-strip interpretation", "AD interpretation", "Disk diffusion category", "Recorded WGS prediction")


EST_MCF, AD_MCF = "E-strip McFarland", "AD McFarland"
DATE_FIELDS = ("E-strip applied", "E-strip read", "AD inoculated", "AD read", "E-strip McFarland date", "AD McFarland date")
MEASURE = re.compile(r"^\s*([<>≤≥]=?|=)?\s*(\d+(?:\.\d*)?|\.\d+)(?:[eE]([+-]?\d+))?\s*$")
CODING = {"penA", "ponA1", "porB1b", "gyrA", "parC", "parE", "mtrR", "mtrC", "folP", "rpoB", "rpoD",
              "rplD", "rplV", "rpsE", "rpsJ", "blaTEM", "tetM", "ermA", "ermB", "ermC", "ermF", "ereA", "ereB", "mefA"}


ORIGINAL_AA_FOLDER = "Assembly/AA"
AA_COMPARISON_FOLDER_NAME = "AA_Comparison"
EXTRACTED_SEQUENCES_FOLDER_NAME = "Extracted_Sequences"
GENE_AA_PRODUCT_SYNONYMS: Dict[str, List[str]] = {
    "penA": [
        "penicillin-binding protein 2",
        "penicillin binding protein 2",
        "peptidoglycan d,d-transpeptidase pbp2",
        "pbp2",
        "pbp-2",
    ],
    "ponA1": [
        "penicillin-binding protein 1a",
        "penicillin binding protein 1a",
        "penicillin-binding protein 1",
        "pbp1a",
        "pbp-1a",
        "pona",
    ],
    "porB1b": [
        "trimeric porin porb.ib",
        "trimeric porin porb",
        "porb.ib",
        "porb1b",
        "porin porb",
        "outer membrane protein porb",
    ],
    "blaTEM": [
        "tem family class a beta-lactamase",
        "tem family",
        "class a beta-lactamase",
        "beta-lactamase tem",
        "blatem",
        "tem-1",
        "tem1",
    ],
    "parC": [
        "dna topoisomerase iv subunit a",
        "topoisomerase iv subunit a",
        "parc",
    ],
    "parE": [
        "dna topoisomerase iv subunit b",
        "topoisomerase iv subunit b",
        "pare",
    ],
    "gyrA": [
        "dna gyrase subunit a",
        "gyrase subunit a",
        "gyra",
    ],
    "mtrR": [
        "multidrug efflux system transcriptional repressor mtrr",
        "transcriptional repressor mtrr",
        "mtrr",
    ],
    "mtrC": [
        "multidrug efflux rnd transporter periplasmic adaptor subunit mtrc",
        "periplasmic adaptor subunit mtrc",
        "mtrc",
    ],
    "mtrD": [
        "multidrug efflux rnd transporter permease subunit mtrd",
        "permease subunit mtrd",
        "mtrd",
    ],
    "rpsJ": [
        "30s ribosomal protein s10",
        "ribosomal protein s10",
        "rpsj",
    ],
    "folP": [
        "dihydropteroate synthase",
        "dihydropteroate",
        "folp",
    ],
    "rpsE": [
        "30s ribosomal protein s5",
        "ribosomal protein s5",
        "rpse",
    ],
    "rpoB": [
        "dna-directed rna polymerase subunit beta",
        "rna polymerase subunit beta",
        "rpob",
    ],
    "rpoD": [
        "rna polymerase sigma factor rpod",
        "sigma factor rpod",
        "rpod",
    ],
    "rplD": [
        "50s ribosomal protein l4",
        "ribosomal protein l4",
        "rpld",
    ],
    "rplV": [
        "50s ribosomal protein l22",
        "ribosomal protein l22",
        "rplv",
    ],
    "macAB_promoter": [],
    "mtrR_promoter": [],
    "16S_rDNA": [],
    "23S_rDNA": [],
    "ermA": ["ermA", "rRNA methyltransferase ermA", "23s rrna methyltransferase erma"],
    "ermB": ["ermB", "rRNA methyltransferase ermB", "23s rrna methyltransferase ermb"],
    "ermC": ["ermC", "rRNA methyltransferase ermC", "23s rrna methyltransferase ermc"],
    "ermF": ["ermF", "rRNA methyltransferase ermF", "23s rrna methyltransferase ermf"],
    "ereA": ["ereA", "erythromycin esterase", "erythromycin esterase erea"],
    "ereB": ["ereB", "erythromycin esterase ereb"],
    "mefA": ["mefA", "macrolide efflux", "mef(a)", "macrolide-efflux"],
}


DEFAULT_GENES = [
    "penA", "ponA1", "porB1b", "blaTEM", "parC", "parE", "gyrA", "mtrR", "mtrR_promoter", "mtrC", "mtrD",
    "macAB_promoter", "rpsJ", "folP", "rpsE", "16S_rDNA", "rpoB", "rpoD", "rplD", "rplV", "23S_rDNA",
    "ermA", "ermB", "ermC", "ermF", "ereA", "ereB", "mefA",
    "norM_promoter", "tetM",
]


GENE_FOLDER_NAMES = {
    "penA":           "1. penA",
    "ponA1":          "2. ponA1",
    "porB1b":         "3. porB1b",
    "blaTEM":         "4. blaTEM",
    "parC":           "5. parC",
    "parE":           "6. parE",
    "gyrA":           "7. gyrA",
    "mtrR":           "8. mtrR",
    "mtrR_promoter":  "9. mtrR_promoter",
    "mtrC":           "10. mtrC",
    "mtrD":           "11. mtrD",
    "macAB_promoter": "12. macAB_promoter",
    "rpsJ":           "13. rpsJ",
    "folP":           "14. folP",
    "rpsE":           "15. rpsE",
    "16S_rDNA":       "16. 16S_rDNA",
    "rpoB":           "17. rpoB",
    "rpoD":           "18. rpoD",
    "rplD":           "19. rplD",
    "rplV":           "20. rplV",
    "23S_rDNA":       "21. 23S_rDNA",
    "ermA":           "22. ermA",
    "ermB":           "23. ermB",
    "ermC":           "24. ermC",
    "ermF":           "25. ermF",
    "ereA":           "26. ereA",
    "ereB":           "27. ereB",
    "mefA":           "28. mefA",
    "norM_promoter":  "29. norM_promoter",
    "tetM":           "30. tetM",
}


GENE_CONFIG = {
    "penA":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ponA1":           {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "porB1b":          {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "blaTEM":          {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "gyrA":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "parC":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "parE":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "mtrR":            {"gene_type": "CDS",       "multi_ref": False, "disruption": True},
    "mtrR_promoter":   {"gene_type": "noncoding", "multi_ref": True,  "disruption": True},
    "mtrC":            {"gene_type": "CDS",       "multi_ref": False, "disruption": True},
    "mtrD":            {"gene_type": "CDS",       "multi_ref": True,  "disruption": True},
    "macAB_promoter":  {"gene_type": "noncoding", "multi_ref": True,  "disruption": True},
    "norM_promoter":   {"gene_type": "noncoding", "multi_ref": False, "disruption": False},
    "23S_rDNA":        {"gene_type": "rRNA",      "multi_ref": False, "disruption": False},
    "16S_rDNA":        {"gene_type": "rRNA",      "multi_ref": False, "disruption": False},
    "tetM":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rpsJ":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rpoB":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "folP":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rpsE":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rpoD":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rplD":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "rplV":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ermA":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ermB":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ermC":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ermF":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ereA":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "ereB":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
    "mefA":            {"gene_type": "CDS",       "multi_ref": False, "disruption": False},
}


DISRUPTION_MIN_DELETION_AA = 15
DISRUPTION_MIN_DELETION_EVENTS = 3
DISRUPTION_MAX_LENGTH_CHANGE = 0.20


GENE_SHEET_ORDER = [
    "penA", "ponA1", "porB1b", "blaTEM",
    "parC", "parE", "gyrA",
    "mtrR", "mtrR_promoter", "mtrC", "mtrD",
    "macAB_promoter",
    "rpsJ", "folP", "rpsE", "16S_rDNA", "rpoB", "rpoD",
    "rplD", "rplV", "23S_rDNA",
    "ermA", "ermB", "ermC", "ermF", "ereA", "ereB", "mefA",
    "norM_promoter", "tetM",
]
