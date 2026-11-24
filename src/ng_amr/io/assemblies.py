"""Discovery of genome assembly FASTA files."""
from __future__ import annotations

from collections import defaultdict

from ng_amr.core.constants import (
    SEQUENCE_EXT,
)
from ng_amr.core.ids import (
    assembly_id,
)


def discover_assemblies(base, audit, mapping=None):
    result = defaultdict(dict)
    nt = base / "Assembly" / "NT"
    for kind, directory in (("annotated", nt), ("plane", nt/"Plane_Fasta")):
        groups = defaultdict(list)
        if directory.is_dir():
            for path in sorted(directory.iterdir()):
                if path.is_file() and path.suffix.lower() in SEQUENCE_EXT:
                    groups[(mapping or {}).get(path.name,assembly_id(path.stem))].append(path)
        for key, paths in groups.items():
            result[key][kind + "_files"] = paths
            if len(paths) == 1:
                result[key][kind] = str(paths[0])
            else:
                audit.add("assemblies", "collision", f"{kind} {key}: {paths}")
    return result
