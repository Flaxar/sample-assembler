"""STEP import helpers for the Sample Assembler UI."""
from __future__ import annotations

import re
from pathlib import Path

from build123d import import_step


def _component_names(filepath: str) -> list[str]:
    """Read assembly occurrence names because build123d does not expose labels."""
    text = Path(filepath).read_text(encoding="latin-1")
    names = re.findall(
        r"NEXT_ASSEMBLY_USAGE_OCCURRENCE\('([^:']+):[^']*'",
        text,
    )
    return list(dict.fromkeys(names))


def import_step_components(filepath: str) -> dict[str, object]:
    """Import all solids in a STEP file, keyed by their component names."""
    imported = import_step(filepath)
    solids = list(imported)
    names = _component_names(filepath)

    components = {}
    for index, shape in enumerate(solids):
        name = names[index] if index < len(names) else f"Shape {index + 1}"
        if name in components:
            name = f"{name} ({index + 1})"
        components[name] = shape
    return components
