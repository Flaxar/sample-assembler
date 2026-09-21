import re
from pathlib import Path

from build123d import import_step


def _step_component_names(filepath):
    """Read assembly occurrence names because build123d does not expose STEP labels."""
    text = Path(filepath).read_text(encoding="latin-1")
    names = re.findall(
        r"NEXT_ASSEMBLY_USAGE_OCCURRENCE\('([^:']+):[^']*'",
        text,
    )
    return list(dict.fromkeys(names))


def __import_part_shape(file_path):
    """Deprecated. Use import_multi_body_step() instead."""
    shapes = import_multi_body_step(file_path)
    return next(iter(shapes.values()), None)


def import_multi_body_step(filepath):
    """
    Import a STEP file and return a mapping of component names to solid shapes.

    Assembly occurrence names are recovered from STEP metadata where available.
    Unnamed solids receive deterministic ``Shape N`` names.
    """
    imported = import_step(filepath)
    solids = list(imported)
    names = _step_component_names(filepath)

    shapes_dict = {}
    for index, shape in enumerate(solids):
        name = names[index] if index < len(names) else f"Shape {index + 1}"
        if name in shapes_dict:
            name = f"{name} ({index + 1})"
        shapes_dict[name] = shape
        print(f"Found part named: {name}")

    return shapes_dict
