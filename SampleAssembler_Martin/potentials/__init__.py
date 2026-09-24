"""Curated interatomic-potential library for Sample Assembler MD.

See ``registry.py`` for the loader and ``../MD_PLAYBOOK.md`` for the design and
per-material potential choices. Each potential is a parameter file (or files)
plus a JSON sidecar carrying its citation, element set, required units /
atom_style, and validity notes.
"""
from .registry import POTENTIALS, Potential, load_potentials  # noqa: F401
