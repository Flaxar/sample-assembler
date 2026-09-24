"""Compatibility facade for the ASE-backed LAMMPS data exporter.

The public module name is retained so embedding code and saved projects do not
need to change. Geometry and atom creation now happen in :mod:`ase_writer`;
LAMMPS is no longer needed to generate the sample.
"""
from __future__ import annotations

from dataclasses import dataclass

from ase import Atoms
from ase_writer import BuildResult, WriterOptions, build_atoms, export_lammps, export_xyz
from ase_writer import build_type_table, estimate_atoms


@dataclass
class MDOptions:
    """Retained for project-file compatibility; MD is no longer run here."""

    potential_key: str = ""
    workflow: str = "thermalize"


def build_script(fills, opts: WriterOptions, title="Sample Assembler build", md=None):
    """Return a lightweight preview; coordinates are built only on export."""
    enabled = [fill for fill in fills if fill.enabled]
    result = BuildResult(
        atoms=Atoms(),
        type_table=build_type_table(enabled),
        est_atoms=sum(estimate_atoms(fill) for fill in enabled),
        per_fill_atoms=[(fill.name, estimate_atoms(fill)) for fill in enabled],
    )
    if md is not None:
        result.warnings.append("MD options are ignored: ASE export does not run LAMMPS.")
    return result


__all__ = [
    "BuildResult", "MDOptions", "WriterOptions", "build_script",
    "export_lammps", "export_xyz",
]
