"""Compatibility facade for the ASE-backed LAMMPS data exporter.

The public module name is retained so embedding code and saved projects do not
need to change. Geometry and atom creation now happen in :mod:`ase_writer`;
LAMMPS is not needed to generate the sample; it can optionally relax an
exported data file through the GUI workflow.
"""
from __future__ import annotations

from dataclasses import dataclass

from ase import Atoms
from ase_writer import BuildResult, WriterOptions, build_atoms, export_lammps, export_xyz
from ase_writer import build_type_table, estimate_atoms


@dataclass
class MDOptions:
    """LAMMPS workflow options retained for project-file compatibility."""

    potential_key: str = ""
    workflow: str = "relax"
    temperature: float = 300.0
    timestep_fs: float = 0.0
    ensemble: str = "nvt"
    equilibration_steps: int = 20000
    snapshots: int = 10
    snapshot_interval: int = 500
    seed: int = 90210
    backend: str = "lammps"
    mace_model: str = "medium"
    mace_fmax: float = 0.05


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
        result.warnings.append("MD options are applied when the exported LAMMPS input is run.")
    return result


__all__ = [
    "BuildResult", "MDOptions", "WriterOptions", "build_script",
    "export_lammps", "export_xyz",
]
