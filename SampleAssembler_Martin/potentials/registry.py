#!/usr/bin/env python3
"""
Curated interatomic-potential library for Sample Assembler MD.

Parallel to ``materials.py``: where that module is a registry of crystallography
(element-keyed unit cells), this is a registry of **interatomic potentials** —
each one a LAMMPS parameter file (or files) plus a JSON sidecar carrying its
citation, element set, required ``units``/``atom_style``, recommended timestep
and validity notes. ``lammps_writer`` consults it to emit a *real* MD potential
(``pair_style``/``pair_coeff``) instead of the clean-build ``pair_style zero``.

Layout::

    potentials/
      registry.py                  (this file)
      <family>/<name>.<ext>        parameter file(s)
      <family>/<name>.json         sidecar metadata (one per potential)

A sidecar is discovered automatically (recursive ``*.json`` scan), so adding a
potential = drop the file + its ``.json`` next to it. See ``../MD_PLAYBOOK.md``
for the sidecar schema, the per-material choices and the provenance rules.

THE ONE SUBTLE CORRECTNESS POINT: for the many-body styles (``tersoff``, ``sw``,
``vashishta``, ``eam/alloy``, ``meam``, ``snap``) the trailing element list on the
``pair_coeff * * <file> <e1> <e2> ...`` line maps **positionally** to LAMMPS atom
types 1..N. ``lammps_writer`` owns that global ordering (``build_type_table``,
sorted by mass). ``Potential.pair_coeff_line(type_table)`` therefore walks the
type table in id order and emits, per type, either the covered element symbol or
``NULL`` if this potential does not describe that type.
"""
from __future__ import annotations

import glob
import json
import os
from dataclasses import dataclass, field
from typing import List, Optional

_HERE = os.path.dirname(os.path.abspath(__file__))


def _lmp_path(p: str) -> str:
    """Absolute path with forward slashes (LAMMPS accepts these on Windows)."""
    return os.path.abspath(p).replace("\\", "/")


@dataclass
class Potential:
    key: str
    name: str
    elements: List[str]                       # elements this potential describes
    pair_style: str                           # full LAMMPS pair_style argument string
    param_files: List[str]                    # absolute paths (resolved at load)
    pair_coeff_form: str = "manybody"         # 'manybody' = `* * <files> <elems>`
    required_units: str = "metal"
    required_atom_style: str = "atomic"
    requires_qeq: bool = False
    qeq_fix: Optional[str] = None             # e.g. 'fix qeq all qeq/reaxff 1 0.0 10.0 1e-6 reaxff'
    needs_packages: List[str] = field(default_factory=list)
    recommended_timestep_fs: float = 1.0
    default_temperature_K: float = 300.0
    applicable_materials: List[str] = field(default_factory=list)
    validity_notes: str = ""
    confidence: str = "med"
    citation: str = ""
    source_url: str = ""
    source_note: str = ""
    license: str = ""
    sidecar_path: str = ""

    # ── coverage helpers ────────────────────────────────────────────────────────
    def covers(self, element: str) -> bool:
        return element in self.elements

    def covers_all(self, elements) -> bool:
        return all(e in self.elements for e in elements)

    def uncovered(self, elements) -> List[str]:
        return [e for e in elements if e not in self.elements]

    # ── LAMMPS line rendering ───────────────────────────────────────────────────
    def pair_style_line(self) -> str:
        return f"pair_style {self.pair_style}"

    def pair_coeff_line(self, type_table: dict) -> str:
        """`pair_coeff * * <files> <elem-in-type-order>` with NULL masking.

        ``type_table`` is ``{element: type_id}`` (1-indexed). Elements are emitted
        in ascending type-id order; a type not covered by this potential becomes
        ``NULL`` so LAMMPS knows that type carries no force from this style.
        """
        if self.pair_coeff_form != "manybody":
            raise NotImplementedError(
                f"pair_coeff_form '{self.pair_coeff_form}' not yet supported "
                f"(potential '{self.key}'). Only 'manybody' is implemented.")
        ordered = [el for el, _ in sorted(type_table.items(), key=lambda kv: kv[1])]
        mapped = [el if el in self.elements else "NULL" for el in ordered]
        files = " ".join(_lmp_path(p) for p in self.param_files)
        return f"pair_coeff * * {files} {' '.join(mapped)}".rstrip()

    def missing_files(self) -> List[str]:
        return [p for p in self.param_files if not os.path.isfile(p)]

    def short_citation(self) -> str:
        return self.citation.split(".")[0] if self.citation else self.key


# ── known sidecar fields (anything else in the JSON is ignored gracefully) ──────
_FIELDS = {
    "key", "name", "elements", "pair_style", "pair_coeff_form", "required_units",
    "required_atom_style", "requires_qeq", "qeq_fix", "needs_packages",
    "recommended_timestep_fs", "default_temperature_K", "applicable_materials",
    "validity_notes", "confidence", "citation", "source_url", "source_note",
    "license",
}


def _from_sidecar(path: str) -> Potential:
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    d_dir = os.path.dirname(path)
    files = [os.path.join(d_dir, fn) for fn in d.get("param_files", [])]
    kwargs = {k: v for k, v in d.items() if k in _FIELDS}
    return Potential(param_files=files, sidecar_path=path, **kwargs)


def load_potentials(root: str = _HERE) -> dict:
    """Scan ``root`` recursively for ``*.json`` sidecars → {key: Potential}."""
    pots: dict = {}
    for js in sorted(glob.glob(os.path.join(root, "**", "*.json"), recursive=True)):
        try:
            pot = _from_sidecar(js)
        except (KeyError, json.JSONDecodeError, TypeError) as exc:
            # A malformed sidecar must not take down the whole library.
            import warnings
            warnings.warn(f"Skipping potential sidecar {js}: {exc}")
            continue
        if pot.key in pots:
            import warnings
            warnings.warn(f"Duplicate potential key '{pot.key}' "
                          f"({js} vs {pots[pot.key].sidecar_path}); keeping first.")
            continue
        pots[pot.key] = pot
    return pots


def potentials_for_material(material_key: str) -> List[Potential]:
    """All curated potentials whose sidecar lists this material as applicable."""
    return [p for p in POTENTIALS.values() if material_key in p.applicable_materials]


def potentials_covering(elements) -> List[Potential]:
    """All curated potentials whose element set covers every element in ``elements``."""
    els = list(elements)
    return [p for p in POTENTIALS.values() if p.covers_all(els)]


POTENTIALS = load_potentials()


if __name__ == "__main__":
    # quick introspection: `python -m potentials.registry` from SampleAssembler/
    if not POTENTIALS:
        print("No potentials curated yet.")
    for key, p in POTENTIALS.items():
        miss = p.missing_files()
        flag = "  !! MISSING FILES" if miss else ""
        print(f"{key:28s} [{p.pair_style:10s}] elements={p.elements} "
              f"units={p.required_units} conf={p.confidence}{flag}")
        if miss:
            for m in miss:
                print(f"     missing: {m}")
