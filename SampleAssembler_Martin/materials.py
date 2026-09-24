#!/usr/bin/env python3
"""
Curated crystallography library for the Sample Assembler.

Each Material provides the cell data used by the ASE writer. Crystalline
materials are described by a unit cell (a1,a2,a3 in units of the scale `a`,
with a3_z = c/a for hexagonal/tetragonal cells) plus a fractional-coordinate
basis tagged with the element at each site. Amorphous materials carry a mass
density which the writer converts to a target atom count.

EVERYTHING is keyed by element symbol. LAMMPS atom *types* are global to a
simulation, so `lammps_writer` builds one global element->type table across all
the materials actually used and resolves each basis site to that global id.

Lattice rotation: a crystal is rotated by rotating its a1/a2/a3 vectors with a
3x3 matrix R (see geometry.rotation_matrix); the *fractional* basis is left
unchanged, which is the correct way to reorient a lattice. Arbitrary SO(3)
rotations (e.g. randomly oriented nanoparticles) are therefore supported
directly, in addition to LAMMPS `orient` integer-Miller reorientation.

Confidence tags ('high'/'med'/'low') and notes are attached for the fact-check
pass. Values are conventional-cell, ~300 K.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Optional

N_AVOGADRO = 6.02214076e23
# g/cm^3 -> atoms/Å^3 :  n = rho * N_A / M  [1/cm^3] ; 1 cm^3 = 1e24 Å^3
_CM3_PER_A3 = 1.0e-24

# ── Standard atomic masses (amu) ───────────────────────────────────────────────
ELEMENT_MASS = {
    "H": 1.008,  "B": 10.811, "C": 12.011, "N": 14.007, "O": 15.999,
    "Al": 26.982, "Si": 28.0855, "P": 30.974, "S": 32.06, "Ti": 47.867,
    "Cu": 63.546, "Ga": 69.723, "Ge": 72.630, "As": 74.9216, "Sr": 87.62,
    "Mo": 95.95, "In": 114.818, "Hf": 178.49, "W": 183.84, "Au": 196.967,
}

# ── Reusable fractional bases ──────────────────────────────────────────────────
_FCC4 = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
_BCC2 = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]
_FCC_OFFSET_QUARTER = [(0.25, 0.25, 0.25), (0.75, 0.75, 0.25),
                       (0.75, 0.25, 0.75), (0.25, 0.75, 0.75)]
_FCC_OFFSET_HALF_X = [(0.5, 0.0, 0.0), (0.0, 0.5, 0.0),
                      (0.0, 0.0, 0.5), (0.5, 0.5, 0.5)]
_DIAMOND8 = _FCC4 + _FCC_OFFSET_QUARTER
_FLUORITE_ANION = _FCC_OFFSET_QUARTER + [
    (0.75, 0.75, 0.75), (0.25, 0.25, 0.75), (0.25, 0.75, 0.25), (0.75, 0.25, 0.25)]

_S3_2 = math.sqrt(3.0) / 2.0          # 0.8660254…  (120° hexagonal)
_HEX_A1 = (1.0, 0.0, 0.0)
_HEX_A2 = (-0.5, _S3_2, 0.0)
_CUBE_A1, _CUBE_A2, _CUBE_A3 = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)


@dataclass
class AlloyInfo:
    """A substitutional solid solution: a fraction `x` of `host` sites become `guest`."""
    host: str
    guest: str
    default_x: float = 0.5
    # Vegard-style a(x) in Å (coeffs for a0 + a1*x + a2*x^2); None -> use base a.
    a_coeffs: Optional[tuple] = None

    def a_of_x(self, x: float, base_a: float) -> float:
        if self.a_coeffs is None:
            return base_a
        c0, c1, c2 = self.a_coeffs
        return c0 + c1 * x + c2 * x * x


@dataclass
class Material:
    key: str
    name: str
    formula: str
    category: str            # 'semiconductor' | 'dielectric' | 'metal' | '2D' | 'amorphous'
    kind: str                # 'crystal' | 'amorphous'
    confidence: str = "high"
    notes: str = ""

    # crystal fields
    a: float = 0.0
    c: Optional[float] = None
    a1: tuple = _CUBE_A1
    a2: tuple = _CUBE_A2
    a3: tuple = _CUBE_A3
    basis: list = field(default_factory=list)          # [(x,y,z), …] fractional
    basis_elements: list = field(default_factory=list)  # parallel element symbols
    alloy: Optional[AlloyInfo] = None
    layered: bool = False     # 2D material — keep region ~1 cell thick in z

    # amorphous fields
    density_g_cm3: float = 0.0
    stoich: list = field(default_factory=list)          # [(element, count_per_formula), …]
    min_sep: float = 1.4     # Å: min separation enforced during random placement

    # ── element / mass introspection ──────────────────────────────────────────
    def elements(self) -> set:
        els = set()
        if self.kind == "crystal":
            els.update(self.basis_elements)
            if self.alloy:
                els.add(self.alloy.host)
                els.add(self.alloy.guest)
        else:
            els.update(e for e, _ in self.stoich)
        return els

    # ── crystal: build the LAMMPS lattice line ─────────────────────────────────
    def lattice_line(self, scale: float, R=None, origin=None) -> str:
        """`lattice custom <scale> a1 .. a2 .. a3 .. [origin ..] basis ..`.

        ``R`` is an optional 3x3 rotation applied to a1/a2/a3. ``origin`` is an
        optional (ox,oy,oz) triple of *lattice-spacing fractions* (0<=o<1) that
        rigidly shifts the whole lattice — see ``WriterOptions.lattice_origin``
        for why a rotated lattice needs one.

        The a-vector components are emitted at full double precision: a rotated
        lattice has irrational direction cosines, and truncating them to 6 dp
        accumulates ~1e-5 A of drift across a 50 A box, which is enough to push
        a boundary atom plane outside the box.
        """
        a1, a2, a3 = self.a1, self.a2, self.a3
        if R is not None:
            a1, a2, a3 = _apply(R, a1), _apply(R, a2), _apply(R, a3)
        parts = [f"lattice custom {scale:.10g}",
                 _vec("a1", a1), _vec("a2", a2), _vec("a3", a3)]
        if origin is not None:
            ox, oy, oz = origin
            parts.append(f"origin {ox:.6f} {oy:.6f} {oz:.6f}")
        for (x, y, z) in self.basis:
            parts.append(f"basis {x:.6f} {y:.6f} {z:.6f}")
        return " ".join(parts)

    def c_over_a(self) -> Optional[float]:
        return None if self.c is None else self.c / self.a

    # ── amorphous: number density ───────────────────────────────────────────────
    def formula_mass(self) -> float:
        return sum(ELEMENT_MASS[e] * n for e, n in self.stoich)

    def number_density_formula(self) -> float:
        """Formula units per Å^3 from mass density."""
        if self.kind != "amorphous":
            return 0.0
        n_per_cm3 = self.density_g_cm3 * N_AVOGADRO / self.formula_mass()
        return n_per_cm3 * _CM3_PER_A3

    def atoms_per_a3(self) -> float:
        if self.kind != "amorphous":
            return 0.0
        per_formula = self.number_density_formula()
        return per_formula * sum(n for _, n in self.stoich)


def _vec(name, v):
    return f"{name} {v[0]:.12f} {v[1]:.12f} {v[2]:.12f}"


def _apply(R, v):
    return (R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
            R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
            R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2])


# ── Builders for the common crystal families ───────────────────────────────────
def _diamond(el):
    return list(_DIAMOND8), [el] * 8


def _fcc(el):
    return list(_FCC4), [el] * 4


def _bcc(el):
    return list(_BCC2), [el] * 2


def _zincblende(cation, anion):
    return list(_DIAMOND8), [cation] * 4 + [anion] * 4


def _rocksalt(cation, anion):
    return _FCC4 + _FCC_OFFSET_HALF_X, [cation] * 4 + [anion] * 4


def _fluorite(cation, anion):
    return _FCC4 + _FLUORITE_ANION, [cation] * 4 + [anion] * 8


# ═══════════════════════════════════════════════════════════════════════════════
#  MATERIAL LIBRARY
# ═══════════════════════════════════════════════════════════════════════════════
#
#  ▸▸ THIS LIBRARY IS DESIGNED TO GROW — MORE MATERIALS WILL BE ADDED. ◂◂
#
#  To add a new material, append one `add(Material(...))` below and it appears in
#  the GUI material picker automatically (grouped by `category`). Recipes:
#    • cubic single element : basis,elements = _diamond/_fcc/_bcc("El")
#    • binary compound      : _zincblende / _rocksalt / _fluorite(cation, anion)
#    • hexagonal/tetragonal : a1=_HEX_A1, a2=_HEX_A2, a3=(0,0,c/a) + explicit basis
#    • solid solution       : alloy=AlloyInfo(host, guest, x, Vegard coeffs)
#    • amorphous            : kind="amorphous", density_g_cm3 + stoich=[(El,count)]
#  Add any new element's mass to ELEMENT_MASS. Set `confidence`/`notes` honestly.
#
#  PLANNED (not yet implemented): more III-V/II-VI compounds, Mo/Co/Ru interconnect
#  metals, MoS2/WS2 2D channels, monoclinic HfO2, Si3N4, and additional sample
#  presets beyond the crossgrating. Keep this structure open for them.
# ═══════════════════════════════════════════════════════════════════════════════
def _build_library():
    M = {}

    def add(mat):
        M[mat.key] = mat

    # ---- Group IV semiconductors ----
    b, e = _diamond("Si")
    add(Material("Si", "Silicon", "Si", "semiconductor", "crystal",
                 a=5.4309, basis=b, basis_elements=e,
                 notes="Diamond cubic Fd-3m, a=5.4309 Å @300K."))

    b, e = _diamond("Ge")
    add(Material("Ge", "Germanium", "Ge", "semiconductor", "crystal",
                 a=5.6575, basis=b, basis_elements=e,
                 notes="Diamond cubic, a=5.6575 Å. Higher-mobility channel."))

    b, e = _diamond("Si")
    add(Material("SiGe", "Silicon-Germanium (alloy)", "Si(1-x)Ge(x)",
                 "semiconductor", "crystal",
                 a=5.4309, basis=b, basis_elements=e,
                 alloy=AlloyInfo("Si", "Ge", 0.30, (5.4309, 0.1992, 0.0272)),
                 notes="Diamond cubic Si with a random Ge fraction x (set type/fraction). "
                       "Vegard a(x)=5.4309+0.1992x+0.0272x^2 (Dismukes 1964)."))

    b, e = _zincblende("Si", "C")
    add(Material("SiC", "Silicon Carbide (3C / β-SiC)", "SiC",
                 "semiconductor", "crystal",
                 a=4.3596, basis=b, basis_elements=e,
                 notes="Cubic zinc-blende 3C-SiC (β-SiC), F-43m, a=4.3596 Å @300K. "
                       "IV-IV compound; Si and C on the two interpenetrating FCC sublattices."))

    # ---- Carbon allotropes ----
    b, e = _diamond("C")
    add(Material("Cdiamond", "Diamond (C)", "C", "semiconductor", "crystal",
                 a=3.5668, basis=b, basis_elements=e,
                 notes="Diamond cubic carbon, a=3.5668 Å."))

    add(Material("Graphite", "Graphite (C)", "C", "2D", "crystal",
                 a=2.461, c=6.708, a1=_HEX_A1, a2=_HEX_A2, a3=(0.0, 0.0, 6.708 / 2.461),
                 basis=[(0.0, 0.0, 0.0), (1/3, 2/3, 0.0),
                        (0.0, 0.0, 0.5), (2/3, 1/3, 0.5)],
                 basis_elements=["C", "C", "C", "C"], confidence="med",
                 notes="Hexagonal AB (Bernal) stacking, a=2.461, c=6.708 Å. "
                       "Interlayer 3.354 Å; stacking detail approximate."))

    add(Material("Graphene", "Graphene (single layer)", "C", "2D", "crystal",
                 a=2.461, c=3.354, a1=_HEX_A1, a2=_HEX_A2, a3=(0.0, 0.0, 3.354 / 2.461),
                 basis=[(0.0, 0.0, 0.0), (1/3, 2/3, 0.0)],
                 basis_elements=["C", "C"], layered=True,
                 notes="Honeycomb, a=2.461 Å, 2-atom basis. Keep region ~1 cell thick in z."))

    # ---- Oxide / perovskite ----
    add(Material("SrTiO3", "Strontium Titanate", "SrTiO3", "dielectric", "crystal",
                 a=3.905, basis=[(0.0, 0.0, 0.0), (0.5, 0.5, 0.5),
                                 (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)],
                 basis_elements=["Sr", "Ti", "O", "O", "O"],
                 notes="Cubic perovskite Pm-3m, a=3.905 Å. Sr corner, Ti body-center, O face-centers."))

    # ---- Nitrides ----
    add(Material("hBN", "Hexagonal Boron Nitride", "BN", "dielectric", "crystal",
                 a=2.504, c=6.661, a1=_HEX_A1, a2=_HEX_A2, a3=(0.0, 0.0, 6.661 / 2.504),
                 basis=[(0.0, 0.0, 0.0), (1/3, 2/3, 0.0),
                        (0.0, 0.0, 0.5), (1/3, 2/3, 0.5)],
                 basis_elements=["B", "N", "N", "B"], confidence="med",
                 notes="P6_3/mmc, a=2.504, c=6.661 Å, AA' stacking (B over N). "
                       "NOTE: uses c/a in a3 (corrects c-spacing vs the old TEM.lmp)."))

    # ---- III-V ----
    b, e = _zincblende("Ga", "As")
    add(Material("GaAs", "Gallium Arsenide", "GaAs", "semiconductor", "crystal",
                 a=5.653, basis=b, basis_elements=e,
                 notes="Zinc-blende F-43m, a=5.653 Å. Standard III-V TEM specimen."))

    add(Material("GaN_wz", "Gallium Nitride (wurtzite)", "GaN", "semiconductor", "crystal",
                 a=3.189, c=5.185, a1=_HEX_A1, a2=_HEX_A2, a3=(0.0, 0.0, 5.185 / 3.189),
                 basis=[(1/3, 2/3, 0.0), (2/3, 1/3, 0.5),
                        (1/3, 2/3, 0.377), (2/3, 1/3, 0.877)],
                 basis_elements=["Ga", "Ga", "N", "N"],
                 notes="Wurtzite P6_3mc, a=3.189, c=5.185 Å, internal u=0.377."))

    b, e = _zincblende("Ga", "N")
    add(Material("GaN_zb", "Gallium Nitride (zinc-blende)", "GaN", "semiconductor", "crystal",
                 a=4.50, basis=b, basis_elements=e, confidence="med",
                 notes="Metastable zinc-blende GaN, a≈4.50 Å."))

    b, e = _zincblende("In", "P")
    add(Material("InP", "Indium Phosphide", "InP", "semiconductor", "crystal",
                 a=5.8687, basis=b, basis_elements=e,
                 notes="Zinc-blende, a=5.8687 Å."))

    b, e = _zincblende("In", "As")
    add(Material("InAs", "Indium Arsenide", "InAs", "semiconductor", "crystal",
                 a=6.0583, basis=b, basis_elements=e,
                 notes="Zinc-blende, a=6.0583 Å."))

    b, e = _zincblende("Ga", "As")  # cation = Ga, substitute In
    add(Material("InGaAs", "Indium Gallium Arsenide (alloy)", "In(x)Ga(1-x)As",
                 "semiconductor", "crystal",
                 a=5.653, basis=b, basis_elements=e,
                 alloy=AlloyInfo("Ga", "In", 0.53, (5.6533, 0.405, 0.0)),
                 notes="GaAs with In substituting Ga on the cation sublattice. "
                       "x=0.53 -> a≈5.868 Å (lattice-matched to InP). a(x)=5.6533+0.405x."))

    # ---- α-Ga (elemental metal) ----
    add(Material("Ga_alpha", "Gallium (α-Ga metal)", "Ga", "metal", "crystal",
                 a=4.5258, c=4.5198, a1=_CUBE_A1, a2=(0.0, 7.6602 / 4.5258, 0.0),
                 a3=(0.0, 0.0, 4.5198 / 4.5258),
                 basis=[(0.0, 0.1549, 0.0810), (0.0, 0.8451, 0.9190),
                        (0.0, 0.1549, 0.4190), (0.0, 0.8451, 0.5810),
                        (0.5, 0.6549, 0.0810), (0.5, 0.3451, 0.9190),
                        (0.5, 0.6549, 0.4190), (0.5, 0.3451, 0.5810)],
                 basis_elements=["Ga"] * 8, confidence="low",
                 notes="Orthorhombic Cmce α-Ga, a=4.5258 b=7.6602 c=4.5198 Å, Ga at "
                       "8f(0,0.1549,0.0810): orbit (0,y,z),(0,-y,-z),(0,y,1/2-z),"
                       "(0,-y,1/2+z) + C-centering (1/2,1/2,0). Low confidence — verify."))

    # ---- Metals (gate / interconnect / test particle) ----
    b, e = _fcc("Au")
    add(Material("Au", "Gold", "Au", "metal", "crystal",
                 a=4.0782, basis=b, basis_elements=e,
                 notes="FCC, a=4.0782 Å. (111) spacing 0.235 nm — the TEM resolution-test fringe."))

    b, e = _fcc("Cu")
    add(Material("Cu", "Copper", "Cu", "metal", "crystal",
                 a=3.6149, basis=b, basis_elements=e,
                 notes="FCC, a=3.6149 Å. Interconnect."))

    b, e = _bcc("W")
    add(Material("W", "Tungsten", "W", "metal", "crystal",
                 a=3.1652, basis=b, basis_elements=e,
                 notes="BCC, a=3.1652 Å. Gate/contact metal."))

    b, e = _rocksalt("Ti", "N")
    add(Material("TiN", "Titanium Nitride", "TiN", "metal", "crystal",
                 a=4.2417, basis=b, basis_elements=e,
                 notes="Rocksalt (NaCl) Fm-3m, a=4.2417 Å. Metal-gate liner."))

    # ---- High-k dielectric ----
    b, e = _fluorite("Hf", "O")
    add(Material("HfO2", "Hafnia (cubic fluorite approx.)", "HfO2", "dielectric", "crystal",
                 a=5.08, basis=b, basis_elements=e, confidence="low",
                 notes="Fluorite (CaF2) approximation, a≈5.08 Å. The stable phase is "
                       "monoclinic baddeleyite; fluorite is a clean lattice for modelling."))

    # ---- Amorphous ----
    add(Material("aC", "Amorphous Carbon", "C", "amorphous", "amorphous",
                 density_g_cm3=2.0, stoich=[("C", 1)],
                 notes="Evaporated a-C support film, ρ≈2.0 g/cm³ (NIST) -> 0.1003 atoms/Å³."))

    add(Material("aSiO2", "Amorphous SiO2 (gate oxide)", "SiO2", "amorphous", "amorphous",
                 density_g_cm3=2.20, stoich=[("Si", 1), ("O", 2)], confidence="med",
                 notes="Fused-silica-like, ρ≈2.20 g/cm³."))

    add(Material("aAl2O3", "Amorphous Al2O3 (ALD)", "Al2O3", "amorphous", "amorphous",
                 density_g_cm3=3.00, stoich=[("Al", 2), ("O", 3)], confidence="med",
                 notes="ALD/amorphous alumina spacer, ρ≈3.0 g/cm³ (crystalline corundum is 3.95)."))

    return M


MATERIALS = _build_library()


def get_material(key: str) -> Material:
    return MATERIALS[key]


def materials_by_category():
    cats = {}
    for m in MATERIALS.values():
        cats.setdefault(m.category, []).append(m)
    for v in cats.values():
        v.sort(key=lambda m: m.name)
    return cats
