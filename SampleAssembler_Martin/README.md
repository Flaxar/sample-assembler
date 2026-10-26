# Sample Assembler

A compact, CAD-style desktop app for designing TEM specimens and exporting them
as **LAMMPS data files**. Geometry and atom coordinates are generated directly
with build123d-compatible containment checks and ASE; LAMMPS is not required
for export. Define arbitrary geometry, fill each region with a material from a
crystallography library (with lattice rotation and alloy fractions), and delete
overlaps before writing the data file.

It is a **standalone** app (the multislice GUI's modules are *not* a
dependency); it only reuses the same Obsidian dark theme for visual continuity.

## Run

```bash
cd SampleAssembler
python sample_assembler.py
```

Dependencies: `numpy`, `matplotlib`, `PyQt5`, `build123d`, and `ase`. STEP files
can be imported from panel 1 with **Import .step**; imported component names are
then available in the component selector below the button.

## Workflow

1. **Add region + material** (panel 1): pick a material, choose a shape
   (box / sphere / hemisphere / cylinder / cone) and its parameters in Å, set an
   optional lattice rotation and (for alloys) the guest fraction, optionally tick
   **Carve** so the region deletes whatever atoms already occupy its footprint,
   then *Add region to sample*.
2. **Sample regions** (panel 2): the stacking list, in build order. Toggle the
   checkbox to enable/disable, reorder with ↑/↓, remove, or clear. Order matters:
   later fills carve / overlap earlier ones.
3. **Build options** (panel 3): auto/explicit box, margin, periodic/shrink
   boundary, overlap deletion + cutoff, amorphous min-separation, dump filename.
4. **MD simulation** (panel 4): choose either the LAMMPS backend with a
   compatible curated potential, or the optional **MACE MLIP (ASE)** backend.
   The thermalize/anneal workflows
   create low-temperature velocities, gradually ramp to the selected
   temperature, cool to 1 K, minimize with conjugate gradients, and write a
   relaxed XYZ snapshot. The `relax`
   workflow skips the thermal stage and performs only the minimization.
5. **Export** (panel 5): write the ASE LAMMPS data file and save/load the design
   as a `.saproj` JSON project.

The right pane shows a live 3D preview, the generated `.lmp`, and an **Info & Log**
tab (global atom-type table, per-region atom-count estimate, warnings, run log).

When multiple compatible crystalline regions use similar lattice parameters, the
builder uses the smallest parameter for all of them. This is useful for strained
systems such as Si and SiGe: thermalization can then relax the shared lattice
instead of starting from an expanded alloy cell. Unrelated crystal structures
are not merged, and amorphous materials are unaffected.

## Materials library (`materials.py`)

> **This library is designed to grow** — more materials, elements, and sample
> presets are planned (see the `PLANNED` note in `materials.py`). Adding a
> material is a single `add(Material(...))` entry and it appears in the GUI
> automatically.

Curated, `~300 K` conventional-cell data, keyed by element so mixed-material
samples share one global LAMMPS atom-type table:

| Category | Materials |
|---|---|
| Semiconductor | Si, Ge, SiGe (alloy), GaAs, GaN (wurtzite & zinc-blende), InP, InAs, InGaAs (alloy), diamond C |
| Dielectric | SrTiO₃, hBN, HfO₂ (fluorite approx.) |
| Metal | Au, Cu, W, TiN, α-Ga |
| 2D | Graphene, Graphite |
| Amorphous | a-C, a-SiO₂, a-Al₂O₃ (density → atom count) |

Lattice rotation: **Euler**, **axis-angle**, **uniform random SO(3)**, or cubic
**Miller** orientation (`orient x/y/z` directions). A crystal is reoriented by
rotating its `a1/a2/a3` cell vectors; the fractional basis is unchanged.

## "Crossgrating" preset

One click builds the classic **gold-on-amorphous-carbon** TEM resolution/
astigmatism test specimen: an a-C support film decorated with randomly sized,
randomly oriented, optionally clustered Au nanoparticles (Au (111) = 0.235 nm
fringes). This is a faithful port of the original `lammps_TEM.py` placement
algorithm. Defaults are kept tractable (400×400×80 Å, ~40 NPs ≈ 1 M atoms);
a full-density a-C film over a large field of view is many millions of atoms
(the Info tab warns when a region is huge and falls back to fast placement).

## Generated LAMMPS data file

ASE writes a standard `lammps-data` file containing the global atom types,
masses, simulation cell, and all generated atomic coordinates. Crystalline
coordinates come from the material basis and rotated cell vectors; amorphous
coordinates are density-based random placements. Carving and overlap removal
are performed before ASE writes the file.

## Conventions

- All coordinates are **Å** (`units box`), matching the multislice convention.
- **Z axis:** the beam enters at `z_min` and exits at `z_max`. The crossgrating
  places the carbon film at `z ≈ 0` and nanoparticles above it.
- The exported coordinates are a clean build: equilibrium lattice positions for
  crystals and density-based random positions for amorphous materials. Optional
  LAMMPS minimization is a separate post-export step; it does not change the
  original data file.

## Interatomic potentials & MD (`potentials/`, `MD_PLAYBOOK.md`)

The `potentials/` folder is a curated library — a LAMMPS parameter file plus a JSON
sidecar (citation, elements, units, atom_style, confidence) per potential,
auto-loaded by `potentials/registry.py`. The library includes a Si-Ge Tersoff
potential for structural relaxation and the GaN potentials. **`MD_PLAYBOOK.md`**
is the reference for which potential / thermostat / ensemble / timestep to use
per material. Adding a potential =
drop the file + a sidecar JSON into `potentials/<family>/`.

MACE is optional because it requires PyTorch and a large model download. Install
`mace-torch` in the active Python environment to enable the MACE backend:

```bash
python -m pip install mace-torch
```

MACE uses ASE's LBFGS optimizer and writes `<sample>.traj` plus the relaxed XYZ
file. The LAMMPS backend remains available for curated non-ML potentials.

## Tools

- **Tools ▸ Split multi-timestep `.xyz`** — splits a multi-snapshot LAMMPS dump
  into per-timestep files for frozen-phonon multislice (ported from
  `xyz_parser.py`).

## Files

| File | Role |
|---|---|
| `sample_assembler.py` | PyQt5 GUI + ASE export (entry point) |
| `materials.py` | crystallography library and fractional bases |
| `geometry.py` | region primitives + rotation math + `Fill` model |
| `ase_writer.py` | build123d-aware atom generation and ASE LAMMPS export |
| `lammps_writer.py` | compatibility facade for `ase_writer.py` |
| `potentials/` | curated interatomic-potential library (param files + JSON sidecars) + `registry.py` |
| `MD_PLAYBOOK.md` | MD reference: per-material potential / thermostat / timestep / workflow guidance |
| `crossgrating.py` | gold-on-carbon preset (port of `lammps_TEM.py`) |
| `preview3d.py` | matplotlib 3D region preview |
| `theme.py` | self-contained Obsidian theme |
