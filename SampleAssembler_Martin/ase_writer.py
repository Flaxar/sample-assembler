"""Build atom coordinates with build123d-compatible geometry and write ASE files."""
from __future__ import annotations

import math
from itertools import product
from dataclasses import dataclass, field

import numpy as np
from ase import Atoms
from ase.io import write

import geometry as geo
from geometry import Box, Sphere, Hemisphere, Cylinder, Cone, ImportedRegion, Fill
from materials import MATERIALS, ELEMENT_MASS


@dataclass
class WriterOptions:
    """Options shared by the GUI and the ASE exporter."""

    margin: float = 0.0
    box_mode: str = "auto"
    box: tuple | None = None
    delete_overlap: bool = True
    overlap_cutoff: float = 1.2
    enforce_min_sep: bool = True
    random_maxtry: int = 50
    min_sep_atom_cap: int = 150000
    dump_file: str = "sample.lmp"
    sort_ids: bool = True
    seed_base: int = 12345
    boundary: str = "p p p"


@dataclass
class BuildResult:
    atoms: Atoms
    type_table: dict
    est_atoms: int
    per_fill_atoms: list
    warnings: list[str] = field(default_factory=list)

    @property
    def script(self):
        """Keep the existing preview pane useful without generating an input script."""
        return _describe_atoms(self)


def _region_volume(region):
    if isinstance(region, Box):
        return abs((region.xhi - region.xlo) * (region.yhi - region.ylo) *
                   (region.zhi - region.zlo))
    if isinstance(region, Sphere):
        return 4 / 3 * math.pi * region.r ** 3
    if isinstance(region, Hemisphere):
        return 2 / 3 * math.pi * region.r ** 3
    if isinstance(region, Cylinder):
        return math.pi * region.r ** 2 * abs(region.hi - region.lo)
    if isinstance(region, Cone):
        h = abs(region.hi - region.lo)
        return math.pi * h / 3 * (region.radlo ** 2 +
                                  region.radlo * region.radhi + region.radhi ** 2)
    if isinstance(region, ImportedRegion):
        return math.prod((region.aabb()[1] - region.aabb()[0],
                          region.aabb()[3] - region.aabb()[2],
                          region.aabb()[5] - region.aabb()[4]))
    return 0.0


def estimate_atoms(fill):
    mat = MATERIALS[fill.material_key]
    if mat.kind == "amorphous":
        return int(round(mat.atoms_per_a3() * _region_volume(fill.region)))
    a1, a2, a3 = mat.a1, mat.a2, mat.a3
    det = (a1[0] * (a2[1] * a3[2] - a2[2] * a3[1])
           - a1[1] * (a2[0] * a3[2] - a2[2] * a3[0])
           + a1[2] * (a2[0] * a3[1] - a2[1] * a3[0]))
    volume = abs(det) * mat.a ** 3
    return int(round(len(mat.basis) * _region_volume(fill.region) / volume))


def build_type_table(fills):
    elements = set()
    for fill in fills:
        if fill.enabled:
            elements.update(MATERIALS[fill.material_key].elements())
    ordered = sorted(elements, key=lambda e: (ELEMENT_MASS.get(e, 0.0), e))
    return {element: index + 1 for index, element in enumerate(ordered)}


def _rotation(fill, seed):
    return geo.resolve_rotation(fill.rotation, rng=np.random.default_rng(seed))


def _cell_vectors(mat, fill, scale):
    rotation = _rotation(fill, seed=fill.seed)
    vectors = np.asarray([mat.a1, mat.a2, mat.a3], dtype=float) * scale
    if rotation is not None:
        vectors = vectors @ np.asarray(rotation, dtype=float).T
    return vectors


def _inside(region, point, bounds=None):
    x, y, z = point
    if bounds is not None:
        if not (bounds[0] - 1e-9 <= x <= bounds[1] + 1e-9 and
                bounds[2] - 1e-9 <= y <= bounds[3] + 1e-9 and
                bounds[4] - 1e-9 <= z <= bounds[5] + 1e-9):
            return False
    if isinstance(region, ImportedRegion):
        from build123d import Vector
        return region.shape.is_inside(Vector(float(x), float(y), float(z)), tolerance=1e-6)
    if isinstance(region, Box):
        return region.xlo - 1e-9 <= x <= region.xhi + 1e-9 and \
            region.ylo - 1e-9 <= y <= region.yhi + 1e-9 and \
            region.zlo - 1e-9 <= z <= region.zhi + 1e-9
    if isinstance(region, Sphere):
        if not (region.clip_xlo is None or
                region.clip_xlo <= x <= region.clip_xhi and
                region.clip_ylo <= y <= region.clip_yhi):
            return False
        return (x - region.cx) ** 2 + (y - region.cy) ** 2 + \
            (z - region.cz) ** 2 <= region.r ** 2 + 1e-9
    if isinstance(region, Hemisphere):
        distance = (x - region.cx) ** 2 + (y - region.cy) ** 2 + \
            (z - region.cz) ** 2
        coordinate = {"x": x, "y": y, "z": z}[region.axis]
        center = {"x": region.cx, "y": region.cy, "z": region.cz}[region.axis]
        half = coordinate >= center if region.direction == "+" else coordinate <= center
        return distance <= region.r ** 2 + 1e-9 and half
    if isinstance(region, Cylinder):
        coordinate = {"x": x, "y": y, "z": z}[region.axis]
        radial = {"x": (y - region.c1) ** 2 + (z - region.c2) ** 2,
                  "y": (x - region.c1) ** 2 + (z - region.c2) ** 2,
                  "z": (x - region.c1) ** 2 + (y - region.c2) ** 2}[region.axis]
        return min(region.lo, region.hi) <= coordinate <= max(region.lo, region.hi) \
            and radial <= region.r ** 2 + 1e-9
    if isinstance(region, Cone):
        coordinate = {"x": x, "y": y, "z": z}[region.axis]
        radial = {"x": (y - region.c1) ** 2 + (z - region.c2) ** 2,
                  "y": (x - region.c1) ** 2 + (z - region.c2) ** 2,
                  "z": (x - region.c1) ** 2 + (y - region.c2) ** 2}[region.axis]
        lo, hi = sorted((region.lo, region.hi))
        if not lo - 1e-9 <= coordinate <= hi + 1e-9:
            return False
        t = (coordinate - region.lo) / (region.hi - region.lo)
        radius = region.radlo + (region.radhi - region.radlo) * t
        return radial <= radius ** 2 + 1e-9
    raise TypeError(f"Unsupported region type: {type(region).__name__}")


def _lattice_points(fill, box, rng):
    mat = MATERIALS[fill.material_key]
    scale = mat.a
    if mat.alloy is not None and fill.alloy_x is not None:
        scale = mat.alloy.a_of_x(fill.alloy_x, mat.a)
    vectors = _cell_vectors(mat, fill, scale)
    inverse = np.linalg.inv(vectors.T)
    # Scan only the fill's own AABB. Scanning the union box for every STEP
    # component makes a multi-component assembly needlessly expensive.
    region_box = fill.aabb()
    corners = np.array(list(product((region_box[0], region_box[1]),
                                    (region_box[2], region_box[3]),
                                    (region_box[4], region_box[5]))), dtype=float)
    # Use all box corners in fractional coordinates; this also covers rotated cells.
    fractional = corners @ inverse.T
    low = np.floor(fractional.min(axis=0)).astype(int) - 1
    high = np.ceil(fractional.max(axis=0)).astype(int) + 1
    origin = np.zeros(3)
    positions, symbols = [], []
    bounds = fill.aabb()
    for i in range(low[0], high[0] + 1):
        for j in range(low[1], high[1] + 1):
            for k in range(low[2], high[2] + 1):
                base = origin + i * vectors[0] + j * vectors[1] + k * vectors[2]
                for basis, element in zip(mat.basis, mat.basis_elements):
                    point = base + np.asarray(basis) @ vectors
                    if _inside(fill.region, point, bounds):
                        positions.append(point)
                        symbols.append(element)
    if mat.alloy is not None and fill.alloy_x:
        host, guest = mat.alloy.host, mat.alloy.guest
        for index, symbol in enumerate(symbols):
            if symbol == host and rng.random() < fill.alloy_x:
                symbols[index] = guest
    return positions, symbols


def _random_points(fill, box, rng, warnings, options):
    mat = MATERIALS[fill.material_key]
    count = estimate_atoms(fill)
    positions, symbols = [], []
    if count > options.min_sep_atom_cap:
        warnings.append(f"'{fill.name}': ~{count:,} amorphous atoms; min-separation skipped.")
        enforce = False
    else:
        enforce = options.enforce_min_sep
    minimum = mat.min_sep
    bounds = fill.aabb()
    placed_grid = {}
    for _ in range(count):
        for _attempt in range(options.random_maxtry if enforce else 1):
            point = rng.uniform([bounds[0], bounds[2], bounds[4]],
                                [bounds[1], bounds[3], bounds[5]])
            if not _inside(fill.region, point, bounds):
                continue
            if enforce:
                key = tuple(np.floor(point / minimum).astype(int))
                neighbors = (
                    placed_grid.get((key[0] + dx, key[1] + dy, key[2] + dz), [])
                    for dx in (-1, 0, 1)
                    for dy in (-1, 0, 1)
                    for dz in (-1, 0, 1)
                )
                if any(np.linalg.norm(point - old) < minimum
                       for bucket in neighbors for old in bucket):
                    continue
            positions.append(point)
            if enforce:
                placed_grid.setdefault(key, []).append(point)
            break
    formula = [element for element, amount in mat.stoich for _ in range(amount)]
    symbols = [formula[index % len(formula)] for index in range(len(positions))]
    rng.shuffle(symbols)
    return positions, symbols


def _remove_overlaps(positions, symbols, cutoff):
    kept_positions, kept_symbols = [], []
    grid = {}
    for point, symbol in zip(positions, symbols):
        key = tuple(np.floor(np.asarray(point) / cutoff).astype(int))
        nearby = (grid.get((key[0] + dx, key[1] + dy, key[2] + dz), [])
                  for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1))
        if any(np.linalg.norm(point - kept_positions[i]) < cutoff
               for bucket in nearby for i in bucket):
            continue
        index = len(kept_positions)
        kept_positions.append(point)
        kept_symbols.append(symbol)
        grid.setdefault(key, []).append(index)
    return kept_positions, kept_symbols


def build_atoms(fills, options, progress=None):
    progress = progress or (lambda _value, _message: None)
    enabled = [fill for fill in fills if fill.enabled]
    type_table = build_type_table(enabled)
    if options.box_mode == "explicit" and options.box:
        box = tuple(options.box)
    else:
        bounds = geo.union_aabb(enabled)
        box = tuple(value + (-options.margin if i % 2 == 0 else options.margin)
                    for i, value in enumerate(bounds))
    positions, symbols, warnings, per_fill = [], [], [], []
    for index, fill in enumerate(enabled):
        progress(5 + int(75 * index / max(len(enabled), 1)),
                 f"Generating atoms for {fill.name}...")
        mat = MATERIALS[fill.material_key]
        per_fill.append((fill.name, estimate_atoms(fill)))
        rng = np.random.default_rng(options.seed_base + index + 1)
        if mat.kind == "crystal":
            new_positions, new_symbols = _lattice_points(fill, box, rng)
        else:
            new_positions, new_symbols = _random_points(fill, box, rng, warnings, options)
        if fill.carve and positions:
            keep = [not _inside(fill.region, point, fill.aabb()) for point in positions]
            positions = [point for point, keep_point in zip(positions, keep) if keep_point]
            symbols = [symbol for symbol, keep_point in zip(symbols, keep) if keep_point]
        positions.extend(new_positions)
        symbols.extend(new_symbols)
        progress(5 + int(75 * (index + 1) / max(len(enabled), 1)),
                 f"Finished {fill.name}.")
    if options.delete_overlap and positions:
        progress(85, "Removing overlapping atoms...")
        positions, symbols = _remove_overlaps(positions, symbols, options.overlap_cutoff)
    shift = np.array([box[0], box[2], box[4]], dtype=float)
    atoms = Atoms(symbols=symbols, positions=np.asarray(positions, dtype=float) - shift,
                  cell=np.diag([box[1] - box[0], box[3] - box[2], box[5] - box[4]]),
                  pbc=False)
    return BuildResult(atoms, type_table, len(atoms), per_fill, warnings)


def export_lammps(path, fills, options, progress=None):
    result = build_atoms(fills, options, progress)
    progress = progress or (lambda _value, _message: None)
    progress(95, "Writing LAMMPS data file...")
    specorder = [element for element, _ in
                 sorted(result.type_table.items(), key=lambda item: item[1])]
    write(path, result.atoms, format="lammps-data", atom_style="atomic",
          units="metal", specorder=specorder, masses=True)
    progress(100, "LAMMPS data export complete.")
    return result


def export_xyz(path, fills, options, progress=None):
    """Generate the sample and write it as a standard ASE XYZ file."""
    result = build_atoms(fills, options, progress)
    progress = progress or (lambda _value, _message: None)
    progress(95, "Writing XYZ file...")
    write(path, result.atoms, format="xyz")
    progress(100, "XYZ export complete.")
    return result


def write_minimization_input(path, data_path, dump_path, potential, type_table,
                             boundary="p p p", overlap_cutoff=1.2,
                             workflow="relax", temperature=800.0,
                             equilibration_steps=5000, timestep_fs=1.0,
                             seed=87287):
    """Write a LAMMPS input which cleans overlaps and minimizes the exported data.

    The data file and dump are resolved relative to the input file so the
    generated workflow can be moved as one self-contained set.
    """
    import os

    input_dir = os.path.dirname(os.path.abspath(path))
    data_name = os.path.basename(data_path)
    dump_name = os.path.basename(dump_path)
    potential_file = potential.param_files[0]
    if not os.path.isfile(potential_file):
        raise FileNotFoundError(f"Potential file not found: {potential_file}")
    potential_name = os.path.basename(potential_file)
    if os.path.dirname(os.path.abspath(potential_file)) != input_dir:
        import shutil
        copied = os.path.join(input_dir, potential_name)
        if os.path.abspath(copied) != os.path.abspath(potential_file):
            shutil.copy2(potential_file, copied)

    ordered_elements = [
        element for element, _ in sorted(type_table.items(), key=lambda item: item[1])
    ]
    mapped = [element if potential.covers(element) else "NULL"
              for element in ordered_elements]
    if any(element == "NULL" for element in mapped):
        missing = [element for element in ordered_elements if not potential.covers(element)]
        raise ValueError(f"Selected potential does not cover: {', '.join(missing)}")

    lines = [
        "# Generated by Sample Assembler",
        "units metal",
        "atom_style atomic",
        f"boundary {boundary}",
        f"read_data {data_name}",
        "",
        f"pair_style {potential.pair_style}",
        f"pair_coeff * * {potential_name} {' '.join(mapped)}",
        "",
        f"delete_atoms overlap {overlap_cutoff:.6g} all all",
        "thermo 500",
    ]
    if workflow in ("thermalize", "anneal"):
        timestep_ps = timestep_fs / 1000.0
        lines.extend([
            f"timestep {timestep_ps:.8g}",
            f"velocity all create {temperature:.8g} {seed:d} loop geom",
            f"fix anneal_heat all nvt temp {temperature:.8g} {temperature:.8g} 0.1",
            f"run {equilibration_steps:d}",
            "unfix anneal_heat",
            f"fix anneal_cool all nvt temp {temperature:.8g} 1.0 0.1",
            f"run {equilibration_steps:d}",
            "unfix anneal_cool",
        ])
    lines.extend([
        "",
        "min_style cg",
        "minimize 1.0e-6 1.0e-8 1000 10000",
        "",
        f"dump 1 all xyz 1 {dump_name}",
        f"dump_modify 1 element {' '.join(mapped)} sort id",
        "run 0",
        "undump 1",
    ])
    with open(path, "w", encoding="utf-8", newline="\n") as stream:
        stream.write("\n".join(lines) + "\n")


def _describe_atoms(result):
    lines = ["# ASE-generated LAMMPS data preview",
             f"# atoms: {len(result.atoms):,}",
             "# atom types: " + ", ".join(
                 f"{element}={number}" for element, number in result.type_table.items())]
    for name, count in result.per_fill_atoms:
        lines.append(f"# {name}: ~{count:,}")
    lines.extend(f"# warning: {warning}" for warning in result.warnings)
    return "\n".join(lines)
