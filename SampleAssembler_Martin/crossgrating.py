#!/usr/bin/env python3
"""
"Crossgrating" gold-on-carbon preset — the classic TEM resolution / astigmatism
test specimen: evaporated amorphous-carbon support film decorated with randomly
sized, randomly oriented, optionally clustered gold nanoparticles. Au (111)
fringes (0.235 nm) and the wide orientation spread make it the standard probe of
instrument resolution and alignment.

This is a faithful port of the original lammps_TEM.py placement algorithm
(clustered positions with a soft-overlap constraint, uniform-random SO(3)
orientation per particle via quaternions), repackaged as a list of Fill objects
so it flows through lammps_writer like any other assembled sample.

Defaults reflect real gold-on-carbon test specimens: a-C film a few–tens of nm
thick, Au islands ~2–10 nm (here 3–10 nm diameter), partially embedded.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from geometry import Box, Sphere, Fill, random_rotation_matrix


@dataclass
class CrossgratingParams:
    box_x: float = 400.0           # Å  (full-density a-C scales as V; large FOV -> millions of atoms)
    box_y: float = 400.0
    carbon_thickness: float = 80.0
    n_particles: int = 40
    r_min: float = 15.0            # Å (3 nm diameter)
    r_max: float = 40.0            # Å (8 nm diameter)
    cluster_strength: float = 0.6  # 0 = uniform, 1 = highly clustered
    overlap_factor: float = 0.4    # <1 allows soft overlap (0.7–0.8 ~ realistic)
    protrusion: float = 0.5        # fraction of each NP sitting above the C surface
    seed: int = 42
    max_tries: int = 5000
    carve_carbon: bool = True      # gold removes the carbon beneath it


def _place_positions(p: CrossgratingParams, radii, rng):
    """Clustered NP placement with a soft min-separation constraint (ported)."""
    positions = []
    n_clusters = max(1, int(p.n_particles * p.cluster_strength))
    cluster_spread = p.box_x * 0.15
    center = np.array([p.box_x / 2.0, p.box_y / 2.0])
    centers = center + rng.normal(scale=cluster_spread, size=(n_clusters, 2))

    n_overlapped = 0
    for i in range(p.n_particles):
        r_i = radii[i]
        placed = False
        for _ in range(p.max_tries):
            if rng.random() < p.cluster_strength:
                c = centers[rng.integers(len(centers))]
                pos = c + rng.normal(scale=p.box_x * 0.05, size=2)
            else:
                pos = rng.random(2) * np.array([p.box_x, p.box_y])
            pos[0] %= p.box_x
            pos[1] %= p.box_y
            ok = True
            for j, p_old in enumerate(positions):
                if np.linalg.norm(pos - p_old) < p.overlap_factor * (r_i + radii[j]):
                    ok = False
                    break
            if ok:
                positions.append(pos)
                placed = True
                break
        if not placed:
            positions.append(pos)
            n_overlapped += 1
    return np.array(positions), n_overlapped


def build_crossgrating(p: CrossgratingParams):
    """Return (fills, stats). fills[0] is the a-C film; the rest are Au NPs."""
    rng = np.random.default_rng(p.seed)
    radii = rng.uniform(p.r_min, p.r_max, p.n_particles)
    positions, n_overlapped = _place_positions(p, radii, rng)

    fills = [Fill(name="a-C support film",
                  region=Box(0.0, p.box_x, 0.0, p.box_y, 0.0, p.carbon_thickness),
                  material_key="aC", seed=p.seed)]

    for i in range(p.n_particles):
        x, y = positions[i]
        r = radii[i]
        # center z: protrusion=0 -> buried (top flush), 0.5 -> half above, 1 -> fully above
        cz = p.carbon_thickness - r + 2.0 * r * p.protrusion
        R = random_rotation_matrix(rng)
        # Clip each NP to the carbon footprint so particles straddling an edge are
        # cut flush (z left free → protrusion preserved); the sample's XY extent
        # then always equals the a-C substrate instead of being enlarged by an
        # overhanging sphere.
        fills.append(Fill(
            name=f"Au NP {i + 1}",
            region=Sphere(float(x), float(y), float(cz), float(r),
                          clip_xlo=0.0, clip_xhi=p.box_x,
                          clip_ylo=0.0, clip_yhi=p.box_y),
            material_key="Au",
            rotation={"mode": "matrix", "R": R.tolist()},
            seed=p.seed + i + 1,
            carve=p.carve_carbon,
        ))

    # stats (nearest-neighbour, mean overlap) — mirrors the original script
    stats = {"n_particles": p.n_particles, "n_overlapped": n_overlapped}
    if len(positions) > 1:
        nn, ov = [], []
        for i in range(len(positions)):
            d = np.linalg.norm(positions - positions[i], axis=1)
            d = d[d > 0]
            j = int(np.argmin(d))
            nn.append(d[j])
            # j here indexes the filtered array; recompute true neighbour radius safely
        stats["mean_nn_distance"] = float(np.mean(nn))
        stats["std_nn_distance"] = float(np.std(nn))
    stats["mean_radius"] = float(np.mean(radii))
    return fills, stats
