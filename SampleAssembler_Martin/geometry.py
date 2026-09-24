#!/usr/bin/env python3
"""
Geometry primitives and lattice-rotation math for the Sample Assembler.

Regions emit LAMMPS ``region`` commands (always ``units box`` — coordinates are
in Å, matching the multislice convention). A hemisphere has no LAMMPS primitive,
so it is built as ``region intersect`` of a sphere and a half-space block, per
the LAMMPS docs.

A ``Fill`` couples a region with a material and per-fill options (rotation,
alloy fraction, RNG seed). The writer turns a list of Fills into the .lmp body.

Rotation: a lattice is reoriented by rotating its a1/a2/a3 vectors with a 3x3
matrix (see materials.Material.lattice_line). This module builds that matrix
from Euler angles, axis-angle, a uniform-random SO(3) sample (quaternion method,
matching the original lammps_TEM.py), or three orthogonal cubic Miller
directions.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np


# ═══════════════════════════════════════════════════════════════════════════════
#  Rotation matrices  (return 3x3, indexable as R[i][j])
# ═══════════════════════════════════════════════════════════════════════════════
def rotation_matrix_euler(rx_deg, ry_deg, rz_deg):
    """Euler rotation Rz@Ry@Rx (degrees) — extrinsic X→Y→Z (≡ intrinsic Z→Y→X)."""
    rx, ry, rz = map(math.radians, (rx_deg, ry_deg, rz_deg))
    cx, sx = math.cos(rx), math.sin(rx)
    cy, sy = math.cos(ry), math.sin(ry)
    cz, sz = math.cos(rz), math.sin(rz)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def rotation_matrix_axis_angle(axis, angle_deg):
    """Rodrigues rotation about an arbitrary axis (degrees)."""
    ax = np.asarray(axis, dtype=float)
    n = np.linalg.norm(ax)
    if n == 0:
        return np.eye(3)
    ax = ax / n
    th = math.radians(angle_deg)
    c, s = math.cos(th), math.sin(th)
    x, y, z = ax
    return np.array([
        [c + x*x*(1-c),   x*y*(1-c)-z*s, x*z*(1-c)+y*s],
        [y*x*(1-c)+z*s,   c + y*y*(1-c), y*z*(1-c)-x*s],
        [z*x*(1-c)-y*s,   z*y*(1-c)+x*s, c + z*z*(1-c)],
    ])


def random_rotation_matrix(rng=None):
    """Uniform random SO(3) via quaternions (same scheme as lammps_TEM.py)."""
    rng = rng if rng is not None else np.random.default_rng()
    u1, u2, u3 = rng.random(3)
    q1 = math.sqrt(1 - u1) * math.sin(2 * math.pi * u2)
    q2 = math.sqrt(1 - u1) * math.cos(2 * math.pi * u2)
    q3 = math.sqrt(u1) * math.sin(2 * math.pi * u3)
    q4 = math.sqrt(u1) * math.cos(2 * math.pi * u3)
    return np.array([
        [1 - 2*(q3*q3 + q4*q4), 2*(q2*q3 - q1*q4),     2*(q2*q4 + q1*q3)],
        [2*(q2*q3 + q1*q4),     1 - 2*(q2*q2 + q4*q4), 2*(q3*q4 - q1*q2)],
        [2*(q2*q4 - q1*q3),     2*(q3*q4 + q1*q2),     1 - 2*(q2*q2 + q3*q3)],
    ])


def rotation_from_miller(x_hkl, y_hkl, z_hkl):
    """Rotation aligning cubic directions [x_hkl]->box-x etc. (cubic only).

    The three Miller directions must be mutually orthogonal and form a
    right-handed set, otherwise R is not a proper rotation and would shear/
    reflect the lattice — so we validate rather than silently misbuild.
    """
    ux = np.asarray(x_hkl, float); ux = ux / np.linalg.norm(ux)
    uy = np.asarray(y_hkl, float); uy = uy / np.linalg.norm(uy)
    uz = np.asarray(z_hkl, float); uz = uz / np.linalg.norm(uz)
    if (abs(ux @ uy) > 1e-6 or abs(ux @ uz) > 1e-6 or abs(uy @ uz) > 1e-6):
        raise ValueError("Miller orient directions must be mutually orthogonal "
                         "(e.g. x=[1 -1 0] y=[1 1 -2] z=[1 1 1]).")
    R = np.array([ux, uy, uz])
    if np.linalg.det(R) < 0:
        raise ValueError("Miller orient triplet is left-handed; swap two axes "
                         "to make it right-handed.")
    return R


def resolve_rotation(spec: dict, rng=None):
    """Turn a rotation spec dict into a 3x3 matrix, or None for identity."""
    if not spec:
        return None
    mode = spec.get("mode", "none")
    if mode == "none":
        return None
    if mode == "euler":
        return rotation_matrix_euler(spec.get("rx", 0), spec.get("ry", 0), spec.get("rz", 0))
    if mode == "axis":
        return rotation_matrix_axis_angle(spec.get("axis", (0, 0, 1)), spec.get("angle", 0))
    if mode == "random":
        return random_rotation_matrix(rng)
    if mode == "miller":
        return rotation_from_miller(spec.get("x", (1, 0, 0)),
                                    spec.get("y", (0, 1, 0)),
                                    spec.get("z", (0, 0, 1)))
    if mode == "matrix":
        R = spec.get("R")
        return np.asarray(R, dtype=float) if R is not None else None
    return None


# ═══════════════════════════════════════════════════════════════════════════════
#  Region primitives
# ═══════════════════════════════════════════════════════════════════════════════
@dataclass
class Region:
    """Base region. Subclasses define geometry + LAMMPS emission."""
    def lammps_regions(self, rid: str):
        """Return (list_of_command_lines, fill_region_id)."""
        raise NotImplementedError

    def aabb(self):
        """Axis-aligned bounding box (xlo,xhi,ylo,yhi,zlo,zhi)."""
        raise NotImplementedError

    def describe(self) -> str:
        raise NotImplementedError


class ImportedRegion(Region):
    """A build123d solid selected from an imported STEP file.

    LAMMPS has no native arbitrary CAD-solid region, so script generation uses
    the solid's axis-aligned bounding box while the exact solid remains
    available to the 3D preview.
    """

    def __init__(self, name, shape):
        self.name = name
        self.shape = shape

    def lammps_regions(self, rid):
        xlo, xhi, ylo, yhi, zlo, zhi = self.aabb()
        return ([f"region {rid} block {xlo:g} {xhi:g} "
                 f"{ylo:g} {yhi:g} {zlo:g} {zhi:g} units box"], rid)

    def aabb(self):
        bounds = self.shape.bounding_box()
        return (bounds.min.X, bounds.max.X, bounds.min.Y, bounds.max.Y,
                bounds.min.Z, bounds.max.Z)

    def describe(self):
        return f"STEP component '{self.name}' (bounding-box fill)"


@dataclass
class Box(Region):
    xlo: float; xhi: float; ylo: float; yhi: float; zlo: float; zhi: float

    def lammps_regions(self, rid):
        return ([f"region {rid} block {self.xlo:g} {self.xhi:g} "
                 f"{self.ylo:g} {self.yhi:g} {self.zlo:g} {self.zhi:g} units box"], rid)

    def aabb(self):
        return (self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi)

    def describe(self):
        return (f"Box  [{self.xlo:g},{self.xhi:g}]×[{self.ylo:g},{self.yhi:g}]"
                f"×[{self.zlo:g},{self.zhi:g}] Å")


@dataclass
class Sphere(Region):
    cx: float; cy: float; cz: float; r: float
    # Optional XY clip window. When all four are set the sphere is intersected
    # with this footprint over the full z range, so a particle straddling the box
    # edge is cut flush instead of protruding — this keeps the sample's XY extent
    # equal to the substrate (used by the crossgrating preset). None on any of the
    # four ⇒ no clipping (a plain sphere). z is never clipped, so protrusion above
    # the surface is preserved.
    clip_xlo: Optional[float] = None
    clip_xhi: Optional[float] = None
    clip_ylo: Optional[float] = None
    clip_yhi: Optional[float] = None

    def clip_bounds(self):
        """Return (xlo, xhi, ylo, yhi) if XY-clipped, else None."""
        b = (self.clip_xlo, self.clip_xhi, self.clip_ylo, self.clip_yhi)
        return b if None not in b else None

    def lammps_regions(self, rid):
        sphere = f"sphere {self.cx:g} {self.cy:g} {self.cz:g} {self.r:g} units box"
        clip = self.clip_bounds()
        if clip is None:
            return ([f"region {rid} {sphere}"], rid)
        # sphere ∩ XY block (INF z keeps the full vertical extent), per LAMMPS docs.
        sid, bid = f"{rid}_s", f"{rid}_b"
        cxlo, cxhi, cylo, cyhi = clip
        return ([
            f"region {sid} {sphere}",
            f"region {bid} block {cxlo:g} {cxhi:g} {cylo:g} {cyhi:g} INF INF units box",
            f"region {rid} intersect 2 {sid} {bid}",
        ], rid)

    def aabb(self):
        xlo, xhi = self.cx - self.r, self.cx + self.r
        ylo, yhi = self.cy - self.r, self.cy + self.r
        clip = self.clip_bounds()
        if clip is not None:
            cxlo, cxhi, cylo, cyhi = clip
            xlo, xhi = max(xlo, cxlo), min(xhi, cxhi)
            ylo, yhi = max(ylo, cylo), min(yhi, cyhi)
        return (xlo, xhi, ylo, yhi, self.cz - self.r, self.cz + self.r)

    def describe(self):
        clip = "" if self.clip_bounds() is None else "  (XY-clipped)"
        return f"Sphere  c=({self.cx:g},{self.cy:g},{self.cz:g}) r={self.r:g} Å{clip}"


@dataclass
class Hemisphere(Region):
    cx: float; cy: float; cz: float; r: float
    axis: str = "z"          # flat-face normal
    direction: str = "+"     # keep the '+' or '-' half along axis

    def lammps_regions(self, rid):
        sid, bid = f"{rid}_s", f"{rid}_b"
        lines = [f"region {sid} sphere {self.cx:g} {self.cy:g} {self.cz:g} {self.r:g} units box"]
        xlo, xhi, ylo, yhi, zlo, zhi = (self.cx - self.r, self.cx + self.r,
                                        self.cy - self.r, self.cy + self.r,
                                        self.cz - self.r, self.cz + self.r)
        c = {"x": self.cx, "y": self.cy, "z": self.cz}[self.axis]
        if self.axis == "x":
            xlo, xhi = (c, xhi) if self.direction == "+" else (xlo, c)
        elif self.axis == "y":
            ylo, yhi = (c, yhi) if self.direction == "+" else (ylo, c)
        else:
            zlo, zhi = (c, zhi) if self.direction == "+" else (zlo, c)
        lines.append(f"region {bid} block {xlo:g} {xhi:g} {ylo:g} {yhi:g} "
                     f"{zlo:g} {zhi:g} units box")
        lines.append(f"region {rid} intersect 2 {sid} {bid}")
        return (lines, rid)

    def aabb(self):
        xlo, xhi, ylo, yhi, zlo, zhi = (self.cx - self.r, self.cx + self.r,
                                        self.cy - self.r, self.cy + self.r,
                                        self.cz - self.r, self.cz + self.r)
        c = {"x": self.cx, "y": self.cy, "z": self.cz}[self.axis]
        if self.axis == "x":
            xlo, xhi = (c, xhi) if self.direction == "+" else (xlo, c)
        elif self.axis == "y":
            ylo, yhi = (c, yhi) if self.direction == "+" else (ylo, c)
        else:
            zlo, zhi = (c, zhi) if self.direction == "+" else (zlo, c)
        return (xlo, xhi, ylo, yhi, zlo, zhi)

    def describe(self):
        return (f"Hemisphere  c=({self.cx:g},{self.cy:g},{self.cz:g}) r={self.r:g} "
                f"flat⊥{self.axis}({self.direction})")


@dataclass
class Cylinder(Region):
    axis: str; c1: float; c2: float; r: float; lo: float; hi: float

    def lammps_regions(self, rid):
        return ([f"region {rid} cylinder {self.axis} {self.c1:g} {self.c2:g} "
                 f"{self.r:g} {self.lo:g} {self.hi:g} units box"], rid)

    def aabb(self):
        if self.axis == "x":
            return (self.lo, self.hi, self.c1 - self.r, self.c1 + self.r,
                    self.c2 - self.r, self.c2 + self.r)
        if self.axis == "y":
            return (self.c1 - self.r, self.c1 + self.r, self.lo, self.hi,
                    self.c2 - self.r, self.c2 + self.r)
        return (self.c1 - self.r, self.c1 + self.r, self.c2 - self.r,
                self.c2 + self.r, self.lo, self.hi)

    def describe(self):
        return (f"Cylinder  {self.axis}-axis c=({self.c1:g},{self.c2:g}) "
                f"r={self.r:g} {self.axis}:[{self.lo:g},{self.hi:g}]")


@dataclass
class Cone(Region):
    axis: str; c1: float; c2: float; radlo: float; radhi: float; lo: float; hi: float

    def lammps_regions(self, rid):
        return ([f"region {rid} cone {self.axis} {self.c1:g} {self.c2:g} "
                 f"{self.radlo:g} {self.radhi:g} {self.lo:g} {self.hi:g} units box"], rid)

    def aabb(self):
        rmax = max(self.radlo, self.radhi)
        if self.axis == "x":
            return (self.lo, self.hi, self.c1 - rmax, self.c1 + rmax,
                    self.c2 - rmax, self.c2 + rmax)
        if self.axis == "y":
            return (self.c1 - rmax, self.c1 + rmax, self.lo, self.hi,
                    self.c2 - rmax, self.c2 + rmax)
        return (self.c1 - rmax, self.c1 + rmax, self.c2 - rmax,
                self.c2 + rmax, self.lo, self.hi)

    def describe(self):
        return (f"Cone  {self.axis}-axis c=({self.c1:g},{self.c2:g}) "
                f"r:{self.radlo:g}->{self.radhi:g} {self.axis}:[{self.lo:g},{self.hi:g}]")


# ═══════════════════════════════════════════════════════════════════════════════
#  Fill = region + material + options
# ═══════════════════════════════════════════════════════════════════════════════
@dataclass
class Fill:
    name: str
    region: Region
    material_key: str
    rotation: dict = field(default_factory=lambda: {"mode": "none"})
    alloy_x: Optional[float] = None     # override AlloyInfo.default_x
    seed: int = 12345
    enabled: bool = True
    carve: bool = False     # delete pre-existing atoms in this footprint first

    def aabb(self):
        return self.region.aabb()


def union_aabb(fills):
    """Bounding box that contains every enabled fill (for create_box sizing)."""
    boxes = [f.aabb() for f in fills if f.enabled]
    if not boxes:
        return (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    xlo = min(b[0] for b in boxes); xhi = max(b[1] for b in boxes)
    ylo = min(b[2] for b in boxes); yhi = max(b[3] for b in boxes)
    zlo = min(b[4] for b in boxes); zhi = max(b[5] for b in boxes)
    return (xlo, xhi, ylo, yhi, zlo, zhi)
