#!/usr/bin/env python3
"""
Lightweight matplotlib 3D preview of the assembled regions.

Draws each Fill's region as a translucent/wireframe primitive colored by the
material category. When the sample has many spherical particles (e.g. the
crossgrating preset's ~140 nanoparticles) the spheres are rendered as a single
scatter of centers sized by radius instead of full wireframes, to stay fast.
"""
from __future__ import annotations

import numpy as np
import matplotlib.colors as mcolors
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from geometry import Box, Sphere, Hemisphere, Cylinder, Cone
from materials import MATERIALS
import theme

# Bright, saturated category colors chosen to read clearly on the dark 3D panes.
_CAT_COLOR = {
    "semiconductor": "#5cc0ff",
    "dielectric":    "#7fe08a",
    "metal":         "#ffae66",
    "2D":            "#cf9bff",
    "amorphous":     "#d0d0d0",
}
_SPHERE_WIREFRAME_LIMIT = 40


class PreviewCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(5, 4), tight_layout=True)
        super().__init__(self.fig)
        self.setParent(parent)
        self.ax = self.fig.add_subplot(111, projection="3d")
        self._init_axes()

    def _init_axes(self):
        theme.apply_mpl_style(self.fig, self.ax)
        t = theme.get_theme()
        # Darken the three 3D "walls" (panes) — set_facecolor alone does NOT touch
        # them, so geometry was washed out against matplotlib's default light panes.
        pane_rgba = mcolors.to_rgba(t["rule_hard"])     # near-black, high contrast
        grid_color = t["rule_soft"]
        for axis in (self.ax.xaxis, self.ax.yaxis, self.ax.zaxis):
            try:
                axis.set_pane_color(pane_rgba)
            except Exception:
                try:
                    axis.pane.set_facecolor(pane_rgba)
                    axis.pane.set_alpha(1.0)
                except Exception:
                    pass
            try:
                axis.pane.set_edgecolor(grid_color)
            except Exception:
                pass
            try:
                axis._axinfo["grid"]["color"] = grid_color
                axis._axinfo["grid"]["linewidth"] = 0.4
            except Exception:
                pass
        self.ax.set_xlabel("x (Å)")
        self.ax.set_ylabel("y (Å)")
        self.ax.set_zlabel("z (Å)")
        try:
            self.ax.set_box_aspect((1, 1, 1))
        except Exception:
            pass

    def render(self, fills, box, cad_components=None, selected_cad_name=None):
        self.ax.clear()
        self._init_axes()
        enabled = [f for f in fills if f.enabled]
        cad_components = cad_components or {}

        for name, shape in cad_components.items():
            self._draw_cad_shape(
                shape,
                selected=(name == selected_cad_name),
            )

        sphere_fills = [f for f in enabled if isinstance(f.region, Sphere)]
        many_spheres = len(sphere_fills) > _SPHERE_WIREFRAME_LIMIT

        for f in enabled:
            col = _CAT_COLOR.get(MATERIALS[f.material_key].category, "#cccccc")
            reg = f.region
            if isinstance(reg, Sphere) and many_spheres:
                continue  # batched below
            self._draw_region(reg, col)

        if many_spheres:
            xs = [f.region.cx for f in sphere_fills]
            ys = [f.region.cy for f in sphere_fills]
            zs = [f.region.cz for f in sphere_fills]
            rr = np.array([f.region.r for f in sphere_fills])
            cols = [_CAT_COLOR.get(MATERIALS[f.material_key].category, "#cccccc")
                    for f in sphere_fills]
            self.ax.scatter(xs, ys, zs, s=(rr * 0.6) ** 2, c=cols,
                            alpha=0.9, edgecolors="none", depthshade=False)

        self._draw_box_outline(box)
        limits = self._combined_limits(box, cad_components)
        self._set_limits(limits)
        self.draw()

    def _draw_cad_shape(self, shape, selected=False):
        """Draw a build123d solid from its triangulated faces."""
        vertices = []
        triangles = []
        for face in shape.faces():
            face_vertices, face_triangles = face.tessellate(0.5, 0.5)
            offset = len(vertices)
            vertices.extend((v.X, v.Y, v.Z) for v in face_vertices)
            triangles.extend(
                [(vertices[offset + a], vertices[offset + b], vertices[offset + c])
                 for a, b, c in face_triangles]
            )
        if not triangles:
            return
        color = theme.color("accent_soft") if selected else theme.color("info")
        collection = Poly3DCollection(
            triangles,
            facecolors=mcolors.to_rgba(color, 0.42 if selected else 0.16),
            edgecolors=color,
            linewidths=1.4 if selected else 0.45,
        )
        self.ax.add_collection3d(collection)

    def _combined_limits(self, box, cad_components):
        limits = list(box)
        for shape in cad_components.values():
            bounds = shape.bounding_box()
            limits[0] = min(limits[0], bounds.min.X)
            limits[1] = max(limits[1], bounds.max.X)
            limits[2] = min(limits[2], bounds.min.Y)
            limits[3] = max(limits[3], bounds.max.Y)
            limits[4] = min(limits[4], bounds.min.Z)
            limits[5] = max(limits[5], bounds.max.Z)
        return tuple(limits)

    # ── primitive drawing ──────────────────────────────────────────────────────
    def _draw_region(self, reg, col):
        if isinstance(reg, Box):
            self._wire_box(reg.xlo, reg.xhi, reg.ylo, reg.yhi, reg.zlo, reg.zhi,
                           col, face_alpha=0.20, edge_lw=1.4)
        elif isinstance(reg, Sphere):
            self._wire_sphere(reg.cx, reg.cy, reg.cz, reg.r, col,
                              clip=reg.clip_bounds())
        elif isinstance(reg, Hemisphere):
            self._wire_sphere(reg.cx, reg.cy, reg.cz, reg.r, col,
                              hemi_axis=reg.axis, hemi_dir=reg.direction)
        elif isinstance(reg, Cylinder):
            self._wire_cylinder(reg, col)
        elif isinstance(reg, Cone):
            self._wire_cone(reg, col)

    def _wire_box(self, xlo, xhi, ylo, yhi, zlo, zhi, col, face_alpha=0.0, edge_lw=1.4):
        # face translucency and edge opacity are set independently — a single
        # collection-wide alpha would also fade the edges (and hide outlines).
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        v = np.array([[xlo, ylo, zlo], [xhi, ylo, zlo], [xhi, yhi, zlo], [xlo, yhi, zlo],
                      [xlo, ylo, zhi], [xhi, ylo, zhi], [xhi, yhi, zhi], [xlo, yhi, zhi]])
        faces = [[v[0], v[1], v[2], v[3]], [v[4], v[5], v[6], v[7]],
                 [v[0], v[1], v[5], v[4]], [v[2], v[3], v[7], v[6]],
                 [v[1], v[2], v[6], v[5]], [v[0], v[3], v[7], v[4]]]
        pc = Poly3DCollection(faces, facecolors=mcolors.to_rgba(col, face_alpha),
                              edgecolors=col, linewidths=edge_lw)
        self.ax.add_collection3d(pc)

    def _wire_sphere(self, cx, cy, cz, r, col, hemi_axis=None, hemi_dir="+", clip=None):
        u = np.linspace(0, 2 * np.pi, 16)
        v = np.linspace(0, np.pi, 12)
        x = cx + r * np.outer(np.cos(u), np.sin(v))
        y = cy + r * np.outer(np.sin(u), np.sin(v))
        z = cz + r * np.outer(np.ones_like(u), np.cos(v))
        if hemi_axis is not None:
            coord = {"x": (x, cx), "y": (y, cy), "z": (z, cz)}[hemi_axis]
            arr, c = coord
            mask = (arr >= c) if hemi_dir == "+" else (arr <= c)
            x = np.where(mask, x, np.nan); y = np.where(mask, y, np.nan)
            z = np.where(mask, z, np.nan)
        if clip is not None:  # XY-clipped sphere — cut the wireframe at the footprint
            cxlo, cxhi, cylo, cyhi = clip
            m = (x >= cxlo) & (x <= cxhi) & (y >= cylo) & (y <= cyhi)
            x = np.where(m, x, np.nan); y = np.where(m, y, np.nan)
            z = np.where(m, z, np.nan)
        self.ax.plot_wireframe(x, y, z, color=col, linewidth=0.8, alpha=0.9)

    def _wire_cylinder(self, reg, col):
        t = np.linspace(0, 2 * np.pi, 20)
        h = np.linspace(reg.lo, reg.hi, 2)
        T, H = np.meshgrid(t, h)
        a, b = reg.r * np.cos(T), reg.r * np.sin(T)
        if reg.axis == "z":
            X, Y, Z = reg.c1 + a, reg.c2 + b, H
        elif reg.axis == "y":
            X, Y, Z = reg.c1 + a, H, reg.c2 + b
        else:
            X, Y, Z = H, reg.c1 + a, reg.c2 + b
        self.ax.plot_wireframe(X, Y, Z, color=col, linewidth=0.8, alpha=0.9)

    def _wire_cone(self, reg, col):
        t = np.linspace(0, 2 * np.pi, 20)
        h = np.linspace(0, 1, 6)
        T, S = np.meshgrid(t, h)
        rad = reg.radlo + (reg.radhi - reg.radlo) * S
        axpos = reg.lo + (reg.hi - reg.lo) * S
        a, b = rad * np.cos(T), rad * np.sin(T)
        if reg.axis == "z":
            X, Y, Z = reg.c1 + a, reg.c2 + b, axpos
        elif reg.axis == "y":
            X, Y, Z = reg.c1 + a, axpos, reg.c2 + b
        else:
            X, Y, Z = axpos, reg.c1 + a, reg.c2 + b
        self.ax.plot_wireframe(X, Y, Z, color=col, linewidth=0.8, alpha=0.9)

    def _draw_box_outline(self, box):
        self._wire_box(*box, col=theme.color("accent"), face_alpha=0.0, edge_lw=1.0)

    def _set_limits(self, box):
        xlo, xhi, ylo, yhi, zlo, zhi = box
        self.ax.set_xlim(xlo, xhi)
        self.ax.set_ylim(ylo, yhi)
        self.ax.set_zlim(zlo, zhi)
        # equal aspect from data ranges
        try:
            self.ax.set_box_aspect((max(xhi - xlo, 1e-6),
                                    max(yhi - ylo, 1e-6),
                                    max(zhi - zlo, 1e-6)))
        except Exception:
            pass
