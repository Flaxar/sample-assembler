#!/usr/bin/env python3
"""
Sample Assembler — a compact CAD-style tool for designing TEM specimens and
exporting ASE-generated LAMMPS data files.

Define arbitrary geometry (boxes, spheres, hemispheres, cylinders, cones), fill
each region with a material from the crystallography library (with optional
lattice rotation and alloy fraction), carve / delete overlaps, then export a
LAMMPS data file. A one-click "Crossgrating" preset reproduces the classic
gold-on-amorphous-carbon TEM resolution test specimen.

  Modules:  theme · materials · geometry · lammps_writer · crossgrating · preview3d

Run:  python sample_assembler.py
Smoke test (offscreen):
  QT_QPA_PLATFORM=offscreen python -c "from sample_assembler import MainWindow; \
      from PyQt5.QtWidgets import QApplication; QApplication([]); MainWindow()"
"""
from __future__ import annotations

import os
import sys
import glob
import json
import shutil
import subprocess
import importlib.util
from dataclasses import asdict

from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QComboBox, QDoubleSpinBox, QSpinBox, QLineEdit,
    QCheckBox, QGroupBox, QStackedWidget, QListWidget, QListWidgetItem,
    QPlainTextEdit, QTabWidget, QSplitter, QScrollArea, QFileDialog, QMessageBox,
    QDialog, QDialogButtonBox, QFormLayout, QRadioButton, QButtonGroup, QProgressBar,
)

import theme
import geometry as geo
from geometry import Box, Sphere, Hemisphere, Cylinder, Cone, ImportedRegion, Fill
from materials import MATERIALS, get_material, materials_by_category
import lammps_writer as lw
from lammps_writer import WriterOptions, MDOptions
from ase_writer import export_lammps, export_xyz, write_minimization_input, relax_with_mace
from crossgrating import CrossgratingParams, build_crossgrating
from preview3d import PreviewCanvas
from step_import import import_step_components

try:
    from potentials.registry import POTENTIALS, potentials_covering
except Exception:  # pragma: no cover - library optional
    POTENTIALS = {}

    def potentials_covering(elements):
        return []

APP_NAME = "Sample Assembler"

_SHAPE_ORDER = ["Box", "Sphere", "Hemisphere", "Cylinder", "Cone"]
_REGION_CLASSES = {"Box": Box, "Sphere": Sphere, "Hemisphere": Hemisphere,
                   "Cylinder": Cylinder, "Cone": Cone}

# field specs per shape: (key, default).  axis/direction keys become combos.
_SHAPE_FIELDS = {
    "Box":        [("xlo", 0), ("xhi", 200), ("ylo", 0), ("yhi", 200), ("zlo", 0), ("zhi", 200)],
    "Sphere":     [("cx", 100), ("cy", 100), ("cz", 100), ("r", 40)],
    "Hemisphere": [("cx", 100), ("cy", 100), ("cz", 100), ("r", 40),
                   ("axis", "z"), ("direction", "+")],
    "Cylinder":   [("axis", "z"), ("c1", 100), ("c2", 100), ("r", 40), ("lo", 0), ("hi", 200)],
    "Cone":       [("axis", "z"), ("c1", 100), ("c2", 100), ("radlo", 40), ("radhi", 10),
                   ("lo", 0), ("hi", 200)],
}


def find_lammps():
    for exe in ("lmp.exe", "lammps.exe", "lmp_serial.exe", "lmp_mpi.exe",
                "lmp_win_no-mpi.exe", "lmp", "lammps", "lmp_serial", "lmp_mpi"):
        path = shutil.which(exe)
        if path:
            return path
    return None


def mace_available():
    if importlib.util.find_spec("mace") is None:
        return False
    try:
        probe = subprocess.run(
            [sys.executable, "-c",
             "import torch; from mace.calculators import mace_mp"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            timeout=30, check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return probe.returncode == 0


# ── helpers ─────────────────────────────────────────────────────────────────────
def _dspin(value=0.0, lo=-1e6, hi=1e6, step=10.0, decimals=3):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(decimals)
    s.setSingleStep(step)
    s.setValue(value)
    return s


def make_button(label, kind=None):
    b = QPushButton(label)
    if kind:
        b.setProperty("kind", kind)
    return b


# ═══════════════════════════════════════════════════════════════════════════════
#  LAMMPS runner thread
# ═══════════════════════════════════════════════════════════════════════════════
class RunnerThread(QThread):
    output_line = pyqtSignal(str)         # one line of LAMMPS stdout/stderr (live)
    finished_run = pyqtSignal(int)        # returncode

    def __init__(self, exe, script_path):
        super().__init__()
        self.exe = exe
        self.script_path = script_path
        self._proc = None

    def run(self):
        try:
            self._proc = subprocess.Popen(
                [self.exe, "-in", os.path.basename(self.script_path)],
                cwd=os.path.dirname(os.path.abspath(self.script_path)),
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1,           # line-buffered streaming
            )
            for line in self._proc.stdout:      # emit each line as LAMMPS prints it
                self.output_line.emit(line.rstrip("\n"))
            self._proc.wait()
            self.finished_run.emit(self._proc.returncode)
        except Exception as exc:  # noqa: BLE001
            self.output_line.emit(f"Run failed: {exc}")
            self.finished_run.emit(-1)

    def stop(self):
        if self._proc and self._proc.poll() is None:
            self._proc.terminate()


class MaceRunnerThread(QThread):
    output_line = pyqtSignal(str)
    finished_relaxation = pyqtSignal(int, str)

    def __init__(self, input_path, output_path, model, fmax):
        super().__init__()
        self.input_path = input_path
        self.output_path = output_path
        self.model = model
        self.fmax = fmax

    def run(self):
        try:
            import json
            code = (
                "from ase_writer import relax_with_mace; "
                "relax_with_mace("
                f"{json.dumps(self.input_path)}, {json.dumps(self.output_path)}, "
                f"model={json.dumps(self.model)}, fmax={self.fmax!r}, "
                f"trajectory_path={json.dumps(os.path.splitext(self.output_path)[0] + '.traj')})"
            )
            self._proc = subprocess.Popen(
                [sys.executable, "-c", code],
                cwd=os.path.dirname(os.path.abspath(__file__)),
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1,
            )
            for line in self._proc.stdout:
                self.output_line.emit(line.rstrip())
            returncode = self._proc.wait()
            if returncode:
                raise RuntimeError(f"MACE process exited with {returncode}")
        except Exception as exc:  # noqa: BLE001
            self.finished_relaxation.emit(-1, str(exc))
        else:
            self.finished_relaxation.emit(0, self.output_path)

    def stop(self):
        if getattr(self, "_proc", None) and self._proc.poll() is None:
            self._proc.terminate()


# ═══════════════════════════════════════════════════════════════════════════════
#  Crossgrating preset dialog
# ═══════════════════════════════════════════════════════════════════════════════
class CrossgratingDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Crossgrating — gold on amorphous carbon")
        self.setModal(True)
        form = QFormLayout(self)
        self.f = {}
        specs = [
            ("box_x", "Box X (Å)", 400, 0), ("box_y", "Box Y (Å)", 400, 0),
            ("carbon_thickness", "a-C thickness (Å)", 80, 0),
            ("n_particles", "Number of Au NPs", 40, 1),
            ("r_min", "Min NP radius (Å)", 15, 0), ("r_max", "Max NP radius (Å)", 40, 0),
            ("cluster_strength", "Cluster strength (0–1)", 0.6, 2),
            ("overlap_factor", "Overlap factor (<1 overlaps)", 0.4, 2),
            ("protrusion", "Protrusion above surface (0–1)", 0.5, 2),
            ("seed", "Random seed", 42, 1),
        ]
        for key, label, default, dec in specs:
            if dec == 0 or dec == 1:
                w = QSpinBox(); w.setRange(0, 10_000_000); w.setValue(int(default))
                if dec == 1 and key not in ("seed", "n_particles"):
                    pass
            else:
                w = QDoubleSpinBox(); w.setRange(0, 1e6); w.setDecimals(dec)
                w.setSingleStep(0.05 if default < 2 else 10); w.setValue(float(default))
            self.f[key] = w
            form.addRow(label, w)
        self.carve = QCheckBox("Gold carves the carbon beneath it"); self.carve.setChecked(True)
        form.addRow(self.carve)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        form.addRow(bb)

    def params(self):
        g = lambda k: self.f[k].value()
        return CrossgratingParams(
            box_x=float(g("box_x")), box_y=float(g("box_y")),
            carbon_thickness=float(g("carbon_thickness")),
            n_particles=int(g("n_particles")),
            r_min=float(g("r_min")), r_max=float(g("r_max")),
            cluster_strength=float(g("cluster_strength")),
            overlap_factor=float(g("overlap_factor")),
            protrusion=float(g("protrusion")), seed=int(g("seed")),
            carve_carbon=self.carve.isChecked(),
        )


# ═══════════════════════════════════════════════════════════════════════════════
#  Main window
# ═══════════════════════════════════════════════════════════════════════════════
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.resize(1320, 860)
        self.fills = []
        self.cad_components = {}
        self.cad_file_path = None
        self.selected_cad_component = None
        self.opts = WriterOptions()
        self.lammps_exe = find_lammps()
        self._last_saved_path = None
        self._last_mace_input = None
        self._last_mace_output = None
        self._last_data_path = None
        self._runner = None
        # Embedding hook (SEB-MS "Specimens" workspace): host sets a callback
        # via set_transfer_callback(cb); cb(dump_paths, elements) receives the
        # last run's .xyz path(s) + mass-sorted element list. Standalone runs
        # never show the Transfer button.
        self.transfer_callback = None
        self._last_dump_paths = []
        self._last_type_table = {}

        self._build_ui()
        self._refresh()

    # ── UI construction ─────────────────────────────────────────────────────────
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(10, 10, 10, 10)
        root.setSpacing(0)

        # left | right are a draggable horizontal splitter so the controls panel
        # is resizable (no fixed width clipping labels under the right pane).
        self.main_splitter = QSplitter(Qt.Horizontal)
        root.addWidget(self.main_splitter)

        # left: scrollable controls (resizable, with a horizontal scrollbar if needed)
        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setMinimumWidth(360)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        left_inner = QWidget()
        left_inner.setMinimumWidth(400)        # content's natural width; scrollbar kicks in below this
        self.left = QVBoxLayout(left_inner)
        self.left.setSpacing(10)
        left_scroll.setWidget(left_inner)
        self.main_splitter.addWidget(left_scroll)

        self._build_add_region_group()
        self._build_fills_group()
        self._build_options_group()
        self._build_md_group()
        self._build_actions_group()
        self.left.addStretch(1)

        # Growing-library note — the materials list and presets will expand over time.
        footer = QLabel("⚙  This tool is evolving: more materials, elements, and sample "
                        "presets are planned. The library is built to extend easily.")
        footer.setProperty("role", "hint"); footer.setWordWrap(True)
        self.left.addWidget(footer)

        # right: preview + tabs
        right = QSplitter(Qt.Vertical)
        self.preview = PreviewCanvas()
        right.addWidget(self.preview)

        self.tabs = QTabWidget()
        self.script_view = QPlainTextEdit(); self.script_view.setReadOnly(True)
        self.script_view.setFont(QFont("Consolas", 9))
        self.info_view = QPlainTextEdit(); self.info_view.setReadOnly(True)
        self.info_view.setFont(QFont("Consolas", 9))
        info_page = QWidget()
        info_layout = QVBoxLayout(info_page)
        info_layout.setContentsMargins(0, 0, 0, 0)
        info_layout.addWidget(self.info_view)
        self.export_progress = QProgressBar()
        self.export_progress.setRange(0, 100)
        self.export_progress.setValue(0)
        self.export_progress.setFormat("Ready")
        info_layout.addWidget(self.export_progress)
        self.console_view = QPlainTextEdit(); self.console_view.setReadOnly(True)
        self.console_view.setFont(QFont("Consolas", 9))
        self.console_view.setMaximumBlockCount(50000)   # bound memory on long runs
        self.console_view.setPlaceholderText("LAMMPS stdout/stderr appears here when you Run.")
        self.tabs.addTab(self.script_view, "Generated .lmp")
        self.tabs.addTab(info_page, "Info & Log")
        self.tabs.addTab(self.console_view, "LAMMPS Console")
        right.addWidget(self.tabs)
        right.setSizes([520, 320])
        self.main_splitter.addWidget(right)

        self.main_splitter.setStretchFactor(0, 0)   # left keeps its size
        self.main_splitter.setStretchFactor(1, 1)   # right grows with the window
        self.main_splitter.setSizes([440, 880])

        self._build_menu()
        self.setStyleSheet(theme.build_stylesheet())
        from PyQt5.QtWidgets import QStyleFactory
        QApplication.setStyle(QStyleFactory.create("Fusion"))

    def _build_add_region_group(self):
        g = QGroupBox("1 · Add region + material")
        v = QVBoxLayout(g)

        # Material and alloy controls apply to either region source.
        mat_row = QHBoxLayout()
        mat_row.addWidget(QLabel("Material"))
        self.material_combo = QComboBox()
        self._populate_materials()
        self.material_combo.currentIndexChanged.connect(self._on_material_changed)
        mat_row.addWidget(self.material_combo, 1)
        v.addLayout(mat_row)

        self.material_note = QLabel("")
        self.material_note.setProperty("role", "hint")
        self.material_note.setWordWrap(True)
        v.addWidget(self.material_note)

        # alloy fraction (shown only for alloys)
        self.alloy_row = QWidget()
        ar = QHBoxLayout(self.alloy_row); ar.setContentsMargins(0, 0, 0, 0)
        ar.addWidget(QLabel("Alloy guest fraction x"))
        self.alloy_x = _dspin(0.3, 0.0, 1.0, 0.05, 2)
        ar.addWidget(self.alloy_x, 1)
        v.addWidget(self.alloy_row)

        self.region_method_group = QButtonGroup(self)
        self.region_method_group.setExclusive(True)

        # STEP component source
        step_group = QGroupBox("Imported STEP component")
        step_layout = QVBoxLayout(step_group)
        self.step_method_button = QRadioButton("Use imported .step component")
        self.region_method_group.addButton(self.step_method_button)
        step_layout.addWidget(self.step_method_button)

        self.import_step_btn = make_button("Import .step")
        self.import_step_btn.clicked.connect(self._import_step)
        step_layout.addWidget(self.import_step_btn)

        self.cad_component_combo = QComboBox()
        self.cad_component_combo.setPlaceholderText("No STEP components imported")
        self.cad_component_combo.currentIndexChanged.connect(self._on_cad_component_changed)
        step_layout.addWidget(self.cad_component_combo)

        self.cad_component_note = QLabel("")
        self.cad_component_note.setProperty("role", "hint")
        self.cad_component_note.setWordWrap(True)
        step_layout.addWidget(self.cad_component_note)
        v.addWidget(step_group)

        # Primitive shape source
        shape_group = QGroupBox("Primitive shape")
        shape_layout = QVBoxLayout(shape_group)
        self.shape_method_button = QRadioButton("Use primitive shape")
        self.shape_method_button.setChecked(True)
        self.region_method_group.addButton(self.shape_method_button)
        shape_layout.addWidget(self.shape_method_button)

        shape_row = QHBoxLayout()
        shape_row.addWidget(QLabel("Shape"))
        self.shape_combo = QComboBox(); self.shape_combo.addItems(_SHAPE_ORDER)
        self.shape_combo.currentIndexChanged.connect(lambda i: self.shape_stack.setCurrentIndex(i))
        shape_row.addWidget(self.shape_combo, 1)
        shape_layout.addLayout(shape_row)

        self.shape_stack = QStackedWidget()
        self.shape_widgets = {}
        for shape in _SHAPE_ORDER:
            self.shape_stack.addWidget(self._build_shape_page(shape))
        shape_layout.addWidget(self.shape_stack)
        v.addWidget(shape_group)

        self.step_method_button.toggled.connect(self._update_region_method)
        self._update_region_method()

        # Lattice rotation follows material selection and region source selection.
        rot_row = QHBoxLayout()
        rot_row.addWidget(QLabel("Lattice rotation"))
        self.rot_combo = QComboBox()
        self.rot_combo.addItems(["None", "Euler (deg)", "Axis-angle", "Random (SO3)", "Miller (cubic)"])
        self.rot_combo.currentIndexChanged.connect(lambda i: self.rot_stack.setCurrentIndex(i))
        rot_row.addWidget(self.rot_combo, 1)
        v.addLayout(rot_row)
        v.addWidget(self._build_rotation_stack())

        self.carve_check = QCheckBox("Carve: delete pre-existing atoms in this footprint first")
        v.addWidget(self.carve_check)

        add_btn = make_button("Add region to sample", "primary")
        add_btn.clicked.connect(self._add_fill)
        v.addWidget(add_btn)
        self.left.addWidget(g)

    def _update_region_method(self):
        use_step = self.step_method_button.isChecked()
        self.cad_component_combo.setEnabled(use_step)
        self.cad_component_note.setEnabled(use_step)
        self.import_step_btn.setEnabled(use_step)
        self.shape_combo.setEnabled(not use_step)
        self.shape_stack.setEnabled(not use_step)

    def _import_step(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Import STEP file", "",
            "STEP files (*.step *.stp);;All files (*)",
        )
        if not path:
            return
        try:
            components = import_step_components(path)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "STEP import failed", str(exc))
            return
        if not components:
            QMessageBox.warning(self, "STEP import", "The selected file contains no solid components.")
            return

        self.cad_file_path = path
        self.cad_components = components
        self.cad_component_combo.blockSignals(True)
        self.cad_component_combo.clear()
        self.cad_component_combo.addItems(list(components))
        self.cad_component_combo.blockSignals(False)
        self.cad_component_combo.setCurrentIndex(0)
        self._on_cad_component_changed(0)
        self._log(f"Imported {len(components)} STEP component(s): {os.path.basename(path)}")
        self._refresh()

    def _on_cad_component_changed(self, index):
        if index < 0 or index >= self.cad_component_combo.count():
            self.selected_cad_component = None
            self.cad_component_note.clear()
            return
        name = self.cad_component_combo.itemText(index)
        shape = self.cad_components[name]
        self.selected_cad_component = shape
        bounds = shape.bounding_box()
        self.cad_component_note.setText(
            f"{name}\n"
            f"Bounds: X[{bounds.min.X:g}, {bounds.max.X:g}]  "
            f"Y[{bounds.min.Y:g}, {bounds.max.Y:g}]  "
            f"Z[{bounds.min.Z:g}, {bounds.max.Z:g}] Å"
        )
        if hasattr(self, "preview"):
            self._refresh()

    def _build_shape_page(self, shape):
        page = QWidget()
        grid = QGridLayout(page); grid.setContentsMargins(0, 0, 0, 0)
        widgets = {}
        for i, (key, default) in enumerate(_SHAPE_FIELDS[shape]):
            grid.addWidget(QLabel(key), i, 0)
            if key == "axis":
                w = QComboBox(); w.addItems(["x", "y", "z"]); w.setCurrentText(str(default))
            elif key == "direction":
                w = QComboBox(); w.addItems(["+", "-"]); w.setCurrentText(str(default))
            else:
                w = _dspin(float(default))
            widgets[key] = w
            grid.addWidget(w, i, 1)
        self.shape_widgets[shape] = widgets
        return page

    def _build_rotation_stack(self):
        self.rot_stack = QStackedWidget()
        # None
        self.rot_stack.addWidget(QWidget())
        # Euler
        ew = QWidget(); eg = QGridLayout(ew); eg.setContentsMargins(0, 0, 0, 0)
        self.euler = {}
        for i, k in enumerate(("rx", "ry", "rz")):
            eg.addWidget(QLabel(k), i, 0)
            self.euler[k] = _dspin(0.0, -360, 360, 5, 1)
            eg.addWidget(self.euler[k], i, 1)
        self.rot_stack.addWidget(ew)
        # Axis-angle
        aw = QWidget(); ag = QGridLayout(aw); ag.setContentsMargins(0, 0, 0, 0)
        self.axisang = {}
        for i, (k, d) in enumerate((("ax", 0), ("ay", 0), ("az", 1), ("angle", 30))):
            ag.addWidget(QLabel(k), i, 0)
            self.axisang[k] = _dspin(float(d), -360, 360, 1, 2)
            ag.addWidget(self.axisang[k], i, 1)
        self.rot_stack.addWidget(aw)
        # Random
        rw = QWidget(); rl = QVBoxLayout(rw); rl.setContentsMargins(0, 0, 0, 0)
        rl.addWidget(QLabel("Uniform random orientation (seeded per region)."))
        self.rot_stack.addWidget(rw)
        # Miller
        mw = QWidget(); mg = QGridLayout(mw); mg.setContentsMargins(0, 0, 0, 0)
        self.miller = {}
        for i, (k, d) in enumerate((("x", "1 0 0"), ("y", "0 1 0"), ("z", "0 0 1"))):
            mg.addWidget(QLabel(f"orient {k}"), i, 0)
            self.miller[k] = QLineEdit(d)
            mg.addWidget(self.miller[k], i, 1)
        self.rot_stack.addWidget(mw)
        return self.rot_stack

    def _build_fills_group(self):
        g = QGroupBox("2 · Sample regions (stacking order)")
        v = QVBoxLayout(g)
        self.fills_list = QListWidget()
        self.fills_list.itemChanged.connect(self._on_item_changed)
        v.addWidget(self.fills_list)
        row = QHBoxLayout()
        for label, slot in (("↑", self._move_up), ("↓", self._move_down),
                            ("Remove", self._remove_fill), ("Clear", self._clear_fills)):
            b = make_button(label); b.clicked.connect(slot); row.addWidget(b)
        v.addLayout(row)
        cg = make_button("Crossgrating preset (Au on a-C)…")
        cg.clicked.connect(self._open_crossgrating)
        v.addWidget(cg)
        preset_hint = QLabel("More sample presets will be added in future versions.")
        preset_hint.setProperty("role", "hint"); preset_hint.setWordWrap(True)
        v.addWidget(preset_hint)
        self.left.addWidget(g)

    def _build_options_group(self):
        g = QGroupBox("3 · Build options")
        form = QFormLayout(g)
        self.auto_box = QCheckBox("Auto box from regions"); self.auto_box.setChecked(True)
        self.auto_box.stateChanged.connect(self._refresh)
        form.addRow(self.auto_box)
        self.margin = _dspin(0.0, 0, 1e5, 5, 2); self.margin.valueChanged.connect(self._refresh)
        form.addRow("Box margin (Å)", self.margin)
        self.boundary = QComboBox(); self.boundary.addItems(["p p p", "s s s", "f f f", "p p f", "p p s"])
        self.boundary.currentIndexChanged.connect(self._refresh)
        form.addRow("Boundary", self.boundary)
        self.del_overlap = QCheckBox("Delete overlapping atoms"); self.del_overlap.setChecked(True)
        self.del_overlap.stateChanged.connect(self._refresh)
        form.addRow(self.del_overlap)
        self.overlap_cut = _dspin(1.2, 0.1, 10, 0.1, 2); self.overlap_cut.valueChanged.connect(self._refresh)
        form.addRow("Overlap cutoff (Å)", self.overlap_cut)
        self.unify_similar_lattices = QCheckBox("Use smallest spacing for similar lattices")
        self.unify_similar_lattices.setChecked(True)
        self.unify_similar_lattices.toggled.connect(self._on_lattice_unification_toggled)
        self.unify_similar_lattices.toggled.connect(self._refresh)
        form.addRow(self.unify_similar_lattices)
        self.lattice_similarity_threshold = _dspin(5.0, 0.0, 100.0, 0.5, 1)
        self.lattice_similarity_threshold.valueChanged.connect(self._refresh)
        form.addRow("Similarity threshold (%)", self.lattice_similarity_threshold)
        self._on_lattice_unification_toggled(True)
        self.min_sep = QCheckBox("Enforce min-separation for amorphous (slow for huge films)")
        self.min_sep.setChecked(True); self.min_sep.stateChanged.connect(self._refresh)
        form.addRow(self.min_sep)
        self.dump_name = QLineEdit("sample.xyz"); self.dump_name.textChanged.connect(self._refresh)
        form.addRow("Dump file", self.dump_name)
        self.seed_base = QSpinBox(); self.seed_base.setRange(1, 10_000_000); self.seed_base.setValue(12345)
        self.seed_base.valueChanged.connect(self._refresh)
        form.addRow("Seed base", self.seed_base)
        self.left.addWidget(g)

    def _on_lattice_unification_toggled(self, checked):
        self.lattice_similarity_threshold.setEnabled(checked)

    def _build_md_group(self):
        g = QGroupBox("4 · MD simulation (optional)")
        v = QVBoxLayout(g)

        self.md_enable = QCheckBox("Run LAMMPS relaxation after export")
        self.md_enable.stateChanged.connect(self._on_md_enable_toggled)
        self.md_enable.setEnabled(False)
        self.md_enable.setToolTip("Disabled: ASE creates the LAMMPS data file directly.")
        v.addWidget(self.md_enable)

        self.md_body = QWidget()
        form = QFormLayout(self.md_body); form.setContentsMargins(0, 0, 0, 0)

        self.md_backend = QComboBox()
        self.md_backend.addItem("LAMMPS", "lammps")
        self.md_backend.addItem("MACE MLIP (ASE)", "mace")
        self.md_backend.currentIndexChanged.connect(self._on_backend_changed)
        form.addRow("Relaxation backend", self.md_backend)

        self.md_potential = QComboBox()                      # populated in _update_md_potentials
        self.md_potential.currentIndexChanged.connect(self._refresh)
        form.addRow("Potential", self.md_potential)

        self.md_mace_model = QComboBox()
        self.md_mace_model.addItems(["medium", "large"])
        self.md_mace_model.currentIndexChanged.connect(self._refresh)
        form.addRow("MACE model", self.md_mace_model)
        self.md_mace_fmax = _dspin(0.05, 0.001, 1.0, 0.01, 3)
        self.md_mace_fmax.valueChanged.connect(self._refresh)
        form.addRow("MACE force threshold (eV/Å)", self.md_mace_fmax)

        self.md_hint = QLabel(""); self.md_hint.setProperty("role", "hint"); self.md_hint.setWordWrap(True)
        form.addRow(self.md_hint)

        self.md_workflow = QComboBox()
        for label, key in (("Thermalize → frozen-phonon snapshots", "thermalize"),
                           ("Relax (minimise only)", "relax"),
                           ("Melt-quench (amorphous)", "melt_quench"),
                           ("Anneal (heat · hold · cool)", "anneal")):
            self.md_workflow.addItem(label, key)
        self.md_workflow.currentIndexChanged.connect(self._on_md_workflow_changed)
        form.addRow("Workflow", self.md_workflow)

        self.md_temp = _dspin(300.0, 0.0, 1e5, 10, 1)
        self.md_temp.valueChanged.connect(self._refresh)
        form.addRow("Temperature (K)", self.md_temp)

        self.md_peakT = _dspin(6000.0, 0.0, 2e4, 100, 0)
        self.md_peakT.valueChanged.connect(self._refresh)
        self.md_peakT_label = QLabel("Melt/anneal peak T (K)")
        form.addRow(self.md_peakT_label, self.md_peakT)

        self.md_timestep = _dspin(0.0, 0.0, 10.0, 0.1, 3)
        self.md_timestep.setSpecialValueText("auto")          # 0 -> potential's recommended dt
        self.md_timestep.valueChanged.connect(self._refresh)
        form.addRow("Timestep (fs, 0=auto)", self.md_timestep)

        self.md_ensemble = QComboBox(); self.md_ensemble.addItems(["nvt", "nve", "npt"])
        self.md_ensemble.currentIndexChanged.connect(self._refresh)
        form.addRow("Ensemble", self.md_ensemble)

        self.md_nequil = QSpinBox(); self.md_nequil.setRange(0, 100_000_000); self.md_nequil.setValue(20000)
        self.md_nequil.valueChanged.connect(self._refresh)
        form.addRow("Equilibration steps", self.md_nequil)

        self.md_nsnap = QSpinBox(); self.md_nsnap.setRange(1, 100_000); self.md_nsnap.setValue(10)
        self.md_nsnap.valueChanged.connect(self._refresh)
        self.md_nsnap_label = QLabel("Snapshots")
        form.addRow(self.md_nsnap_label, self.md_nsnap)

        self.md_interval = QSpinBox(); self.md_interval.setRange(1, 10_000_000); self.md_interval.setValue(500)
        self.md_interval.valueChanged.connect(self._refresh)
        self.md_interval_label = QLabel("Snapshot interval (steps)")
        form.addRow(self.md_interval_label, self.md_interval)

        self.md_seed = QSpinBox(); self.md_seed.setRange(1, 100_000_000); self.md_seed.setValue(90210)
        self.md_seed.valueChanged.connect(self._refresh)
        form.addRow("MD seed", self.md_seed)

        v.addWidget(self.md_body)

        if not POTENTIALS:
            self.md_enable.setEnabled(False)
            self.md_enable.setToolTip("No potentials curated in potentials/ — clean build only.")
        self.md_body.setEnabled(False)
        self._update_backend_controls()
        self._apply_md_workflow_visibility()                  # no _refresh (preview not built yet)
        self.left.addWidget(g)

    def _on_md_enable_toggled(self):
        self.md_body.setEnabled(self.md_enable.isChecked())
        self._refresh()

    def _on_backend_changed(self):
        self._update_backend_controls()
        self._refresh()

    def _update_backend_controls(self):
        is_mace = self.md_backend.currentData() == "mace"
        self.md_potential.setEnabled(not is_mace)
        self.md_mace_model.setEnabled(is_mace)
        self.md_mace_fmax.setEnabled(is_mace)
        self.md_enable.setText(
            "Run MACE MLIP relaxation after export" if is_mace
            else "Run LAMMPS relaxation after export")
        if hasattr(self, "run_btn"):
            self.run_btn.setText("Run MACE relaxation" if is_mace else "Run LAMMPS")
        if is_mace and not mace_available():
            self.md_enable.setEnabled(False)
            self.md_enable.setToolTip("Install mace-torch to enable MACE relaxation.")

    def _apply_md_workflow_visibility(self):
        wf = self.md_workflow.currentData()
        is_therm = (wf == "thermalize")
        is_meltanneal = wf in ("melt_quench", "anneal")
        for w in (self.md_nsnap, self.md_nsnap_label, self.md_interval, self.md_interval_label):
            w.setVisible(is_therm)
        for w in (self.md_peakT, self.md_peakT_label):
            w.setVisible(is_meltanneal)

    def _on_md_workflow_changed(self):
        self._apply_md_workflow_visibility()
        self._refresh()

    def _update_md_potentials(self):
        """Repopulate the potential combo with potentials covering the sample's elements
        (preserving the current selection by key). Signals blocked → no recursion."""
        elems = set()
        for f in self.fills:
            if f.enabled:
                elems.update(get_material(f.material_key).elements())
        covering = potentials_covering(elems) if elems and POTENTIALS else []
        pots = covering if covering else list(POTENTIALS.values())
        backend = self.md_backend.currentData()
        self.md_enable.setEnabled(bool(covering) if backend == "lammps" else mace_available())
        if backend == "lammps" and not covering:
            self.md_enable.setChecked(False)
            self.md_enable.setToolTip("No curated potential covers all sample elements.")
        elif backend == "lammps":
            self.md_enable.setToolTip("Run a LAMMPS minimization after exporting the data file.")
        prev = self.md_potential.currentData()
        self.md_potential.blockSignals(True)
        self.md_potential.clear()
        for p in sorted(pots, key=lambda p: p.key):
            self.md_potential.addItem(p.name, p.key)
        i = self.md_potential.findData(prev)
        if i >= 0:
            self.md_potential.setCurrentIndex(i)
        self.md_potential.blockSignals(False)
        self._update_md_hint(elems)

    def _update_md_hint(self, elems):
        key = self.md_potential.currentData()
        p = POTENTIALS.get(key) if key else None
        if p is None:
            self.md_hint.setText("No curated potential covers this sample — clean build.")
            return
        conf = "" if p.confidence == "high" else f"   ⚠ {p.confidence}-confidence"
        missing = sorted(e for e in elems if not p.covers(e))
        cover = "" if not missing else f"   ⚠ does NOT cover {', '.join(missing)} (those types get NULL)"
        self.md_hint.setText(f"{p.citation}{conf}{cover}")

    def _md_options(self):
        if not self.md_enable.isChecked():
            return None
        return MDOptions(
            potential_key=self.md_potential.currentData() or "",
            workflow=self.md_workflow.currentData() or "relax",
            temperature=self.md_temp.value(),
            timestep_fs=self.md_timestep.value(),
            ensemble=self.md_ensemble.currentText(),
            equilibration_steps=self.md_nequil.value(),
            snapshots=self.md_nsnap.value(),
            snapshot_interval=self.md_interval.value(),
            seed=self.md_seed.value(),
            backend=self.md_backend.currentData() or "lammps",
            mace_model=self.md_mace_model.currentText(),
            mace_fmax=self.md_mace_fmax.value(),
        )

    def _build_actions_group(self):
        g = QGroupBox("5 · Generate / Run")
        v = QVBoxLayout(g)
        gen = make_button("Export LAMMPS data…", "primary"); gen.clicked.connect(self._generate)
        v.addWidget(gen)
        xyz = make_button("Export .xyz data…")
        xyz.clicked.connect(self._export_xyz)
        v.addWidget(xyz)
        self.run_btn = make_button("Run LAMMPS")
        self.run_btn.clicked.connect(self._run_relaxation)
        self.run_btn.setEnabled(False)
        self.run_btn.setToolTip("Export a data file with minimization enabled first.")
        v.addWidget(self.run_btn)
        self.stop_relax_btn = make_button("Stop relaxation")
        self.stop_relax_btn.clicked.connect(self._stop_relaxation)
        self.stop_relax_btn.setEnabled(False)
        v.addWidget(self.stop_relax_btn)
        choose_lammps = make_button("Choose LAMMPS executable…")
        choose_lammps.clicked.connect(self._choose_lammps)
        v.addWidget(choose_lammps)
        # Visible only when embedded in SEB-MS (set_transfer_callback)
        self.transfer_btn = make_button("Transfer to SEB-MS →", "primary")
        self.transfer_btn.clicked.connect(self._transfer_to_host)
        self.transfer_btn.setVisible(False)
        self.transfer_btn.setEnabled(False)
        self.transfer_btn.setToolTip("Load the dumped .xyz snapshot(s) into the "
                                     "SEB-MS Simulator's input files")
        v.addWidget(self.transfer_btn)
        status = "Export creates a LAMMPS data file; optional relaxation runs LAMMPS afterward."
        lbl = QLabel(status); lbl.setProperty("role", "hint"); lbl.setWordWrap(True)
        v.addWidget(lbl)
        row = QHBoxLayout()
        sp = make_button("Save project"); sp.clicked.connect(self._save_project)
        lp = make_button("Load project"); lp.clicked.connect(self._load_project)
        row.addWidget(sp); row.addWidget(lp)
        v.addLayout(row)
        self.left.addWidget(g)

    def _choose_lammps(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose LAMMPS executable", "",
            "LAMMPS executable (lmp*.exe lammps*.exe);;Executable files (*.exe);;All files (*)",
        )
        if not path:
            return
        name = os.path.basename(path).lower()
        if not (name.startswith("lmp") or name.startswith("lammps")):
            QMessageBox.warning(self, "Invalid executable",
                                "Choose the LAMMPS executable, for example lmp.exe.")
            return
        self.lammps_exe = path
        self._log(f"Selected LAMMPS executable: {path}")
        if self._last_saved_path and self._md_options() is not None:
            self.run_btn.setEnabled(True)

    def _build_menu(self):
        mb = self.menuBar()
        filem = mb.addMenu("&File")
        filem.addAction("Export LAMMPS data…", self._generate)
        filem.addAction("Save project…", self._save_project)
        filem.addAction("Load project…", self._load_project)
        filem.addSeparator()
        # kept as an attribute so embedded mode (set_transfer_callback) can
        # hide it — closing the embedded pane would just blank the workspace
        self._quit_action = filem.addAction("Quit", self.close)
        toolsm = mb.addMenu("&Tools")
        toolsm.addAction("Split multi-timestep .xyz…", self._split_xyz_tool)
        helpm = mb.addMenu("&Help")
        helpm.addAction("About", self._about)

    # ── material combo ──────────────────────────────────────────────────────────
    def _populate_materials(self):
        cats = materials_by_category()
        first = True
        for cat in ("semiconductor", "dielectric", "metal", "2D", "amorphous"):
            if cat not in cats:
                continue
            if not first:
                self.material_combo.insertSeparator(self.material_combo.count())
            first = False
            for m in cats[cat]:
                self.material_combo.addItem(f"{m.name}  ({m.formula})", m.key)

    def _current_material_key(self):
        return self.material_combo.currentData()

    def _on_material_changed(self):
        key = self._current_material_key()
        if key is None:
            return
        m = get_material(key)
        flag = "" if m.confidence == "high" else f"  ⚠ {m.confidence}-confidence"
        self.material_note.setText(m.notes + flag)
        self.alloy_row.setVisible(m.alloy is not None)
        if m.alloy is not None:
            self.alloy_x.setValue(m.alloy.default_x)

    # ── building Fills from the form ──────────────────────────────────────────────
    def _read_region(self):
        shape = self.shape_combo.currentText()
        w = self.shape_widgets[shape]
        vals = {}
        for key, _ in _SHAPE_FIELDS[shape]:
            widget = w[key]
            vals[key] = widget.currentText() if isinstance(widget, QComboBox) else widget.value()
        return _REGION_CLASSES[shape](**vals)

    def _read_rotation(self):
        mode = self.rot_combo.currentIndex()
        if mode == 1:
            return {"mode": "euler", **{k: v.value() for k, v in self.euler.items()}}
        if mode == 2:
            a = self.axisang
            return {"mode": "axis", "axis": (a["ax"].value(), a["ay"].value(), a["az"].value()),
                    "angle": a["angle"].value()}
        if mode == 3:
            return {"mode": "random"}
        if mode == 4:
            def hkl(s):
                return tuple(float(x) for x in s.split())
            return {"mode": "miller", "x": hkl(self.miller["x"].text()),
                    "y": hkl(self.miller["y"].text()), "z": hkl(self.miller["z"].text())}
        return {"mode": "none"}

    def _add_fill(self):
        key = self._current_material_key()
        mat = get_material(key)
        try:
            if self.step_method_button.isChecked():
                component_name = self.cad_component_combo.currentText()
                if not component_name or component_name not in self.cad_components:
                    raise ValueError("Import a STEP file and select a component first.")
                region = ImportedRegion(
                    component_name,
                    self.cad_components[component_name],
                )
            else:
                region = self._read_region()
            rotation = self._read_rotation()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Invalid input", f"Could not read region/rotation:\n{exc}")
            return
        alloy_x = self.alloy_x.value() if mat.alloy is not None else None
        n = sum(1 for f in self.fills if f.material_key == key) + 1
        fill_name = (f"{mat.name} — {region.name}" if isinstance(region, ImportedRegion)
                     else f"{mat.name} #{n}")
        fill = Fill(name=fill_name, region=region, material_key=key,
                    rotation=rotation, alloy_x=alloy_x, seed=self.seed_base.value(),
                    carve=self.carve_check.isChecked())
        self.fills.append(fill)
        self._sync_list()
        self._refresh()

    # ── fills list management ─────────────────────────────────────────────────────
    def _sync_list(self):
        self.fills_list.blockSignals(True)
        self.fills_list.clear()
        for f in self.fills:
            item = QListWidgetItem(f"{f.name} — {f.region.describe()}")
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if f.enabled else Qt.Unchecked)
            self.fills_list.addItem(item)
        self.fills_list.blockSignals(False)

    def _on_item_changed(self, item):
        idx = self.fills_list.row(item)
        if 0 <= idx < len(self.fills):
            self.fills[idx].enabled = (item.checkState() == Qt.Checked)
            self._refresh()

    def _selected_index(self):
        idx = self.fills_list.currentRow()
        return idx if 0 <= idx < len(self.fills) else None

    def _move_up(self):
        i = self._selected_index()
        if i and i > 0:
            self.fills[i - 1], self.fills[i] = self.fills[i], self.fills[i - 1]
            self._sync_list(); self.fills_list.setCurrentRow(i - 1); self._refresh()

    def _move_down(self):
        i = self._selected_index()
        if i is not None and i < len(self.fills) - 1:
            self.fills[i + 1], self.fills[i] = self.fills[i], self.fills[i + 1]
            self._sync_list(); self.fills_list.setCurrentRow(i + 1); self._refresh()

    def _remove_fill(self):
        i = self._selected_index()
        if i is not None:
            del self.fills[i]; self._sync_list(); self._refresh()

    def _clear_fills(self):
        if self.fills and QMessageBox.question(self, "Clear", "Remove all regions?") == QMessageBox.Yes:
            self.fills = []; self._sync_list(); self._refresh()

    def _open_crossgrating(self):
        dlg = CrossgratingDialog(self)
        if dlg.exec_() == QDialog.Accepted:
            fills, stats = build_crossgrating(dlg.params())
            self.fills = fills
            self._sync_list()
            self._refresh()
            self._log(f"Crossgrating: {stats}")

    # ── options → WriterOptions ──────────────────────────────────────────────────
    def _gather_opts(self):
        opts = WriterOptions()
        opts.boundary = self.boundary.currentText()
        opts.margin = self.margin.value()
        opts.box_mode = "auto" if self.auto_box.isChecked() else "explicit"
        opts.delete_overlap = self.del_overlap.isChecked()
        opts.overlap_cutoff = self.overlap_cut.value()
        opts.unify_similar_lattices = self.unify_similar_lattices.isChecked()
        opts.lattice_similarity_threshold = self.lattice_similarity_threshold.value() / 100.0
        opts.enforce_min_sep = self.min_sep.isChecked()
        opts.dump_file = self.dump_name.text().strip() or "sample.xyz"
        opts.seed_base = self.seed_base.value()
        if opts.box_mode == "explicit":
            opts.box = geo.union_aabb(self.fills)  # seed explicit box from current extents
        return opts

    # ── refresh: script + info + preview ──────────────────────────────────────────
    def _refresh(self):
        self.opts = self._gather_opts()
        self._update_md_potentials()
        try:
            result = lw.build_script(self.fills, self.opts, title=APP_NAME + " build",
                                     md=self._md_options())
            self.script_view.setPlainText(result.script)
            self._render_info(result)
        except Exception as exc:  # noqa: BLE001
            import traceback
            self.script_view.setPlainText(f"# Error building script:\n# {exc}")
            self.info_view.setPlainText(traceback.format_exc())
        try:
            box = (self.opts.box if self.opts.box_mode == "explicit" and self.opts.box
                   else self._auto_box())
            self.preview.render(
                self.fills,
                box,
                cad_components=self.cad_components,
                selected_cad_name=self.cad_component_combo.currentText(),
            )
        except Exception as exc:  # noqa: BLE001
            pass

    def _auto_box(self):
        b = geo.union_aabb(self.fills)
        m = self.margin.value()
        return (b[0] - m, b[1] + m, b[2] - m, b[3] + m, b[4] - m, b[5] + m)

    def _render_info(self, result):
        lines = [f"Estimated atoms: {result.est_atoms:,}", ""]
        lines.append("Atom types (global):")
        for el, tid in sorted(result.type_table.items(), key=lambda kv: kv[1]):
            lines.append(f"  type {tid} = {el}")
        lines.append("")
        lines.append("Per-region estimate:")
        for name, count in result.per_fill_atoms:
            lines.append(f"  {name}: ~{count:,}")
        if result.warnings:
            lines.append("")
            lines.append("Warnings:")
            for w in result.warnings:
                lines.append(f"  ⚠ {w}")
        # keep any prior run log appended below a divider
        prior = self.info_view.toPlainText()
        log_marker = "\n— run log —\n"
        if log_marker in prior:
            lines.append(log_marker.strip())
            lines.append(prior.split(log_marker, 1)[1])
        self.info_view.setPlainText("\n".join(lines))

    def _log(self, text):
        self.info_view.appendPlainText("\n— run log —\n" + text)

    # ── generate / run ────────────────────────────────────────────────────────────
    def _generate(self):
        if not [f for f in self.fills if f.enabled]:
            QMessageBox.information(self, "Nothing to build", "Add at least one region first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export LAMMPS data", "sample.lmp",
            "LAMMPS data (*.lmp *.data);;All files (*)",
        )
        if not path:
            return
        self._start_export_progress("Preparing LAMMPS data export...")
        try:
            result = export_lammps(path, self.fills, self._gather_opts(),
                                   progress=self._export_progress)
        except Exception:
            self._finish_export_progress("Export failed")
            raise
        self._last_data_path = path
        self._last_type_table = dict(result.type_table)
        self._last_saved_path = path
        md = self._md_options()
        if md is not None and md.backend == "mace":
            from ase.io import write as ase_write
            data_dir = os.path.dirname(os.path.abspath(path))
            self._last_mace_input = os.path.splitext(path)[0] + "_initial.xyz"
            self._last_mace_output = os.path.join(
                data_dir, self.dump_name.text().strip() or "relaxed.xyz")
            ase_write(self._last_mace_input, result.atoms, format="xyz")
            self._last_saved_path = self._last_mace_input
            self._last_dump_paths = []
            self.run_btn.setEnabled(mace_available())
            self._log(f"Wrote MACE input structure {self._last_mace_input}")
        elif md is not None:
            potential = POTENTIALS.get(md.potential_key)
            if potential is None:
                raise ValueError(f"Unknown LAMMPS potential: {md.potential_key}")
            data_dir = os.path.dirname(os.path.abspath(path))
            dump_path = os.path.join(data_dir, self.dump_name.text().strip() or "relaxed.xyz")
            input_path = os.path.splitext(path)[0] + "_relax.in"
            write_minimization_input(
                input_path, path, dump_path, potential, result.type_table,
                boundary=self._gather_opts().boundary,
                overlap_cutoff=self._gather_opts().overlap_cutoff,
                workflow=md.workflow,
                temperature=md.temperature if md.temperature > 0 else 800.0,
                equilibration_steps=md.equilibration_steps,
                timestep_fs=md.timestep_fs if md.timestep_fs > 0 else 1.0,
                seed=md.seed,
            )
            self._last_saved_path = input_path
            self._last_dump_paths = []
            self.run_btn.setEnabled(bool(self.lammps_exe))
            self._log(f"Wrote LAMMPS minimization input {input_path}")
        else:
            self._last_dump_paths = [path]
            self.run_btn.setEnabled(False)
        self.transfer_btn.setEnabled(bool(self.transfer_callback))
        self._log(f"Wrote ASE LAMMPS data file {path}  ({len(result.atoms):,} atoms)")
        QMessageBox.information(self, "Export complete",
                                f"LAMMPS data written to:\n{path}\n\n"
                                f"Atoms: {len(result.atoms):,}")
        self._finish_export_progress()

    def _export_xyz(self):
        if not [f for f in self.fills if f.enabled]:
            QMessageBox.information(self, "Nothing to build", "Add at least one region first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export XYZ data", "sample.xyz",
            "XYZ files (*.xyz);;All files (*)",
        )
        if not path:
            return
        self._start_export_progress("Preparing XYZ export...")
        try:
            result = export_xyz(path, self.fills, self._gather_opts(),
                                progress=self._export_progress)
        except Exception:
            self._finish_export_progress("Export failed")
            raise
        self._last_dump_paths = [path]
        self._last_type_table = dict(result.type_table)
        self.transfer_btn.setEnabled(bool(self.transfer_callback))
        self._log(f"Wrote ASE XYZ file {path}  ({len(result.atoms):,} atoms)")
        QMessageBox.information(self, "Export complete",
                                f"XYZ data written to:\n{path}\n\n"
                                f"Atoms: {len(result.atoms):,}")
        self._finish_export_progress()

    def _export_progress(self, value, message):
        self.export_progress.setValue(value)
        self.export_progress.setFormat(message)
        QApplication.processEvents()

    def _start_export_progress(self, message):
        self.export_progress.setValue(0)
        self.export_progress.setFormat(message)
        QApplication.processEvents()

    def _finish_export_progress(self, message="Ready"):
        self.export_progress.setValue(100)
        self.export_progress.setFormat(message)

    def _run_lammps(self):
        if not self.lammps_exe:
            return
        if not self._last_saved_path or not os.path.exists(self._last_saved_path):
            QMessageBox.information(self, "Generate first", "Generate the .lmp file before running.")
            return
        self.run_btn.setEnabled(False)
        self.stop_relax_btn.setEnabled(True)
        cmd = f"{self.lammps_exe} -in {os.path.basename(self._last_saved_path)}"
        self._log(f"Running: {cmd}")
        # show the live console
        self.console_view.clear()
        self.console_view.appendPlainText(f"$ {cmd}\n")
        self.tabs.setCurrentWidget(self.console_view)
        self._runner = RunnerThread(self.lammps_exe, self._last_saved_path)
        self._runner.output_line.connect(self._on_run_output)
        self._runner.finished_run.connect(self._on_run_finished)
        self._runner.start()

    def _run_relaxation(self):
        if self.md_backend.currentData() == "mace":
            if not self._last_mace_input or not self._last_mace_output:
                QMessageBox.information(self, "Generate first",
                                        "Export the sample with MACE selected first.")
                return
            try:
                with open(self._last_mace_input, encoding="utf-8") as stream:
                    atom_count = int(stream.readline().strip())
            except (OSError, ValueError) as exc:
                QMessageBox.warning(self, "MACE input error", str(exc))
                return
            if atom_count > 10000:
                choice = QMessageBox.question(
                    self, "Large MACE relaxation",
                    f"This structure contains {atom_count:,} atoms. "
                    "MACE relaxation on CPU may take a very long time and "
                    "requires substantial memory. Continue?",
                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
                if choice != QMessageBox.Yes:
                    return
            self.run_btn.setEnabled(False)
            self.stop_relax_btn.setEnabled(True)
            self.console_view.clear()
            self.console_view.appendPlainText("Starting MACE relaxation...\n")
            self.tabs.setCurrentWidget(self.console_view)
            self._runner = MaceRunnerThread(
                self._last_mace_input, self._last_mace_output,
                self.md_mace_model.currentText(), self.md_mace_fmax.value())
            self._runner.output_line.connect(self._on_mace_output)
            self._runner.finished_relaxation.connect(self._on_mace_finished)
            self._runner.start()
            return
        self._run_lammps()

    def _stop_relaxation(self):
        if isinstance(self._runner, MaceRunnerThread):
            self._runner.stop()
        elif isinstance(self._runner, RunnerThread):
            self._runner.stop()

    def _on_mace_finished(self, rc, message):
        self.run_btn.setEnabled(True)
        self.stop_relax_btn.setEnabled(False)
        if rc == 0:
            self._last_dump_paths = [message]
            self._log(f"MACE relaxation complete: {message}")
            QMessageBox.information(self, "MACE relaxation complete",
                                    f"Relaxed XYZ written to:\n{message}")
        else:
            self.console_view.appendPlainText(f"MACE failed: {message}")
            QMessageBox.warning(self, "MACE relaxation failed", message)

    def _on_mace_output(self, line):
        self.console_view.appendPlainText(line)
        sb = self.console_view.verticalScrollBar()
        sb.setValue(sb.maximum())

    def _on_run_output(self, line):
        self.console_view.appendPlainText(line)
        sb = self.console_view.verticalScrollBar()
        sb.setValue(sb.maximum())             # follow the tail

    def _on_run_finished(self, rc):
        self.console_view.appendPlainText(f"\n[process exited with code {rc}]")
        self.run_btn.setEnabled(True)
        self.stop_relax_btn.setEnabled(False)
        d = os.path.dirname(self._last_saved_path)
        dump_name = os.path.basename(self.dump_name.text().strip())
        single = os.path.join(d, dump_name)
        stem, ext = os.path.splitext(dump_name)
        # MD 'thermalize' dumps one file per frozen-phonon frame: <stem>_<step><ext>
        frames = sorted(glob.glob(os.path.join(d, f"{stem}_*{ext or '.xyz'}")))
        self._last_dump_paths = []
        if rc == 0 and frames:
            self._last_dump_paths = frames
            self._log(f"Dumped {len(frames)} frozen-phonon frames: "
                      f"{os.path.basename(frames[0])} … {os.path.basename(frames[-1])}")
            QMessageBox.information(self, "Done",
                                    f"LAMMPS finished.\n{len(frames)} snapshot frame(s) in:\n{d}")
        elif rc == 0 and os.path.exists(single):
            self._last_dump_paths = [single]
            self._log(f"Dump written: {single}")
            QMessageBox.information(self, "Done", f"LAMMPS finished.\nDump: {single}")
        else:
            QMessageBox.warning(self, "LAMMPS", f"Exited with code {rc}. See log.")
        self.transfer_btn.setEnabled(
            bool(self.transfer_callback and self._last_dump_paths))

    # ── embedding hook (SEB-MS host) ───────────────────────────────────────────
    def set_transfer_callback(self, cb):
        """Enable embedded mode: cb(dump_paths, elements) is called by the
        'Transfer to SEB-MS' button after a successful LAMMPS run."""
        self.transfer_callback = cb
        self.transfer_btn.setVisible(cb is not None)
        self.transfer_btn.setEnabled(bool(cb and self._last_dump_paths))
        if cb is not None and hasattr(self, "_quit_action"):
            self._quit_action.setVisible(False)   # host owns app lifetime

    def _transfer_to_host(self):
        if not (self.transfer_callback and self._last_dump_paths):
            return
        # Elements in global (mass-sorted) LAMMPS type order — exactly the
        # order the multislice 'Elements' field expects (type 1 first).
        elements = [el for el, _tid in
                    sorted(self._last_type_table.items(), key=lambda kv: kv[1])]
        self.transfer_callback(list(self._last_dump_paths), elements)

    def refresh_theme(self):
        """Re-apply the stylesheet (embedded host calls this on theme switch)."""
        self.setStyleSheet(theme.build_stylesheet())

    # ── project save / load ────────────────────────────────────────────────────────
    def _fill_to_dict(self, f):
        return {"name": f.name, "shape": type(f.region).__name__,
                "region": asdict(f.region), "material_key": f.material_key,
                "rotation": f.rotation, "alloy_x": f.alloy_x, "seed": f.seed,
                "enabled": f.enabled, "carve": f.carve}

    def _fill_from_dict(self, d):
        region = _REGION_CLASSES[d["shape"]](**d["region"])
        return Fill(name=d["name"], region=region, material_key=d["material_key"],
                    rotation=d.get("rotation", {"mode": "none"}), alloy_x=d.get("alloy_x"),
                    seed=d.get("seed", 12345), enabled=d.get("enabled", True),
                    carve=d.get("carve", False))

    def _save_project(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save project", "sample.saproj",
                                              "Sample Assembler project (*.saproj);;JSON (*.json)")
        if not path:
            return
        data = {"fills": [self._fill_to_dict(f) for f in self.fills],
                "options": {"boundary": self.boundary.currentText(),
                            "auto_box": self.auto_box.isChecked(),
                            "margin": self.margin.value(),
                            "delete_overlap": self.del_overlap.isChecked(),
                            "overlap_cutoff": self.overlap_cut.value(),
                            "unify_similar_lattices": self.unify_similar_lattices.isChecked(),
                            "lattice_similarity_threshold": self.lattice_similarity_threshold.value(),
                            "dump_file": self.dump_name.text(),
                            "seed_base": self.seed_base.value()},
                "md": {"enable": self.md_enable.isChecked(),
                       "backend": self.md_backend.currentData(),
                       "potential": self.md_potential.currentData(),
                       "mace_model": self.md_mace_model.currentText(),
                       "mace_fmax": self.md_mace_fmax.value(),
                       "workflow": self.md_workflow.currentData(),
                       "temperature": self.md_temp.value(),
                       "peak_T": self.md_peakT.value(),
                       "timestep": self.md_timestep.value(),
                       "ensemble": self.md_ensemble.currentText(),
                       "n_equil": self.md_nequil.value(),
                       "n_snapshots": self.md_nsnap.value(),
                       "interval": self.md_interval.value(),
                       "seed": self.md_seed.value()}}
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2)
        self._log(f"Saved project: {path}")

    def _load_project(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load project", "",
                                              "Sample Assembler project (*.saproj *.json);;All files (*)")
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as fh:
                data = json.load(fh)
            self.fills = [self._fill_from_dict(d) for d in data.get("fills", [])]
            o = data.get("options", {})
            self.boundary.setCurrentText(o.get("boundary", "p p p"))
            self.auto_box.setChecked(o.get("auto_box", True))
            self.margin.setValue(o.get("margin", 0.0))
            self.del_overlap.setChecked(o.get("delete_overlap", True))
            self.overlap_cut.setValue(o.get("overlap_cutoff", 1.2))
            self.unify_similar_lattices.setChecked(o.get("unify_similar_lattices", True))
            self.lattice_similarity_threshold.setValue(
                o.get("lattice_similarity_threshold", 5.0)
            )
            self.dump_name.setText(o.get("dump_file", "sample.xyz"))
            self.seed_base.setValue(o.get("seed_base", 12345))
            m = data.get("md", {})
            backend = self.md_backend.findData(m.get("backend", "lammps"))
            if backend >= 0:
                self.md_backend.setCurrentIndex(backend)
            self.md_enable.setChecked(m.get("enable", False))
            self.md_body.setEnabled(self.md_enable.isChecked())
            model = self.md_mace_model.findText(m.get("mace_model", "medium"))
            if model >= 0:
                self.md_mace_model.setCurrentIndex(model)
            self.md_mace_fmax.setValue(m.get("mace_fmax", 0.05))
            iw = self.md_workflow.findData(m.get("workflow", "thermalize"))
            if iw >= 0:
                self.md_workflow.setCurrentIndex(iw)
            self.md_temp.setValue(m.get("temperature", 300.0))
            self.md_peakT.setValue(m.get("peak_T", 6000.0))
            self.md_timestep.setValue(m.get("timestep", 0.0))
            ie = self.md_ensemble.findText(m.get("ensemble", "nvt"))
            if ie >= 0:
                self.md_ensemble.setCurrentIndex(ie)
            self.md_nequil.setValue(m.get("n_equil", 20000))
            self.md_nsnap.setValue(m.get("n_snapshots", 10))
            self.md_interval.setValue(m.get("interval", 500))
            self.md_seed.setValue(m.get("seed", 90210))
            self._pending_md_potential = m.get("potential")
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Load failed", str(exc))
            return
        self._sync_list(); self._refresh()
        # potential combo is populated by _refresh → now restore the saved selection
        if getattr(self, "_pending_md_potential", None):
            i = self.md_potential.findData(self._pending_md_potential)
            if i >= 0:
                self.md_potential.setCurrentIndex(i)
        self._log(f"Loaded project: {path}")

    # ── tools ────────────────────────────────────────────────────────────────────
    def _split_xyz_tool(self):
        path, _ = QFileDialog.getOpenFileName(self, "Multi-timestep .xyz", "",
                                              "XYZ dump (*.xyz);;All files (*)")
        if not path:
            return
        outdir = QFileDialog.getExistingDirectory(self, "Output directory for snapshots")
        if not outdir:
            return
        prefix, _ = os.path.splitext(os.path.basename(path))
        n = split_xyz(path, outdir, prefix + "_")
        self._log(f"Split {path} -> {n} snapshots in {outdir}")
        QMessageBox.information(self, "Done", f"Wrote {n} snapshot files to:\n{outdir}")

    def _about(self):
        QMessageBox.about(self, "About " + APP_NAME,
                          f"<b>{APP_NAME}</b><br>Design TEM specimens and export "
                          "ASE-generated LAMMPS data files.<br><br>"
                          "Materials, geometry, lattice rotation, overlap deletion, and a "
                          "gold-on-carbon crossgrating preset.")


# ── multi-timestep .xyz splitter (ported from xyz_parser.py) ───────────────────
def split_xyz(input_file, out_dir=".", prefix="snapshot_"):
    with open(input_file, "r") as f:
        lines = f.readlines()
    blocks, current = [], []
    for line in lines:
        if line.startswith("ITEM: TIMESTEP") and current:
            blocks.append(current); current = []
        current.append(line)
    if current:
        blocks.append(current)
    for i, block in enumerate(blocks, start=1):
        out = os.path.join(out_dir, f"{prefix}{i:05d}.xyz")
        with open(out, "w") as fh:
            fh.writelines(block)
    return len(blocks)


def main():
    app = QApplication(sys.argv)
    from PyQt5.QtWidgets import QStyleFactory
    app.setStyle(QStyleFactory.create("Fusion"))
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
