#!/usr/bin/env python3
"""
Self-contained theme for the Sample Assembler app.

Mirrors the "Obsidian" palette used by the SEB multislice GUI (app.py):
dark gray #2b2b2b background + desaturated orange #d88a55 accent, Fusion style.
This module is standalone (the multislice app's theme.py is NOT a dependency)
so the Sample Assembler can run on its own.

Always pull colors from get_theme() — never hardcode hex in the UI code.
"""

# ── Semantic color tokens (Obsidian) ──────────────────────────────────────────
_OBSIDIAN = {
    "bg_app":        "#2b2b2b",   # window background
    "bg_panel":      "#323232",   # side panels
    "bg_card":       "#383838",   # cards / group boxes
    "bg_input":      "#2f2f2f",   # text fields, combos
    "bg_hover":      "#444444",
    "rule_soft":     "#4a4a4a",   # subtle separators / borders
    "rule_hard":     "#1f1f1f",
    "fg_primary":    "#ececec",   # main text
    "fg_secondary":  "#c4c4c4",
    "fg_tertiary":   "#9a9a9a",
    "fg_quaternary": "#6f6f6f",   # hints / disabled
    "accent":        "#9e7b64",   # desaturated orange
    "accent_soft":   "#b9744a",
    "accent_text":   "#1d1d1d",   # text on accent fills
    "ok":            "#7fae6e",   # success green
    "warn":          "#d8a657",   # warning amber
    "danger":        "#c75d52",   # destructive red
    "info":          "#60a0c0",   # informational blue
}

_THEME = dict(_OBSIDIAN)


def get_theme():
    """Return the active token dict (copy-safe to read; do not mutate)."""
    return _THEME


def color(token):
    """Convenience accessor with a sane fallback."""
    return _THEME.get(token, "#ff00ff")


# ── Stylesheet generation ──────────────────────────────────────────────────────
def build_stylesheet():
    t = _THEME
    return f"""
    QWidget {{
        background-color: {t['bg_app']};
        color: {t['fg_primary']};
        font-family: 'Segoe UI', 'Inter', sans-serif;
        font-size: 13px;
    }}
    QMainWindow, QDialog {{ background-color: {t['bg_app']}; }}

    QLabel {{ background: transparent; color: {t['fg_primary']}; }}
    QLabel[role="title"]     {{ font-size: 17px; font-weight: 600; color: {t['fg_primary']}; }}
    QLabel[role="subtitle"]  {{ font-size: 12px; color: {t['fg_tertiary']}; }}
    QLabel[role="hint"]      {{ font-size: 11px; color: {t['fg_quaternary']}; }}
    QLabel[role="section"]   {{ font-size: 12px; font-weight: 600; color: {t['accent']};
                                text-transform: uppercase; letter-spacing: 1px; }}

    QGroupBox {{
        background-color: {t['bg_card']};
        border: 1px solid {t['rule_soft']};
        border-radius: 8px;
        margin-top: 14px;
        padding: 10px 10px 10px 10px;
        font-weight: 600;
    }}
    QGroupBox::title {{
        subcontrol-origin: margin;
        left: 10px; top: 2px;
        padding: 0 5px;
        color: {t['accent']};
    }}

    QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QPlainTextEdit, QTextEdit {{
        background-color: {t['bg_input']};
        border: 1px solid {t['rule_soft']};
        border-radius: 5px;
        padding: 4px 6px;
        color: {t['fg_primary']};
        selection-background-color: {t['accent']};
        selection-color: {t['accent_text']};
    }}
    QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
        border: 1px solid {t['accent']};
    }}
    QComboBox::drop-down {{ border: none; width: 18px; }}
    QComboBox QAbstractItemView {{
        background-color: {t['bg_input']};
        border: 1px solid {t['rule_soft']};
        selection-background-color: {t['accent_soft']};
        selection-color: {t['accent_text']};
        outline: none;
    }}

    QPushButton {{
        background-color: {t['bg_hover']};
        border: 1px solid {t['rule_soft']};
        border-radius: 6px;
        padding: 6px 14px;
        color: {t['fg_primary']};
    }}
    QPushButton:hover {{ background-color: {t['rule_soft']}; }}
    QPushButton:pressed {{ background-color: {t['rule_hard']}; }}
    QPushButton:disabled {{ color: {t['fg_quaternary']}; background-color: {t['bg_card']}; }}

    QPushButton[kind="primary"] {{
        background-color: {t['accent']};
        color: {t['accent_text']};
        border: 1px solid {t['accent']};
        font-weight: 600;
    }}
    QPushButton[kind="primary"]:hover {{ background-color: {t['accent_soft']}; }}
    QPushButton[kind="danger"] {{
        background-color: {t['danger']};
        color: #fff; border: 1px solid {t['danger']}; font-weight: 600;
    }}

    QListWidget, QTreeWidget, QTableWidget {{
        background-color: {t['bg_input']};
        border: 1px solid {t['rule_soft']};
        border-radius: 6px;
        outline: none;
    }}
    QListWidget::item {{ padding: 5px 6px; border-radius: 4px; }}
    QListWidget::item:selected {{ background-color: {t['accent_soft']}; color: {t['accent_text']}; }}
    QListWidget::item:hover {{ background-color: {t['bg_hover']}; }}

    QTabWidget::pane {{ border: 1px solid {t['rule_soft']}; border-radius: 6px; top: -1px; }}
    QTabBar::tab {{
        background: {t['bg_card']}; color: {t['fg_secondary']};
        padding: 6px 14px; border: 1px solid {t['rule_soft']};
        border-bottom: none; border-top-left-radius: 6px; border-top-right-radius: 6px;
    }}
    QTabBar::tab:selected {{ background: {t['bg_app']}; color: {t['accent']}; font-weight: 600; }}

    QCheckBox, QRadioButton {{ spacing: 6px; background: transparent; }}
    QCheckBox::indicator, QRadioButton::indicator {{ width: 15px; height: 15px; }}
    QCheckBox::indicator {{ border: 1px solid {t['rule_soft']}; border-radius: 3px; background: {t['bg_input']}; }}
    QCheckBox::indicator:checked {{ background: {t['accent']}; border: 1px solid {t['accent']}; }}

    QScrollBar:vertical {{ background: {t['bg_app']}; width: 11px; margin: 0; }}
    QScrollBar::handle:vertical {{ background: {t['rule_soft']}; border-radius: 5px; min-height: 24px; }}
    QScrollBar::handle:vertical:hover {{ background: {t['fg_quaternary']}; }}
    QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; }}
    QScrollBar:horizontal {{ background: {t['bg_app']}; height: 11px; margin: 0; }}
    QScrollBar::handle:horizontal {{ background: {t['rule_soft']}; border-radius: 5px; min-width: 24px; }}

    QSplitter::handle {{ background: {t['rule_soft']}; }}
    QHeaderView::section {{
        background-color: {t['bg_card']}; color: {t['fg_secondary']};
        padding: 4px 6px; border: none; border-right: 1px solid {t['rule_soft']};
    }}
    QToolTip {{
        background-color: {t['bg_card']}; color: {t['fg_primary']};
        border: 1px solid {t['accent_soft']}; padding: 4px;
    }}
    """


# ── matplotlib styling to match the dark theme ─────────────────────────────────
def apply_mpl_style(fig=None, ax=None):
    """Tint a matplotlib figure/axes to the app palette."""
    t = _THEME
    if fig is not None:
        fig.patch.set_facecolor(t["bg_card"])
    if ax is not None:
        ax.set_facecolor(t["bg_input"])
        for spine in ax.spines.values():
            spine.set_color(t["rule_soft"])
        ax.tick_params(colors=t["fg_tertiary"], labelsize=8)
        ax.xaxis.label.set_color(t["fg_secondary"])
        ax.yaxis.label.set_color(t["fg_secondary"])
        if hasattr(ax, "zaxis"):
            ax.zaxis.label.set_color(t["fg_secondary"])
        ax.title.set_color(t["fg_primary"])
