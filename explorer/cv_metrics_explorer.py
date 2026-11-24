#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CV Metrics Explorer (Tkinter + Matplotlib) - stage 2 of the voltammetry pipeline.

Reads the per-cycle CSV files exported by CV Analyzer (stage 1) and compares
the parameters between experiments, anions and cycle number.

- Each CSV file is one experiment.
- Rows are cycles (expects a 'cycle' column; otherwise uses row index + 1).
- Experiment metadata is parsed from the file name, which is the name CV
  Analyzer's export dialog writes:
    Anion_Station_Electrode_Measurement_Laborant.csv
  Example: Carbonate_A_20002_ACT-night_Maya-R..csv

Tabs:
1) Params vs Cycles (one experiment, multiple parameters)
2) Experiments vs Cycles (one parameter, multiple experiments) - the main view
3) Scatter (two parameters, many experiments) - experimental
4) Correlation Matrix (filter by experiment type + anion) - experimental
5) Target Correlations (target vs all others) - experimental
6) Appendix Page (arrange graphs from the other tabs on a printable page)

Dependencies: Python 3.9+ with Tkinter, pandas, numpy, matplotlib, pillow
"""

import os
import re
import json
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import matplotlib.colors as mcolors


import numpy as np
import pandas as pd
import itertools
import pickle
import textwrap

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, ScalarFormatter, MultipleLocator
from matplotlib.lines import Line2D
from matplotlib.backends.backend_agg import FigureCanvasAgg
from PIL import Image



# ========== SECTION 1: COLOR MAPPING ==========
def group_color_map(group_names):
    """
    Deterministic mapping from group name -> matplotlib color (uses default color cycle).
    Returns dict: group -> color
    """
    cycle = matplotlib.rcParams.get("axes.prop_cycle", None)
    colors = None
    if cycle is not None:
        try:
            colors = cycle.by_key().get("color", None)
        except Exception:
            colors = None
    if not colors:
        colors = ["C0","C1","C2","C3","C4","C5","C6","C7","C8","C9"]

    groups_sorted = sorted({g for g in group_names if g is not None})
    cmap = {}
    it = itertools.cycle(colors)
    for g in groups_sorted:
        cmap[g] = next(it)
    return cmap


# ========== SECTION 1B: FAVORITE PARAMETERS (persisted) ==========
try:
    _APP_DIR = os.path.dirname(os.path.abspath(__file__))
except Exception:
    _APP_DIR = os.getcwd()

FAV_PATH = os.path.join(_APP_DIR, "cv_explorer_favorites.json")


def load_favorite_params() -> set:
    try:
        with open(FAV_PATH, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, list):
            return set(data)
    except Exception:
        pass
    return set()


def save_favorite_params(favs: set):
    try:
        with open(FAV_PATH, "w", encoding="utf-8") as f:
            json.dump(sorted(favs), f, ensure_ascii=False, indent=2)
    except Exception as e:
        print(f"[WARN] Could not save favorites: {e}")


FAVORITE_PARAMS = load_favorite_params()
_PARAM_TREE_INSTANCES = []  # all live ParamTree widgets, so favorites stay in sync everywhere

CONFIG_PATH = os.path.join(_APP_DIR, "cv_explorer_config.json")


def load_config() -> dict:
    try:
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict):
            return data
    except Exception:
        pass
    return {}


def save_config(cfg: dict):
    try:
        with open(CONFIG_PATH, "w", encoding="utf-8") as f:
            json.dump(cfg, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print(f"[WARN] Could not save config: {e}")


# ========== SECTION 2: PARAM TREE (CUSTOM) ==========
# Potentials (V vs RHE) at which CV Analyzer samples the forward / backward
# currents. Must match TARGET_E in analyzer/cv_analyzer.py; the column names
# are built the same way (1.1 -> "1p1", 1.375 -> "1p375").
TARGET_E = [1.0, 1.1, 1.2, 1.375, 1.4, 1.45, 1.5, 1.55, 1.6]


def _ekey(E):
    return f"{float(E):g}".replace(".", "p").replace("-", "m")


def _fixed_e_cols(kind):
    return [f"raw_I_{kind}_{_ekey(E)}" for E in TARGET_E]


PARAM_TREE = {
    "Raw": {
        "Onset": {
            "Percentiles": [
                "onset_raw_E_p5", "onset_raw_E_p10",
                "onset_raw_E_p15", "onset_raw_E_p25"
            ],
            "Curvature (k)": [
                "onset_raw_E_k2", "onset_raw_E_k3", "onset_raw_E_k4",
                "onset_raw_E_k5", "onset_raw_E_k8", "onset_raw_E_k12",
                "onset_raw_E_k15", "onset_raw_E_k20"
            ],
            "Junction": [
                "raw_junction_E"
            ]
        },
        "Activity": {
            "I @ E (Forward)": _fixed_e_cols("fwd"),
            "I @ E (Backward)": _fixed_e_cols("bwd"),
            "I @ E (Forward-Backward)": _fixed_e_cols("for_back"),
            "Junction currents": [
                "raw_junction_I_fwd","raw_junction_I_bwd","raw_junction_for_back"
            ]
        },
        "Peaks": {
            "Peak E": ["anodic_E_raw", "cathodic_E_raw"],
            "Peak I": ["anodic_I_raw", "cathodic_I_raw"],
            "I ratio": ["I_ratio_raw"],
            "ΔEp": ["dEp_raw"]
        },
        "Charge": {
            "Q end": ["Q_end_raw"]
        }
    },

    # OER-corrected values. CV Analyzer computes the fixed-potential currents
    # on the raw curve only, so there is no "Activity" group here.
    "Corr": {
        "Onset": {
            "Percentiles": [
                "onset_corr_E_p5", "onset_corr_E_p10",
                "onset_corr_E_p15", "onset_corr_E_p25"
            ],
            "Curvature (k)": [
                "onset_corr_E_k2", "onset_corr_E_k3",
                "onset_corr_E_k4", "onset_corr_E_k5",
                "onset_corr_E_k8", "onset_corr_E_k12",
                "onset_corr_E_k15", "onset_corr_E_k20"
            ]
        },
        "Peaks": {
            "Peak E": ["anodic_E_corr", "cathodic_E_corr"],
            "Peak I": ["anodic_I_corr", "cathodic_I_corr"],
            "I ratio": ["I_ratio_corr"],
            "ΔEp": ["dEp_corr"]
        },
        "Charge": {
            "Q end": ["Q_end_corr"]
        }
    }
}

_TREE_PARAMS = {p for groups in PARAM_TREE.values() for blocks in groups.values()
                for params in blocks.values() for p in params}

# Numeric columns of the loaded files that are not in PARAM_TREE (for example
# from an older or newer Analyzer version). Shown under "Other columns" so
# that no column is hidden from the parameter trees.
LOADED_EXTRA_PARAMS = set()


# ========== SECTION 3: TREE WIDGETS ==========
class ParamTree(ttk.Frame):
    """
    A simple tree widget for selecting parameter names.
    Leaf nodes are the parameter strings (contain '_').
    Includes a "Favorites" section (persisted to disk) for quick access,
    plus a right-click menu / button to toggle favorites for the current
    selection.
    """
    def __init__(self, parent):
        super().__init__(parent)

        self.tree = ttk.Treeview(self, show="tree", selectmode="extended")
        ysb = ttk.Scrollbar(self, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=ysb.set)

        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        ysb.pack(side=tk.RIGHT, fill=tk.Y)

        # Right click -> context menu with "toggle favorite"
        self.tree.bind("<Button-3>", self._on_right_click)
        # Double click on a parameter -> quick toggle favorite
        self.tree.bind("<Double-Button-1>", self._on_double_click)

        _PARAM_TREE_INSTANCES.append(self)
        self._build_tree()

    def destroy(self):
        try:
            _PARAM_TREE_INSTANCES.remove(self)
        except ValueError:
            pass
        super().destroy()

    def _build_tree(self):
        self.tree.delete(*self.tree.get_children(""))

        if FAVORITE_PARAMS:
            fav_id = self.tree.insert("", "end", text="★ Favorites", open=True)
            for p in sorted(FAVORITE_PARAMS):
                self.tree.insert(fav_id, "end", text=p, values=(p,))

        for domain, groups in PARAM_TREE.items():
            dom_id = self.tree.insert("", "end", text=domain, open=True)
            for group, blocks in groups.items():
                g_id = self.tree.insert(dom_id, "end", text=group, open=False)
                for subgroup, params in blocks.items():
                    sg_id = self.tree.insert(g_id, "end", text=subgroup, open=False)
                    for p in params:
                        # Star-prefix params that are already favorited, so they're
                        # recognizable wherever they show up in the tree (not just
                        # inside the Favorites folder).
                        label = f"★ {p}" if p in FAVORITE_PARAMS else p
                        self.tree.insert(sg_id, "end", text=label, values=(p,))

        if LOADED_EXTRA_PARAMS:
            other_id = self.tree.insert("", "end", text="Other columns", open=False)
            for p in sorted(LOADED_EXTRA_PARAMS):
                label = f"★ {p}" if p in FAVORITE_PARAMS else p
                self.tree.insert(other_id, "end", text=label, values=(p,))

    def get_selected_params(self) -> list:
        # Leaf items always carry the real param name in "values"; the displayed
        # text may have a "★ " prefix, so read from values, not text.
        selected = []
        for item in self.tree.selection():
            vals = self.tree.item(item, "values")
            if vals:
                selected.append(vals[0])
        # de-duplicate (a param can appear both under Favorites and its normal spot)
        seen = set()
        out = []
        for p in selected:
            if p not in seen:
                out.append(p)
                seen.add(p)
        return out

    def clear_selection(self):
        self.tree.selection_remove(self.tree.selection())

    def select_all_leaves(self):
        leaves = []
        def walk(node=""):
            for child in self.tree.get_children(node):
                if self.tree.item(child, "values"):
                    leaves.append(child)
                walk(child)
        walk("")
        self.tree.selection_set(leaves)

    def toggle_favorites_for_selection(self):
        """Toggle favorite status for every selected leaf param."""
        changed = False
        for item in self.tree.selection():
            vals = self.tree.item(item, "values")
            if not vals:
                continue
            p = vals[0]
            if p in FAVORITE_PARAMS:
                FAVORITE_PARAMS.discard(p)
            else:
                FAVORITE_PARAMS.add(p)
            changed = True
        if changed:
            save_favorite_params(FAVORITE_PARAMS)
            for inst in list(_PARAM_TREE_INSTANCES):
                try:
                    inst._build_tree()
                except Exception:
                    pass

    def _on_double_click(self, event):
        item = self.tree.identify_row(event.y)
        if not item:
            return
        vals = self.tree.item(item, "values")
        if not vals:
            return  # double-clicked a category header, not a param
        self.tree.selection_set(item)
        self.toggle_favorites_for_selection()
        return "break"

    def _on_right_click(self, event):
        item = self.tree.identify_row(event.y)
        if not item:
            return
        vals = self.tree.item(item, "values")
        if not vals:
            return
        if item not in self.tree.selection():
            self.tree.selection_set(item)

        menu = tk.Menu(self, tearoff=0)
        is_fav = vals[0] in FAVORITE_PARAMS
        label = "Remove from favorites ☆" if is_fav else "Add to favorites ★"
        menu.add_command(label=label, command=self.toggle_favorites_for_selection)
        try:
            menu.tk_popup(event.x_root, event.y_root)
        finally:
            menu.grab_release()


class AnionExperimentTree(ttk.Frame):
    """
    Tree for selecting experiments grouped by anion.
    Leaf nodes are experiment stems (values contain the stem).
    """
    def __init__(self, parent):
        super().__init__(parent)

        self.tree = ttk.Treeview(self, show="tree", selectmode="extended")
        ysb = ttk.Scrollbar(self, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=ysb.set)

        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        ysb.pack(side=tk.RIGHT, fill=tk.Y)

        self._anion_nodes = {}
        self._anion_labels = {}  # node_id -> plain anion name (without checkbox glyph/count)

        # Clicking an anion header toggles selection of every experiment under it,
        # like a group checkbox. Clicking a leaf (an experiment) behaves normally.
        self.tree.bind("<Button-1>", self._on_click)
        self.tree.bind("<<TreeviewSelect>>", self._update_group_checkboxes)

    def set_experiments(self, stems: list, experiments_by_stem: dict):
        prev = set(self.get_selected_stems())

        self.tree.delete(*self.tree.get_children(""))
        self._anion_nodes.clear()
        self._anion_labels.clear()

        groups = {}
        for stem in stems:
            exp = experiments_by_stem.get(stem)
            anion = exp.anion if exp is not None else "Unknown"
            groups.setdefault(anion, []).append(stem)

        for anion in sorted(groups.keys()):
            a_id = self.tree.insert("", "end", text=anion, open=False)
            self._anion_nodes[anion] = a_id
            self._anion_labels[a_id] = anion
            for stem in sorted(groups[anion]):
                self.tree.insert(a_id, "end", text=stem, values=(stem,))

        to_select = []
        def walk(node=""):
            for child in self.tree.get_children(node):
                vals = self.tree.item(child, "values")
                if vals and vals[0] in prev:
                    to_select.append(child)
                walk(child)
        walk("")
        if to_select:
            self.tree.selection_set(to_select)
        self._update_group_checkboxes()

    def _on_click(self, event):
        # Don't hijack clicks on the expand/collapse arrow (+/-) — let the user
        # open the group and then click individual experiments normally.
        elem = self.tree.identify("element", event.x, event.y)
        if "indicator" in elem:
            return

        item = self.tree.identify_row(event.y)
        if not item:
            return
        if item in self._anion_labels:
            # Anion header clicked (on its label/text, not the arrow): toggle-select
            # all experiments in that group, like a checkbox ("all selected" ->
            # deselect all, otherwise -> select all).
            children = self.tree.get_children(item)
            if not children:
                return "break"
            current = set(self.tree.selection())
            all_selected = all(c in current for c in children)
            if all_selected:
                new_sel = current - set(children)
            else:
                new_sel = current | set(children)
            self.tree.selection_set(list(new_sel))
            return "break"
        # Otherwise let the default Treeview click/selection behavior run.

    def _update_group_checkboxes(self, event=None):
        sel = set(self.tree.selection())
        for a_id, anion in self._anion_labels.items():
            children = self.tree.get_children(a_id)
            n = len(children)
            if n == 0:
                continue
            selected_n = sum(1 for c in children if c in sel)
            if selected_n == 0:
                box = "☐"
            elif selected_n == n:
                box = "☑"
            else:
                box = "◐"
            self.tree.item(a_id, text=f"{box} {anion} ({n})")

    def get_selected_stems(self) -> list:
        stems = []
        for item in self.tree.selection():
            txt = self.tree.item(item, "text")
            vals = self.tree.item(item, "values")
            if vals:
                stems.append(vals[0])
            elif "_" in txt:
                stems.append(txt)

        seen=set()
        out=[]
        for s in stems:
            if s not in seen:
                out.append(s); seen.add(s)
        return out

    def clear_selection(self):
        self.tree.selection_remove(self.tree.selection())

    def select_all_leaves(self):
        leaves = []
        def walk(node=""):
            for child in self.tree.get_children(node):
                if self.tree.item(child, "values"):
                    leaves.append(child)
                walk(child)
        walk("")
        self.tree.selection_set(leaves)


# ========== SECTION 4: METADATA PARSING ==========
META_RE = re.compile(
    r"^(?P<anion>[^_]+)_(?P<station>[^_]+)_(?P<electrode>[^_]+)_(?P<experimenttype>[^_]+)_(?P<operator>[^_]+)$"
)

# Experiment type (as written in the file name / shown in the filter) -> wording used in graph titles.
# The measurement types set in CV Analyzer (File -> Naming settings...) are read from the
# same settings file, so a type only has to be described once. The rows below are fallbacks.
# Matching is case-insensitive and by substring; add more rows here if needed.
ANALYZER_SETTINGS_PATH = os.path.join(os.path.expanduser("~"), ".cv_overlay_settings.json")

EXPERIMENT_TYPE_DESCRIPTIONS = [
    ("ACT-night", "activation, no additive"),
    ("STB-day", "stability, with additive"),
    ("UAL-ACT-night", "activation without urea"),
    ("UOS-night", "activation with urea"),
]


def _analyzer_type_descriptions() -> list:
    """[(measurement name, wording), ...] from CV Analyzer's naming settings, if any."""
    try:
        with open(ANALYZER_SETTINGS_PATH, "r", encoding="utf-8") as f:
            meas = (json.load(f).get("naming") or {}).get("meas") or []
    except Exception:
        return []
    out = []
    for row in meas:
        if isinstance(row, (list, tuple)) and len(row) >= 3 and str(row[1]).strip() and str(row[2]).strip():
            out.append((str(row[1]).strip(), str(row[2]).strip()))
    return out


def describe_experiment_type(exp_type: str) -> str:
    t = (exp_type or "").lower()
    rows = _analyzer_type_descriptions() + EXPERIMENT_TYPE_DESCRIPTIONS
    # longest key first, so "UAL-ACT-night" wins over "ACT-night"
    for key, desc in sorted(rows, key=lambda r: -len(r[0])):
        if key.lower() in t:
            return desc
    return ""


def parse_meta_from_stem(stem: str) -> dict:
    m = META_RE.match(stem)
    if not m:
        return {
            "anion": "Unknown",
            "station": "Unknown",
            "electrode": "Unknown",
            "experimenttype": "Unknown",
            "operator": "Unknown",
            "stem": stem,
        }
    d = m.groupdict()
    d["stem"] = stem
    return d


# ========== SECTION 5: UNITS ==========
def infer_unit(param: str):
    """
    Heuristic units inference:
    - E or dEp -> V
    - I -> mA (CV Analyzer keeps EC-Lab's units)
    - Q -> mC (mA * s)
    - ratios -> None (unitless)
    """
    p = (param or "").lower()

    # unitless
    if "ratio" in p:
        return None
    if p.startswith("corr") or "_corr" in p:
        # these are still physical features usually, but correlation results are unitless;
        # here the 'corr' in name means corrected, not correlation. Don't use this.
        pass

    if "dep" in p or p.endswith("_ep") or "ep_" in p:
        return "V"

    # Voltage-like
    if "_e_" in p or p.endswith("_e") or p.startswith("onset_") and "_e_" in p:
        return "V"
    if p.endswith("_e_raw") or p.endswith("_e_corr") or "junction_e" in p:
        return "V"

    # Current-like (CV Analyzer keeps EC-Lab's mA)
    if "_i_" in p or p.endswith("_i") or "junction_i" in p or "junction_for_back" in p:
        return "mA"

    # Charge-like (integral of mA over s = mC)
    if p.startswith("q") or "_q" in p:
        return "mC"

    return None

def label_with_unit(param: str) -> str:
    u = infer_unit(param)
    return f"{param} ({u})" if u else param

def y_label_for_param(param: str) -> str:
    u = infer_unit(param)
    return f"{param} [{u}]" if u else param


# ========== SECTION 6: DATA STORE ==========
class Experiment:
    def __init__(self, path: str, df: pd.DataFrame):
        self.path = path
        self.stem = os.path.splitext(os.path.basename(path))[0]
        self.meta = parse_meta_from_stem(self.stem)

        self.df = df.copy()
        if "cycle" not in self.df.columns:
            self.df.insert(0, "cycle", np.arange(1, len(self.df) + 1))

        self.df["cycle"] = pd.to_numeric(self.df["cycle"], errors="coerce")
        self.df = self.df.dropna(subset=["cycle"])
        self.df["cycle"] = self.df["cycle"].astype(int)

        self.numeric_cols = [
            c for c in self.df.columns
            if c not in ("cycle",) and pd.api.types.is_numeric_dtype(self.df[c])
        ]

    @property
    def anion(self) -> str:
        return self.meta.get("anion", "Unknown")

    @property
    def experimenttype(self) -> str:
        return self.meta.get("experimenttype", "Unknown")

    def label(self) -> str:
        return f"{self.stem}"


# ========== SECTION 7: UI HELPERS ==========
def safe_read_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if df.shape[1] == 1:
        df = pd.read_csv(path, sep=";")
    return df

def set_combobox_values(cb: ttk.Combobox, values: list, keep=None):
    cb["values"] = values
    if keep is not None and keep in values:
        cb.set(keep)
    elif values:
        cb.set(values[0])
    else:
        cb.set("")


# ========== SECTION 8: MAIN APP ==========
class CVExplorer(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("CV Metrics Explorer")
        self.geometry("1280x820")
        self.minsize(1100, 700)

        self.experiments: dict = {}  # stem -> Experiment

        # ---------- Correlation filter state (Matrix tab) ----------
        self.corr_abs_var = tk.BooleanVar(value=True)
        self.corr_low_var = tk.DoubleVar(value=0.70)
        self.corr_high_var = tk.DoubleVar(value=1.00)
        self._corr_last_raw = None
        self._corr_last_title = ""
        self._corr_cbar = None

        # ---------- Target correlations (Special tab) ----------
        self.tcorr_abs_var = tk.BooleanVar(value=True)
        self.tcorr_minabs_var = tk.DoubleVar(value=0.00)
        self.tcorr_topn_var = tk.IntVar(value=30)
        self._tcorr_last_series = None
        self._tcorr_last_title = ""

        # ---------- Plot font sizes (global; live) ----------
        self.fs_title  = tk.IntVar(value=12)
        self.fs_label  = tk.IntVar(value=11)
        self.fs_ticks  = tk.IntVar(value=10)
        self.fs_legend = tk.IntVar(value=9)

        # ---------- Legend visibility (per tab; live) ----------
        self.pvc_legend_var = tk.BooleanVar(value=True)
        self.evc_legend_var = tk.BooleanVar(value=True)
        self.sc_legend_var  = tk.BooleanVar(value=True)

        # ---------- Experiments-vs-Cycles: merge legend by anion ----------
        self.evc_legend_by_anion_var = tk.BooleanVar(value=True)

        # ---------- Experiments-vs-Cycles: fixed Y range ----------
        # True => Y axis spans the min/max of the chosen parameter over all loaded experiments
        # of the SAME experiment type as the plotted ones, so graphs of different anions
        # (same measurement type) are directly comparable.
        self.evc_fixed_y_var = tk.BooleanVar(value=True)
        self._evc_last_auto_title = ""  # lets us tell "auto-filled" from "typed by the user"
        self._evc_replot_id = None      # pending (debounced) automatic re-plot of the Experiments-vs-Cycles tab

        # ---------- X tick step (cycles) for the cycle-axis tabs ----------
        self.pvc_xtick_step_var = tk.StringVar(value="")
        self.evc_xtick_step_var = tk.StringVar(value="")

        # cached legend handles/labels per tab, so the "show legend" toggle
        # can restore a legend without needing to re-plot
        self._pvc_legend_hl = None
        self._evc_legend_hl = None
        self._sc_legend_hl = None

        # ---------- Per-tab label controls (editable; live) ----------
        # empty string => "auto"
        self.pvc_title_var  = tk.StringVar(value="")
        self.pvc_xlabel_var = tk.StringVar(value="Cycle")
        self.pvc_ylabel_var = tk.StringVar(value="Value")
        self.pvc_minor_var  = tk.BooleanVar(value=False)

        self.evc_title_var  = tk.StringVar(value="")
        self.evc_xlabel_var = tk.StringVar(value="Cycle")
        self.evc_ylabel_var = tk.StringVar(value="")
        self.evc_minor_var  = tk.BooleanVar(value=False)

        self.sc_title_var   = tk.StringVar(value="")
        self.sc_xlabel_var  = tk.StringVar(value="")
        self.sc_ylabel_var  = tk.StringVar(value="")
        self.sc_minor_var   = tk.BooleanVar(value=False)

        self.corr_title_var  = tk.StringVar(value="")
        self.corr_xlabel_var = tk.StringVar(value="")
        self.corr_ylabel_var = tk.StringVar(value="")
        self.corr_minor_var  = tk.BooleanVar(value=False)

        self.tcorr_title_var  = tk.StringVar(value="")
        self.tcorr_xlabel_var = tk.StringVar(value="Correlation r")
        self.tcorr_ylabel_var = tk.StringVar(value="")
        self.tcorr_minor_var  = tk.BooleanVar(value=False)

        # ---------- auto titles (fallback) ----------
        self._auto = {
            "pvc": {"title": "", "xlabel": "Cycle", "ylabel": "Value"},
            "evc": {"title": "", "xlabel": "Cycle", "ylabel": ""},
            "sc":  {"title": "", "xlabel": "", "ylabel": ""},
            "corr":{"title": "", "xlabel": "", "ylabel": ""},
            "tcorr":{"title": "", "xlabel": "Correlation r", "ylabel": ""},
        }

        self._build_topbar()
        self._build_notebook()

        self._install_live_style_traces()
        self._refresh_filters_and_lists()
        self._auto_load_last_folder()

    # ========== SECTION 8.1: LIVE STYLE ==========
    def _install_live_style_traces(self):
        vars_to_trace = [
            self.fs_title, self.fs_label, self.fs_ticks, self.fs_legend,
            self.pvc_title_var, self.pvc_xlabel_var, self.pvc_ylabel_var, self.pvc_minor_var,
            self.pvc_legend_var, self.pvc_xtick_step_var,
            self.evc_title_var, self.evc_xlabel_var, self.evc_ylabel_var, self.evc_minor_var,
            self.evc_legend_var, self.evc_xtick_step_var,
            self.sc_title_var, self.sc_xlabel_var, self.sc_ylabel_var, self.sc_minor_var,
            self.sc_legend_var,
            self.corr_title_var, self.corr_xlabel_var, self.corr_ylabel_var, self.corr_minor_var,
            self.tcorr_title_var, self.tcorr_xlabel_var, self.tcorr_ylabel_var, self.tcorr_minor_var,
        ]
        for v in vars_to_trace:
            try:
                v.trace_add("write", lambda *_: self._apply_styles_all())
            except Exception:
                pass

        # This one changes *what* is grouped, not just styling, so it triggers a re-plot
        try:
            self.evc_legend_by_anion_var.trace_add("write", lambda *_: self._evc_on_group_toggle())
            self.evc_fixed_y_var.trace_add("write", lambda *_: self._evc_on_group_toggle())
        except Exception:
            pass

    def _evc_on_group_toggle(self):
        try:
            stems = self.evc_exp_tree.get_selected_stems() if hasattr(self, "evc_exp_tree") else []
        except Exception:
            stems = []
        if stems and self.evc_param_cb.get():
            self.plot_experiments_vs_cycles()

    def _apply_styles_all(self):
        # Apply styles to all axes that exist (if plotted)
        self._apply_style_to_one("pvc", getattr(self, "pvc_ax", None), getattr(self, "pvc_canvas", None), x_is_cycle=True)
        self._apply_style_to_one("evc", getattr(self, "evc_ax", None), getattr(self, "evc_canvas", None), x_is_cycle=True)
        self._apply_style_to_one("sc",  getattr(self, "sc_ax", None),  getattr(self, "sc_canvas", None),  x_is_cycle=False)
        self._apply_style_to_one("corr",getattr(self, "corr_ax", None),getattr(self, "corr_canvas", None),x_is_cycle=False, is_matrix=True)
        self._apply_style_to_one("tcorr",getattr(self, "tcorr_ax", None),getattr(self, "tcorr_canvas", None),x_is_cycle=False)

    def _get_user_or_auto(self, tabkey: str, which: str) -> str:
        # which: title/xlabel/ylabel
        user_var = {
            ("pvc","title"): self.pvc_title_var, ("pvc","xlabel"): self.pvc_xlabel_var, ("pvc","ylabel"): self.pvc_ylabel_var,
            ("evc","title"): self.evc_title_var, ("evc","xlabel"): self.evc_xlabel_var, ("evc","ylabel"): self.evc_ylabel_var,
            ("sc","title"):  self.sc_title_var,  ("sc","xlabel"):  self.sc_xlabel_var,  ("sc","ylabel"):  self.sc_ylabel_var,
            ("corr","title"):self.corr_title_var,("corr","xlabel"):self.corr_xlabel_var,("corr","ylabel"):self.corr_ylabel_var,
            ("tcorr","title"):self.tcorr_title_var,("tcorr","xlabel"):self.tcorr_xlabel_var,("tcorr","ylabel"):self.tcorr_ylabel_var,
        }.get((tabkey, which), None)

        if user_var is None:
            return self._auto.get(tabkey, {}).get(which, "")

        val = user_var.get().strip()
        if val == "":
            return self._auto.get(tabkey, {}).get(which, "")
        return val

    def _minor_enabled(self, tabkey: str) -> bool:
        v = {
            "pvc": self.pvc_minor_var,
            "evc": self.evc_minor_var,
            "sc": self.sc_minor_var,
            "corr": self.corr_minor_var,
            "tcorr": self.tcorr_minor_var,
        }.get(tabkey, None)
        return bool(v.get()) if v is not None else False

    def _disable_offset(self, ax):
        # Prevent "+1.468" offset text. Only on numeric axes: an axis with
        # text tick labels (parameter names in the correlation tabs) keeps
        # its own formatter, otherwise the names turn into 0, 1, 2, ...
        for axis in (ax.xaxis, ax.yaxis):
            if not isinstance(axis.get_major_formatter(), ScalarFormatter):
                continue
            fmt = ScalarFormatter(useOffset=False)
            fmt.set_scientific(False)
            axis.set_major_formatter(fmt)
            try:
                axis.offsetText.set_visible(False)
            except Exception:
                pass

    def _apply_style_to_one(self, tabkey: str, ax, canvas, x_is_cycle: bool, is_matrix: bool=False):
        if ax is None or canvas is None:
            return

        # If nothing plotted yet, don't fight default formatting too hard
        # (still safe to apply)
        title  = self._get_user_or_auto(tabkey, "title")
        xlabel = self._get_user_or_auto(tabkey, "xlabel")
        ylabel = self._get_user_or_auto(tabkey, "ylabel")

        # Apply texts
        if title:
            ax.set_title(title, fontsize=int(self.fs_title.get()))
        ax.set_xlabel(xlabel, fontsize=int(self.fs_label.get()))
        ax.set_ylabel(ylabel, fontsize=int(self.fs_label.get()))

        # Tick sizes
        ax.tick_params(axis="both", which="both", labelsize=int(self.fs_ticks.get()))

        # Minor ticks
        if self._minor_enabled(tabkey) and not is_matrix:
            ax.minorticks_on()
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
        else:
            # matrix heatmaps look messy with minor ticks
            if not is_matrix:
                ax.minorticks_off()

        # Integers for cycle x-axis, with an optional fixed tick step
        if x_is_cycle:
            step_var = {
                "pvc": getattr(self, "pvc_xtick_step_var", None),
                "evc": getattr(self, "evc_xtick_step_var", None),
            }.get(tabkey)
            step_val = None
            if step_var is not None:
                raw = step_var.get().strip()
                if raw:
                    try:
                        step_val = int(raw)
                    except ValueError:
                        step_val = None
            if step_val and step_val > 0:
                ax.xaxis.set_major_locator(MultipleLocator(step_val))
            else:
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        # Disable offset text
        self._disable_offset(ax)

        # Legend show/hide (uses cached handles/labels so it can be restored
        # even after being hidden, without needing to re-plot)
        legend_var = {
            "pvc": getattr(self, "pvc_legend_var", None),
            "evc": getattr(self, "evc_legend_var", None),
            "sc":  getattr(self, "sc_legend_var", None),
        }.get(tabkey)
        hl_attr = {"pvc": "_pvc_legend_hl", "evc": "_evc_legend_hl", "sc": "_sc_legend_hl"}.get(tabkey)

        if legend_var is not None:
            existing = ax.get_legend()
            if not legend_var.get():
                if existing is not None:
                    existing.remove()
            else:
                hl = getattr(self, hl_attr, None) if hl_attr else None
                if hl and hl[0]:
                    handles, labels = hl
                    ax.legend(handles, labels, loc="best")

        # Legend font size (if legend exists)
        leg = ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                t.set_fontsize(int(self.fs_legend.get()))

        canvas.draw_idle()

    def _add_style_controls(self, parent: ttk.Frame, tabkey: str, x_is_cycle: bool=False, show_legend: bool=False):
        """
        Adds:
        - Title/X/Y label entries (empty => auto title/labels)
        - Minor ticks checkbox
        - X tick step entry (only for cycle-axis tabs)
        - Show legend checkbox (+ merge-by-anion for evc)
        - Global font size spinboxes (shared)
        """
        # Label edit box
        box = ttk.LabelFrame(parent, text="Plot labels (live)")
        box.pack(fill=tk.X, pady=(8, 6))

        def row(lbl, var):
            r = ttk.Frame(box); r.pack(fill=tk.X, padx=6, pady=2)
            ttk.Label(r, text=lbl, width=8).pack(side=tk.LEFT)
            e = ttk.Entry(r, textvariable=var)
            e.pack(side=tk.LEFT, fill=tk.X, expand=True)

        title_var = {"pvc": self.pvc_title_var, "evc": self.evc_title_var, "sc": self.sc_title_var,
                     "corr": self.corr_title_var, "tcorr": self.tcorr_title_var}[tabkey]
        xlabel_var = {"pvc": self.pvc_xlabel_var, "evc": self.evc_xlabel_var, "sc": self.sc_xlabel_var,
                      "corr": self.corr_xlabel_var, "tcorr": self.tcorr_xlabel_var}[tabkey]
        ylabel_var = {"pvc": self.pvc_ylabel_var, "evc": self.evc_ylabel_var, "sc": self.sc_ylabel_var,
                      "corr": self.corr_ylabel_var, "tcorr": self.tcorr_ylabel_var}[tabkey]
        minor_var = {"pvc": self.pvc_minor_var, "evc": self.evc_minor_var, "sc": self.sc_minor_var,
                     "corr": self.corr_minor_var, "tcorr": self.tcorr_minor_var}[tabkey]

        row("Title", title_var)
        row("X label", xlabel_var)
        row("Y label", ylabel_var)

        ttk.Checkbutton(box, text="Minor ticks", variable=minor_var).pack(anchor="w", padx=6, pady=(4, 2))

        if x_is_cycle:
            xtick_var = {"pvc": self.pvc_xtick_step_var, "evc": self.evc_xtick_step_var}.get(tabkey)
            if xtick_var is not None:
                r = ttk.Frame(box); r.pack(fill=tk.X, padx=6, pady=2)
                ttk.Label(r, text="X tick step", width=10).pack(side=tk.LEFT)
                ttk.Entry(r, textvariable=xtick_var, width=8).pack(side=tk.LEFT)
                ttk.Label(r, text="(cycles; empty = auto)", foreground="#666").pack(side=tk.LEFT, padx=(6, 0))

        if show_legend:
            legend_var = {"pvc": self.pvc_legend_var, "evc": self.evc_legend_var, "sc": self.sc_legend_var}.get(tabkey)
            if legend_var is not None:
                ttk.Checkbutton(box, text="Show legend", variable=legend_var).pack(anchor="w", padx=6, pady=(4, 2))
            if tabkey == "evc":
                ttk.Checkbutton(
                    box, text="Merge legend by anion", variable=self.evc_legend_by_anion_var
                ).pack(anchor="w", padx=6, pady=(0, 6))

        # Font controls box (shared vars)
        fb = ttk.LabelFrame(parent, text="Font sizes (global; live)")
        fb.pack(fill=tk.X, pady=(0, 6))

        def spin(lbl, var):
            r = ttk.Frame(fb); r.pack(fill=tk.X, padx=6, pady=2)
            ttk.Label(r, text=lbl, width=8).pack(side=tk.LEFT)
            ttk.Spinbox(r, from_=6, to=40, width=6, textvariable=var).pack(side=tk.LEFT)

        spin("Title",  self.fs_title)
        spin("Labels", self.fs_label)
        spin("Ticks",  self.fs_ticks)
        spin("Legend", self.fs_legend)

    # ========== SECTION 8.2: LAYOUT ==========
    def _build_topbar(self):
        top = ttk.Frame(self)
        top.pack(side=tk.TOP, fill=tk.X, padx=10, pady=8)

        ttk.Button(top, text="Load CSV files…", command=self.load_files).pack(side=tk.LEFT)
        ttk.Button(top, text="Load a folder…", command=self.load_folder).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(top, text="↻ Reload last folder", command=self._auto_load_last_folder).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(top, text="Clear", command=self.clear_all).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(top, text="Save (tight)…", command=self.save_current_figure).pack(side=tk.LEFT, padx=(8, 0))


        ttk.Separator(top, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        ttk.Label(top, text="Filter: Anion").pack(side=tk.LEFT)
        self.filter_anion = ttk.Combobox(top, state="readonly", width=18)
        self.filter_anion.pack(side=tk.LEFT, padx=(6, 12))
        self.filter_anion.bind("<<ComboboxSelected>>", lambda e: self._refresh_filters_and_lists())

        ttk.Label(top, text="Filter: Experiment type").pack(side=tk.LEFT)
        self.filter_type = ttk.Combobox(top, state="readonly", width=18)
        self.filter_type.pack(side=tk.LEFT, padx=(6, 12))
        self.filter_type.bind("<<ComboboxSelected>>", lambda e: self._refresh_filters_and_lists())

        self.status_var = tk.StringVar(value="Load one or more CSV files to begin.")
        ttk.Label(top, textvariable=self.status_var).pack(side=tk.RIGHT)

    def _build_notebook(self):
        nb = ttk.Notebook(self)
        nb.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0, 10))
        self.nb = nb

        self.tab_params = self._build_tab_params_vs_cycles(nb)
        self.tab_exps = self._build_tab_exps_vs_cycles(nb)
        self.tab_scatter = self._build_tab_scatter(nb)
        self.tab_corr = self._build_tab_corr(nb)
        self.tab_tcorr = self._build_tab_target_corr(nb)

        nb.add(self.tab_params, text="Params vs Cycles")
        nb.add(self.tab_exps, text="Experiments vs Cycles")
        nb.add(self.tab_scatter, text="Scatter")
        nb.add(self.tab_corr, text="Correlation Matrix")
        nb.add(self.tab_tcorr, text="Target Correlations")

    def _make_plot_area(self, parent: ttk.Frame):
        fig = Figure(figsize=(6, 4), dpi=110, constrained_layout=True)
        ax = fig.add_subplot(111)
        canvas = FigureCanvasTkAgg(fig, master=parent)
        widget = canvas.get_tk_widget()
        toolbar = NavigationToolbar2Tk(canvas, parent)
        toolbar.update()
        widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        return fig, ax, canvas

    def _make_scrollable_left(self, parent: ttk.Frame, width: int = 320):
        """
        Left control panels have grown (labels, font sizes, tick step,
        legend options, favorites...) and can exceed the window height,
        which hides buttons like "Plot" below the visible area.
        This wraps the panel in a scrollable canvas so everything stays
        reachable (mouse wheel or the scrollbar) regardless of window size.
        Returns (container, inner_frame): pack the container, add widgets to inner_frame.
        """
        container = ttk.Frame(parent)
        canvas = tk.Canvas(container, highlightthickness=0, width=width)
        vsb = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        inner = ttk.Frame(canvas)

        inner.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=inner, anchor="nw")
        canvas.configure(yscrollcommand=vsb.set)

        canvas.pack(side=tk.LEFT, fill=tk.Y, expand=False)
        vsb.pack(side=tk.LEFT, fill=tk.Y)

        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        def _bind_wheel(_e=None):
            canvas.bind_all("<MouseWheel>", _on_mousewheel)

        def _unbind_wheel(_e=None):
            canvas.unbind_all("<MouseWheel>")

        canvas.bind("<Enter>", _bind_wheel)
        canvas.bind("<Leave>", _unbind_wheel)

        return container, inner

    # ========== SECTION 9: TABS ==========
    def _build_tab_params_vs_cycles(self, nb: ttk.Notebook):
        root = ttk.Frame(nb)
        left_container, left = self._make_scrollable_left(root)
        right = ttk.Frame(root)
        left_container.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        ttk.Label(left, text="Experiment").pack(anchor="w")
        self.pvc_exp_cb = ttk.Combobox(left, state="readonly", width=42)
        self.pvc_exp_cb.pack(anchor="w", pady=(0, 8))
        self.pvc_exp_cb.bind("<<ComboboxSelected>>", lambda e: self._pvc_on_experiment_changed())

        ttk.Label(left, text="Parameters (multi-select)").pack(anchor="w")
        self.pvc_param_tree = ParamTree(left)
        self.pvc_param_tree.pack(anchor="w", fill=tk.BOTH, expand=True)

        btnrow = ttk.Frame(left)
        btnrow.pack(fill=tk.X, pady=(6, 0))
        ttk.Button(btnrow, text="Select all", command=self.pvc_param_tree.select_all_leaves).pack(side=tk.LEFT)
        ttk.Button(btnrow, text="Clear", command=self.pvc_param_tree.clear_selection).pack(side=tk.LEFT, padx=(6, 0))
        ttk.Button(btnrow, text="★ Toggle favorite", command=self.pvc_param_tree.toggle_favorites_for_selection).pack(side=tk.LEFT, padx=(6, 0))

        self._add_style_controls(left, "pvc", x_is_cycle=True, show_legend=True)

        ttk.Button(left, text="Plot", command=self.plot_params_vs_cycles).pack(anchor="w", pady=(6, 0))

        self.pvc_fig, self.pvc_ax, self.pvc_canvas = self._make_plot_area(right)
        return root

    def _build_tab_exps_vs_cycles(self, nb: ttk.Notebook):
        root = ttk.Frame(nb)
        left_container, left = self._make_scrollable_left(root)
        right = ttk.Frame(root)
        left_container.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        ttk.Label(left, text="Parameter").pack(anchor="w")
        prow = ttk.Frame(left)
        prow.pack(anchor="w", fill=tk.X, pady=(0, 8))
        self.evc_param_cb = ttk.Combobox(prow, state="readonly", width=36)
        self.evc_param_cb.pack(side=tk.LEFT)
        self.evc_param_cb.bind(
            "<<ComboboxSelected>>",
            lambda e: (self._evc_update_fav_btn(), self._evc_schedule_replot())
        )
        self.evc_fav_btn = ttk.Button(prow, text="☆", width=3, command=self._evc_toggle_favorite_current)
        self.evc_fav_btn.pack(side=tk.LEFT, padx=(4, 0))

        ttk.Label(left, text="Experiments (multi-select)").pack(anchor="w")
        self.evc_exp_tree = AnionExperimentTree(left)
        self.evc_exp_tree.pack(anchor="w", fill=tk.BOTH, expand=True)
        # any change of the experiment selection (incl. clicking an anion header) re-plots
        self.evc_exp_tree.tree.bind("<<TreeviewSelect>>", lambda e: self._evc_schedule_replot(), add="+")

        btnrow = ttk.Frame(left)
        btnrow.pack(fill=tk.X, pady=(6, 0))
        ttk.Button(btnrow, text="Select all", command=self.evc_exp_tree.select_all_leaves).pack(side=tk.LEFT)
        ttk.Button(btnrow, text="Clear", command=self.evc_exp_tree.clear_selection).pack(side=tk.LEFT, padx=(6, 0))

        ttk.Checkbutton(
            left, text="Fixed Y range (same experiment type)", variable=self.evc_fixed_y_var
        ).pack(anchor="w", pady=(8, 0))
        ttk.Label(
            left, text="(most extreme value of this parameter among loaded files of the same experiment type)",
            foreground="#666", wraplength=300
        ).pack(anchor="w")

        self._add_style_controls(left, "evc", x_is_cycle=True, show_legend=True)

        ttk.Button(left, text="Plot", command=self.plot_experiments_vs_cycles).pack(anchor="w", pady=(6, 0))

        self.evc_fig, self.evc_ax, self.evc_canvas = self._make_plot_area(right)
        return root

    def _build_tab_scatter(self, nb: ttk.Notebook):
        root = ttk.Frame(nb)
        left_container, left = self._make_scrollable_left(root)
        right = ttk.Frame(root)
        left_container.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        ttk.Label(left, text="X parameter").pack(anchor="w")
        self.sc_x_cb = ttk.Combobox(left, state="readonly", width=42)
        self.sc_x_cb.pack(anchor="w", pady=(0, 8))

        ttk.Label(left, text="Y parameter").pack(anchor="w")
        self.sc_y_cb = ttk.Combobox(left, state="readonly", width=42)
        self.sc_y_cb.pack(anchor="w", pady=(0, 8))

        ttk.Label(left, text="Experiments (multi-select)").pack(anchor="w")
        self.sc_exp_tree = AnionExperimentTree(left)
        self.sc_exp_tree.pack(anchor="w", fill=tk.BOTH, expand=True)

        btnrow = ttk.Frame(left)
        btnrow.pack(fill=tk.X, pady=(6, 0))
        ttk.Button(btnrow, text="Select all", command=self.sc_exp_tree.select_all_leaves).pack(side=tk.LEFT)
        ttk.Button(btnrow, text="Clear", command=self.sc_exp_tree.clear_selection).pack(side=tk.LEFT, padx=(6, 0))

        ttk.Separator(left, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=8)

        ttk.Label(left, text="Point mode").pack(anchor="w")
        self.sc_mode = ttk.Combobox(
            left, state="readonly", width=28,
            values=["Per cycle points", "One point per experiment (mean)", "One point per experiment (last cycle)"]
        )
        self.sc_mode.pack(anchor="w", pady=(0, 8))
        self.sc_mode.set("Per cycle points")

        ttk.Label(left, text="Color by").pack(anchor="w")
        self.sc_color_by = ttk.Combobox(left, state="readonly", width=22, values=["Experiment", "Anion", "Experiment type"])
        self.sc_color_by.pack(anchor="w", pady=(0, 8))
        self.sc_color_by.set("Experiment")

        self._add_style_controls(left, "sc", x_is_cycle=False, show_legend=True)

        ttk.Button(left, text="Plot", command=self.plot_scatter).pack(anchor="w", pady=(6, 0))

        self.sc_fig, self.sc_ax, self.sc_canvas = self._make_plot_area(right)
        return root

    def _build_tab_corr(self, nb: ttk.Notebook):
        root = ttk.Frame(nb)
        left_container, left = self._make_scrollable_left(root)
        right = ttk.Frame(root)
        left_container.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        ttk.Label(left, text="Correlation mode").pack(anchor="w")
        self.corr_mode = ttk.Combobox(
            left, state="readonly", width=34,
            values=["Across experiments (aggregate cycles)", "Across all cycles (stacked)"]
        )
        self.corr_mode.pack(anchor="w", pady=(0, 8))
        self.corr_mode.set("Across experiments (aggregate cycles)")

        ttk.Label(left, text="Aggregation (for 'Across experiments')").pack(anchor="w")
        self.corr_agg = ttk.Combobox(left, state="readonly", width=34, values=["mean", "median", "last"])
        self.corr_agg.pack(anchor="w", pady=(0, 8))
        self.corr_agg.set("mean")

        ttk.Label(left, text="Method").pack(anchor="w")
        self.corr_method = ttk.Combobox(left, state="readonly", width=34, values=["pearson", "spearman"])
        self.corr_method.pack(anchor="w", pady=(0, 8))
        self.corr_method.set("pearson")

        ttk.Separator(left, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=8)

        ttk.Label(left, text="Parameters (multi-select)").pack(anchor="w")
        self.corr_param_tree = ParamTree(left)
        self.corr_param_tree.pack(anchor="w", fill=tk.BOTH, expand=False)

        btnrow = ttk.Frame(left)
        btnrow.pack(fill=tk.X, pady=(6, 0))
        ttk.Button(btnrow, text="Select all", command=self.corr_param_tree.select_all_leaves).pack(side=tk.LEFT)
        ttk.Button(btnrow, text="Clear", command=self.corr_param_tree.clear_selection).pack(side=tk.LEFT, padx=(6, 0))
        ttk.Button(btnrow, text="★ Toggle favorite", command=self.corr_param_tree.toggle_favorites_for_selection).pack(side=tk.LEFT, padx=(6, 0))

        # Style controls (note: minor ticks not recommended for matrix, but checkbox exists; we ignore it in renderer)
        self._add_style_controls(left, "corr", x_is_cycle=False, show_legend=False)

        ttk.Button(left, text="Compute", command=self.plot_correlation_matrix).pack(anchor="w", pady=(6, 6))

        ttk.Separator(left, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=8)

        flt = ttk.LabelFrame(left, text="Filter displayed correlations")
        flt.pack(fill=tk.X, pady=(0, 8))

        ttk.Checkbutton(
            flt, text="Use absolute values (|r|)",
            variable=self.corr_abs_var,
            command=self._corr_on_filter_change
        ).pack(anchor="w", padx=6, pady=(4, 2))

        r1 = ttk.Frame(flt); r1.pack(fill=tk.X, padx=6, pady=(2, 2))
        ttk.Label(r1, text="Low").pack(side=tk.LEFT)
        self.corr_low_scale = tk.Scale(
            r1, from_=0.0, to=1.0, resolution=0.01, orient="horizontal",
            variable=self.corr_low_var, command=lambda _v: self._corr_on_filter_change()
        )
        self.corr_low_scale.pack(side=tk.RIGHT, fill=tk.X, expand=True)

        r2 = ttk.Frame(flt); r2.pack(fill=tk.X, padx=6, pady=(2, 4))
        ttk.Label(r2, text="High").pack(side=tk.LEFT)
        self.corr_high_scale = tk.Scale(
            r2, from_=0.0, to=1.0, resolution=0.01, orient="horizontal",
            variable=self.corr_high_var, command=lambda _v: self._corr_on_filter_change()
        )
        self.corr_high_scale.pack(side=tk.RIGHT, fill=tk.X, expand=True)

        self.corr_range_lbl = ttk.Label(flt, text="Showing: |r| in [0.70, 1.00]")
        self.corr_range_lbl.pack(anchor="w", padx=6, pady=(0, 6))

        pairs_box = ttk.LabelFrame(left, text="Pairs (in range)")
        pairs_box.pack(fill=tk.BOTH, expand=True, pady=(0, 8))

        self.corr_pairs_lb = tk.Listbox(pairs_box, width=44, height=10, exportselection=False)
        self.corr_pairs_lb.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        sb = ttk.Scrollbar(pairs_box, orient="vertical", command=self.corr_pairs_lb.yview)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.corr_pairs_lb.config(yscrollcommand=sb.set)

        self.corr_fig, self.corr_ax, self.corr_canvas = self._make_plot_area(right)
        return root

    def _build_tab_target_corr(self, nb: ttk.Notebook):
        root = ttk.Frame(nb)
        left_container, left = self._make_scrollable_left(root)
        right = ttk.Frame(root)
        left_container.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        ttk.Label(left, text="Mode (observations)").pack(anchor="w")
        self.tcorr_mode = ttk.Combobox(
            left, state="readonly", width=34,
            values=["Across experiments (aggregate cycles)", "Across all cycles (stacked)"]
        )
        self.tcorr_mode.pack(anchor="w", pady=(0, 8))
        self.tcorr_mode.set("Across all cycles (stacked)")

        ttk.Label(left, text="Aggregation (for 'Across experiments')").pack(anchor="w")
        self.tcorr_agg = ttk.Combobox(left, state="readonly", width=34, values=["mean", "median", "last"])
        self.tcorr_agg.pack(anchor="w", pady=(0, 8))
        self.tcorr_agg.set("mean")

        ttk.Label(left, text="Method").pack(anchor="w")
        self.tcorr_method = ttk.Combobox(left, state="readonly", width=34, values=["pearson", "spearman"])
        self.tcorr_method.pack(anchor="w", pady=(0, 8))
        self.tcorr_method.set("pearson")

        ttk.Separator(left, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=8)

        ttk.Label(left, text="Target parameter").pack(anchor="w")
        self.tcorr_target_cb = ttk.Combobox(left, state="readonly", width=34, values=[])
        self.tcorr_target_cb.pack(anchor="w", pady=(0, 8))

        ttk.Checkbutton(
            left, text="Sort by |r| (strength)",
            variable=self.tcorr_abs_var,
            command=self._tcorr_render_last
        ).pack(anchor="w", pady=(0, 6))

        minrow = ttk.LabelFrame(left, text="Display filter")
        minrow.pack(fill=tk.X, pady=(0, 8))
        ttk.Label(minrow, text="Min |r|").pack(anchor="w", padx=6, pady=(4, 0))
        self.tcorr_min_scale = tk.Scale(
            minrow, from_=0.0, to=1.0, resolution=0.01, orient="horizontal",
            variable=self.tcorr_minabs_var, command=lambda _v: self._tcorr_render_last()
        )
        self.tcorr_min_scale.pack(fill=tk.X, padx=6, pady=(0, 6))

        topnrow = ttk.Frame(left)
        topnrow.pack(fill=tk.X, pady=(0, 8))
        ttk.Label(topnrow, text="Top N").pack(side=tk.LEFT)
        self.tcorr_topn_spin = tk.Spinbox(
            topnrow, from_=5, to=200, width=6,
            textvariable=self.tcorr_topn_var, command=self._tcorr_render_last
        )
        self.tcorr_topn_spin.pack(side=tk.LEFT, padx=(8, 0))
        ttk.Label(topnrow, text="(bars shown)").pack(side=tk.LEFT, padx=(6, 0))

        self._add_style_controls(left, "tcorr", x_is_cycle=False, show_legend=False)

        ttk.Button(left, text="Compute", command=self.plot_target_correlations).pack(anchor="w", pady=(6, 6))

        box = ttk.LabelFrame(left, text="Target vs others (in range)")
        box.pack(fill=tk.BOTH, expand=True)
        self.tcorr_lb = tk.Listbox(box, width=44, height=14, exportselection=False)
        self.tcorr_lb.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb = ttk.Scrollbar(box, orient="vertical", command=self.tcorr_lb.yview)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.tcorr_lb.config(yscrollcommand=sb.set)

        self.tcorr_fig, self.tcorr_ax, self.tcorr_canvas = self._make_plot_area(right)
        return root

    # ========== SECTION 10: LOAD / CLEAR ==========
    def load_files(self):
        paths = filedialog.askopenfilenames(
            title="Select CSV files",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not paths:
            return
        self._load_paths(list(paths))

    def load_folder(self):
        folder = filedialog.askdirectory(title="Select a folder containing CSV files")
        if not folder:
            return
        paths = [
            os.path.join(folder, fn) for fn in os.listdir(folder)
            if fn.lower().endswith(".csv")
        ]
        if not paths:
            messagebox.showinfo("No CSV files", "No .csv files found in the selected folder.")
            return
        self._load_paths(paths)

        cfg = load_config()
        cfg["last_folder"] = folder
        save_config(cfg)

    def _auto_load_last_folder(self):
        """Silently (re)load CSVs from the last folder used with 'Load a folder…',
        if any. Runs once automatically on startup, and can be re-triggered
        manually from the toolbar."""
        cfg = load_config()
        folder = cfg.get("last_folder")
        if not folder or not os.path.isdir(folder):
            return
        paths = [
            os.path.join(folder, fn) for fn in os.listdir(folder)
            if fn.lower().endswith(".csv")
        ]
        if not paths:
            return
        self._load_paths(paths)
        self.status_var.set(
            f"Auto-loaded {len(paths)} file(s) from last folder: {folder}"
        )

    def _load_paths(self, paths: list):
        loaded = 0
        skipped = 0
        for p in paths:
            try:
                df = safe_read_csv(p)
                exp = Experiment(p, df)
                self.experiments[exp.stem] = exp
                loaded += 1
            except Exception as e:
                skipped += 1
                print(f"[WARN] Failed to load {p}: {e}")

        self.status_var.set(f"Loaded {loaded} file(s), skipped {skipped}. Total experiments: {len(self.experiments)}")
        self._refresh_filters_and_lists()

    def clear_all(self):
        self.experiments.clear()
        self.status_var.set("Cleared. Load one or more CSV files to begin.")
        self._refresh_filters_and_lists()

        for ax, canvas in [
            (getattr(self, "pvc_ax", None), getattr(self, "pvc_canvas", None)),
            (getattr(self, "evc_ax", None), getattr(self, "evc_canvas", None)),
            (getattr(self, "sc_ax", None), getattr(self, "sc_canvas", None)),
            (getattr(self, "corr_ax", None), getattr(self, "corr_canvas", None)),
            (getattr(self, "tcorr_ax", None), getattr(self, "tcorr_canvas", None)),
        ]:
            if ax is not None:
                ax.clear()
                if canvas is not None:
                    canvas.draw_idle()
                    
    def save_current_figure(self):
        tab = self.nb.tab(self.nb.select(), "text")
    
        mapping = {
            "Params vs Cycles": ("pvc_fig",),
            "Experiments vs Cycles": ("evc_fig",),
            "Scatter": ("sc_fig",),
            "Correlation Matrix": ("corr_fig",),
            "Target Correlations": ("tcorr_fig",),
        }
        fig_name = mapping.get(tab, ("pvc_fig",))[0]
        fig = getattr(self, fig_name, None)
        if fig is None:
            return
    
        suggested = ""
        try:
            suggested = re.sub(r'[\\/:*?"<>|\r\n]+', " ", fig.axes[0].get_title()).strip(" ._")[:100]
        except Exception:
            pass
        path = filedialog.asksaveasfilename(
            defaultextension=".png",
            initialfile=suggested or "figure",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg"), ("All files", "*.*")]
        )
        if not path:
            return
    
        fig.savefig(path, bbox_inches="tight", pad_inches=0.25, dpi=300)


    # ========== SECTION 11: FILTERING + LISTS ==========
    def _get_filtered_experiment_stems(self) -> list:
        stems = sorted(self.experiments.keys())
        fa = self.filter_anion.get()
        ft = self.filter_type.get()

        def ok(stem: str) -> bool:
            exp = self.experiments[stem]
            if fa and fa != "All" and exp.anion != fa:
                return False
            if ft and ft != "All" and exp.experimenttype != ft:
                return False
            return True

        return [s for s in stems if ok(s)]

    def _all_numeric_params(self, stems: list) -> list:
        s = set()
        for stem in stems:
            s |= set(self.experiments[stem].numeric_cols)
        return sorted(s)

    def _refresh_filters_and_lists(self):
        anions = sorted({e.anion for e in self.experiments.values()}) if self.experiments else []
        types = sorted({e.experimenttype for e in self.experiments.values()}) if self.experiments else []

        current_anion = self.filter_anion.get() or "All"
        current_type = self.filter_type.get() or "All"

        set_combobox_values(self.filter_anion, ["All"] + anions, keep=current_anion)
        set_combobox_values(self.filter_type, ["All"] + types, keep=current_type)

        filtered = self._get_filtered_experiment_stems()

        set_combobox_values(self.pvc_exp_cb, filtered, keep=self.pvc_exp_cb.get() if hasattr(self, "pvc_exp_cb") else None)
        self._pvc_on_experiment_changed()

        all_params = self._all_numeric_params(filtered)

        # Columns the parameter trees don't list yet -> "Other columns"
        extra = set(self._all_numeric_params(sorted(self.experiments))) - _TREE_PARAMS
        if extra != LOADED_EXTRA_PARAMS:
            LOADED_EXTRA_PARAMS.clear()
            LOADED_EXTRA_PARAMS.update(extra)
            for inst in list(_PARAM_TREE_INSTANCES):
                try:
                    inst._build_tree()
                except Exception:
                    pass

        if hasattr(self, "evc_param_cb"):
            self._evc_set_param_values(all_params, keep=self._evc_current_param() or None)
        if hasattr(self, "evc_exp_tree"):
            self.evc_exp_tree.set_experiments(filtered, self.experiments)
            self._evc_schedule_replot()  # filters / loaded files changed -> redraw

        set_combobox_values(self.sc_x_cb, all_params, keep=self.sc_x_cb.get() if hasattr(self, "sc_x_cb") else None)
        set_combobox_values(self.sc_y_cb, all_params, keep=self.sc_y_cb.get() if hasattr(self, "sc_y_cb") else None)
        if hasattr(self, "sc_exp_tree"):
            self.sc_exp_tree.set_experiments(filtered, self.experiments)

        if hasattr(self, "tcorr_target_cb"):
            keep = self.tcorr_target_cb.get() if self.tcorr_target_cb.get() in all_params else None
            set_combobox_values(self.tcorr_target_cb, all_params, keep=keep)

    def _pvc_on_experiment_changed(self):
        stem = self.pvc_exp_cb.get()
        if not stem or stem not in self.experiments:
            return

    # ---- Experiments vs Cycles: favorite-parameter dropdown helpers ----
    def _evc_current_param(self) -> str:
        """Real parameter name from the combobox, stripping the '★ ' favorite marker if present."""
        val = self.evc_param_cb.get()
        if val.startswith("★ "):
            return val[2:]
        return val

    def _evc_set_param_values(self, all_params: list, keep: str = None):
        """Fill the parameter dropdown with favorites first (starred), then the rest, A-Z."""
        favs = sorted(p for p in all_params if p in FAVORITE_PARAMS)
        others = sorted(p for p in all_params if p not in FAVORITE_PARAMS)
        display = [f"★ {p}" for p in favs] + others

        keep_display = None
        if keep:
            keep_display = f"★ {keep}" if keep in FAVORITE_PARAMS else keep

        set_combobox_values(self.evc_param_cb, display, keep=keep_display)
        self._evc_update_fav_btn()

    def _evc_update_fav_btn(self, event=None):
        if not hasattr(self, "evc_fav_btn"):
            return
        p = self._evc_current_param()
        self.evc_fav_btn.config(text=("★" if p and p in FAVORITE_PARAMS else "☆"))

    def _evc_toggle_favorite_current(self):
        p = self._evc_current_param()
        if not p:
            return
        if p in FAVORITE_PARAMS:
            FAVORITE_PARAMS.discard(p)
        else:
            FAVORITE_PARAMS.add(p)
        save_favorite_params(FAVORITE_PARAMS)

        all_params = self._all_numeric_params(self._get_filtered_experiment_stems())
        self._evc_set_param_values(all_params, keep=p)

        # keep the parameter trees (Params vs Cycles / Correlation Matrix) in sync
        for inst in list(_PARAM_TREE_INSTANCES):
            try:
                inst._build_tree()
            except Exception:
                pass

    # ========== SECTION 12: PLOTTING ==========
    def plot_params_vs_cycles(self):
        stem = self.pvc_exp_cb.get()
        if not stem or stem not in self.experiments:
            messagebox.showwarning("No experiment", "Select an experiment.")
            return

        params = self.pvc_param_tree.get_selected_params() if hasattr(self, "pvc_param_tree") else []
        if not params:
            messagebox.showwarning("No parameters", "Select one or more parameters.")
            return

        exp = self.experiments[stem]
        df = exp.df

        self.pvc_ax.clear()
        x = df["cycle"].values.astype(int)  # force integer cycle

        plotted = 0
        for p in params:
            if p not in df.columns:
                continue
            y = pd.to_numeric(df[p], errors="coerce").values
            m = ~np.isnan(y)
            if not np.any(m):
                continue
            # dots only:
            self.pvc_ax.scatter(x[m], y[m], label=label_with_unit(p), s=18, alpha=0.9)
            plotted += 1

        if plotted == 0:
            messagebox.showinfo("Nothing to plot", "Selected parameters not found / empty in this experiment.")
            return

        # Auto title/labels
        self._auto["pvc"]["title"] = exp.label()
        self._auto["pvc"]["xlabel"] = "Cycle"
        self._auto["pvc"]["ylabel"] = "Value"

        self._pvc_legend_hl = self.pvc_ax.get_legend_handles_labels()
        self.pvc_ax.grid(True, alpha=0.25)
        
        self.pvc_canvas.draw_idle()

        self._apply_styles_all()

    def _evc_auto_title(self, anions, types, param: str) -> str:
        """Anion name(s) + experiment-type wording, e.g. "Carbonate - activation without urea".
        The type comes from the 'Filter: Experiment type' box; if that is 'All', from the
        experiments being plotted."""
        names = sorted({a for a in anions if a and a != "Unknown"})
        flt = self.filter_type.get()
        if flt and flt != "All":
            type_list = [flt]
        else:
            type_list = sorted({t for t in types if t and t != "Unknown"})
        descs = []
        for t in type_list:
            d = describe_experiment_type(t)
            if d and d not in descs:
                descs.append(d)
        parts = []
        if names:
            parts.append(", ".join(names))
        if descs:
            parts.append(" / ".join(descs))
        return " - ".join(parts) if parts else f"{param} vs Cycle"

    def _evc_schedule_replot(self, *_):
        """Debounced: many selection events can fire at once (e.g. clicking an anion header)."""
        if self._evc_replot_id is not None:
            try:
                self.after_cancel(self._evc_replot_id)
            except Exception:
                pass
        self._evc_replot_id = self.after(60, self._evc_auto_replot)

    def _evc_auto_replot(self):
        self._evc_replot_id = None
        param = self._evc_current_param()
        stems = self.evc_exp_tree.get_selected_stems() if hasattr(self, "evc_exp_tree") else []
        if not param or not stems:
            # nothing selected -> empty graph (no pop-ups for automatic redraws)
            self.evc_ax.clear()
            self._evc_legend_hl = None
            self.evc_canvas.draw_idle()
            return
        self.plot_experiments_vs_cycles(silent=True)

    def _evc_global_ylim(self, param: str, exp_types=None, pad_frac: float = 0.05):
        """(ymin, ymax) of `param` over all loaded experiments whose experiment type is in
        `exp_types` (None = every experiment), independent of the anion filter and of the
        current selection, plus a small margin. Returns None if there are no finite values.
        Restricting to one experiment type keeps different measurement types from stretching each other's axis."""
        lo, hi = np.inf, -np.inf
        for exp in self.experiments.values():
            if exp_types is not None and exp.experimenttype not in exp_types:
                continue
            if param not in exp.df.columns:
                continue
            y = pd.to_numeric(exp.df[param], errors="coerce").to_numpy(dtype=float)
            y = y[np.isfinite(y)]
            if y.size == 0:
                continue
            lo = min(lo, float(y.min()))
            hi = max(hi, float(y.max()))
        if not (np.isfinite(lo) and np.isfinite(hi)):
            return None
        span = hi - lo
        pad = span * pad_frac if span > 0 else (abs(hi) * 0.05 or 1.0)
        return lo - pad, hi + pad

    def plot_experiments_vs_cycles(self, silent: bool = False):
        param = self._evc_current_param()
        if not param:
            messagebox.showwarning("No parameter", "Select a parameter.")
            return

        stems = self.evc_exp_tree.get_selected_stems() if hasattr(self, "evc_exp_tree") else []
        if not stems:
            messagebox.showwarning("No experiments", "Select one or more experiments.")
            return

        self.evc_ax.clear()

        anions = []
        for stem in stems:
            exp = self.experiments.get(stem)
            if exp is not None:
                anions.append(getattr(exp, "anion", ""))
        color_map = group_color_map(anions)

        plotted = 0
        plotted_anions = []
        plotted_types = set()
        for stem in stems:
            exp = self.experiments.get(stem)
            if exp is None:
                continue
            df = exp.df
            if param not in df.columns:
                continue

            x = df["cycle"].values.astype(int)
            y = pd.to_numeric(df[param], errors="coerce").values
            m = ~np.isnan(y)
            if not np.any(m):
                continue

            anion = getattr(exp, "anion", "")
            c = color_map.get(anion, None)
            label = f"{anion} | {stem}" if anion else stem

            self.evc_ax.plot(x[m], y[m], label=label, color=c, alpha=0.95)
            plotted += 1
            if anion not in plotted_anions:
                plotted_anions.append(anion)
            plotted_types.add(getattr(exp, "experimenttype", ""))

        if plotted == 0:
            if not silent:
                messagebox.showinfo("Nothing to plot", "Selected parameter not found in the selected experiments.")
            self.evc_canvas.draw_idle()
            return

        # Default title = anion(s) + wording for the experiment type (see EXPERIMENT_TYPE_DESCRIPTIONS)
        auto_title = self._evc_auto_title(plotted_anions, plotted_types, param)
        self._auto["evc"]["title"] = auto_title
        self._auto["evc"]["xlabel"] = "Cycle"
        self._auto["evc"]["ylabel"] = y_label_for_param(param)

        if bool(self.evc_fixed_y_var.get()):
            _ylim = self._evc_global_ylim(param, plotted_types)
            if _ylim is not None:
                self.evc_ax.set_ylim(*_ylim)

        # Legend: one entry per anion (merged), or one entry per experiment
        if bool(self.evc_legend_by_anion_var.get()):
            handles = [
                Line2D([0], [0], color=color_map.get(a, "C0"), label=(a if a else "Unknown"))
                for a in sorted(plotted_anions)
            ]
            labels = [h.get_label() for h in handles]
        else:
            handles, labels = self.evc_ax.get_legend_handles_labels()

        self._evc_legend_hl = (handles, labels)

        # Put the auto title into the Title textbox so it is visible and editable. If the user
        # typed their own text, keep it; if the box is empty or still holds the previous auto
        # title, refresh it for the new plot.
        cur_title = self.evc_title_var.get().strip()
        if cur_title == "" or cur_title == self._evc_last_auto_title:
            self.evc_title_var.set(auto_title)
        self._evc_last_auto_title = auto_title

        self.evc_ax.grid(True, alpha=0.25)
        
        self.evc_canvas.draw_idle()

        self._apply_styles_all()

    def plot_scatter(self):
        x_param = self.sc_x_cb.get()
        y_param = self.sc_y_cb.get()
        if not x_param or not y_param:
            messagebox.showwarning("Select parameters", "Choose X and Y parameters.")
            return
    
        stems = self.sc_exp_tree.get_selected_stems() if hasattr(self, "sc_exp_tree") else []
        if not stems:
            messagebox.showwarning("No experiments", "Select one or more experiments.")
            return
    
        mode = self.sc_mode.get()
        color_by = self.sc_color_by.get()
    
        self.sc_ax.clear()
    
        # Remove existing colorbar if you previously added one (optional safety)
        if getattr(self, "_sc_cbar", None) is not None:
            try:
                self._sc_cbar.remove()
            except Exception:
                pass
            self._sc_cbar = None
    
        def group_key(exp: Experiment) -> str:
            if color_by == "Anion":
                return exp.anion
            if color_by == "Experiment type":
                return exp.experimenttype
            return exp.stem
    
        # key -> list of (xs, ys, cycles_or_none)
        groups = {}
    
        for stem in stems:
            exp = self.experiments.get(stem)
            if exp is None:
                continue
            df = exp.df
            if x_param not in df.columns or y_param not in df.columns:
                continue
    
            x = pd.to_numeric(df[x_param], errors="coerce")
            y = pd.to_numeric(df[y_param], errors="coerce")
            cyc = df["cycle"]
    
            if mode == "Per cycle points":
                m = x.notna() & y.notna()
                if not m.any():
                    continue
                xs = x[m].to_numpy()
                ys = y[m].to_numpy()
                cs = cyc[m].to_numpy().astype(float)  # cycles for alpha gradient
    
            elif mode == "One point per experiment (mean)":
                xs = np.array([x.mean(skipna=True)])
                ys = np.array([y.mean(skipna=True)])
                cs = None
    
            else:  # "One point per experiment (last cycle)"
                m = x.notna() & y.notna()
                if not m.any():
                    continue
                last_idx = np.where(m.to_numpy())[0][-1]
                xs = np.array([x.iloc[last_idx]])
                ys = np.array([y.iloc[last_idx]])
                cs = None
    
            key = group_key(exp)
            groups.setdefault(key, []).append((xs, ys, cs))
    
        if not groups:
            messagebox.showinfo("Nothing to plot", "No matching data for the chosen parameters.")
            return
    
        # Group colors (keep existing logic)
        color_map = group_color_map(list(groups.keys()))
    
        # Dot size / alpha (fallbacks if not implemented)
        dot_size = float(self.dot_size.get()) if hasattr(self, "dot_size") else 18.0
        alpha_max = float(self.dot_alpha.get()) if hasattr(self, "dot_alpha") else 0.9
        alpha_min = 0.15  # first cycles more transparent
    
        for key, series_list in groups.items():
            xs_all = np.concatenate([t[0] for t in series_list])
            ys_all = np.concatenate([t[1] for t in series_list])
    
            base = color_map.get(key, "C0")
            r, g, b, _ = mcolors.to_rgba(base)
    
            if mode == "Per cycle points":
                cs_all = np.concatenate([t[2] for t in series_list]).astype(float)
    
                cmin = np.nanmin(cs_all)
                cmax = np.nanmax(cs_all)
                if not np.isfinite(cmin) or not np.isfinite(cmax) or cmax == cmin:
                    tnorm = np.ones_like(cs_all)
                else:
                    tnorm = (cs_all - cmin) / (cmax - cmin)
    
                alphas = alpha_min + tnorm * (alpha_max - alpha_min)
    
                colors = np.column_stack([
                    np.full_like(alphas, r),
                    np.full_like(alphas, g),
                    np.full_like(alphas, b),
                    alphas
                ])
    
                self.sc_ax.scatter(xs_all, ys_all, c=colors, s=dot_size, label=key)
    
            else:
                self.sc_ax.scatter(xs_all, ys_all, color=base, s=dot_size, alpha=alpha_max, label=key)
    
        # Titles/labels via your auto system
        self._auto["sc"]["title"] = f"Scatter: {y_param} vs {x_param}"
        self._auto["sc"]["xlabel"] = y_label_for_param(x_param)
        self._auto["sc"]["ylabel"] = y_label_for_param(y_param)
    
        handles = []
        labels = []
        for key in groups.keys():
            base = color_map.get(key, "C0")
            handles.append(
                Line2D(
                    [0], [0],
                    marker="o",
                    linestyle="None",
                    markersize=max(4.0, (dot_size ** 0.5)),  # rough mapping from scatter area->legend marker size
                    markerfacecolor=base,
                    markeredgecolor=base,
                    alpha=1.0
                )
            )
            labels.append(key)
    
        self._sc_legend_hl = (handles, labels)

        self.sc_ax.grid(True, alpha=0.25)
    
        self.sc_canvas.draw_idle()
        self._apply_styles_all()


    # ========== SECTION 13: CORRELATION MATRIX ==========
    def _corr_on_filter_change(self):
        low = float(self.corr_low_var.get())
        high = float(self.corr_high_var.get())
        if low > high:
            self.corr_high_var.set(low)
            high = low

        use_abs = bool(self.corr_abs_var.get())

        if use_abs:
            self.corr_low_scale.config(from_=0.0, to=1.0)
            self.corr_high_scale.config(from_=0.0, to=1.0)
            low = max(0.0, min(1.0, low))
            high = max(0.0, min(1.0, high))
            self.corr_low_var.set(low)
            self.corr_high_var.set(high)
            if hasattr(self, "corr_range_lbl"):
                self.corr_range_lbl.config(text=f"Showing: |r| in [{low:.2f}, {high:.2f}]")
        else:
            self.corr_low_scale.config(from_=-1.0, to=1.0)
            self.corr_high_scale.config(from_=-1.0, to=1.0)
            low = max(-1.0, min(1.0, low))
            high = max(-1.0, min(1.0, high))
            self.corr_low_var.set(low)
            self.corr_high_var.set(high)
            if hasattr(self, "corr_range_lbl"):
                self.corr_range_lbl.config(text=f"Showing: r in [{low:.2f}, {high:.2f}]")

        self._corr_render_last()

    def _corr_render_last(self):
        if self._corr_last_raw is None:
            return

        corr_raw = self._corr_last_raw.copy()
        use_abs = bool(self.corr_abs_var.get())
        low = float(self.corr_low_var.get())
        high = float(self.corr_high_var.get())

        corr_for_filter = corr_raw.abs() if use_abs else corr_raw

        vals = corr_raw.values.astype(float)
        filt_vals = corr_for_filter.values.astype(float)

        n = vals.shape[0]
        mask = np.zeros((n, n), dtype=bool)
        mask |= np.isnan(vals) | np.isnan(filt_vals)
        mask |= (filt_vals < low) | (filt_vals > high)
        mask |= np.eye(n, dtype=bool)

        display_vals = corr_raw.abs().values if use_abs else corr_raw.values
        mvals = np.ma.array(display_vals, mask=mask)

        self.corr_ax.clear()

        if getattr(self, "_corr_cbar", None) is not None:
            try:
                self._corr_cbar.remove()
            except Exception:
                pass
            self._corr_cbar = None

        # matplotlib.cm.get_cmap was removed in matplotlib 3.9
        try:
            cmap = matplotlib.colormaps["viridis"].copy()
        except AttributeError:  # matplotlib < 3.5
            cmap = matplotlib.cm.get_cmap("viridis").copy()
        cmap.set_bad(alpha=0.0)

        if use_abs:
            im = self.corr_ax.imshow(mvals, aspect="auto", vmin=0, vmax=1, cmap=cmap)
        else:
            im = self.corr_ax.imshow(mvals, aspect="auto", vmin=-1, vmax=1, cmap=cmap)

        self._corr_cbar = self.corr_fig.colorbar(im, ax=self.corr_ax, fraction=0.046, pad=0.04)

        cols = list(corr_raw.columns)
        self.corr_ax.set_xticks(range(len(cols)))
        self.corr_ax.set_yticks(range(len(cols)))
        self.corr_ax.set_xticklabels(cols, rotation=90, fontsize=7)
        self.corr_ax.set_yticklabels(cols, fontsize=7)

        self._auto["corr"]["title"] = (self._corr_last_title or "Correlation Matrix")
        self._auto["corr"]["xlabel"] = ""
        self._auto["corr"]["ylabel"] = ""

        
        self.corr_canvas.draw_idle()

        if hasattr(self, "corr_pairs_lb"):
            self.corr_pairs_lb.delete(0, tk.END)

            pairs = []
            for i in range(n):
                for j in range(i + 1, n):
                    r = corr_raw.iat[i, j]
                    if pd.isna(r):
                        continue
                    key = abs(r) if use_abs else r
                    if low <= key <= high:
                        pairs.append((abs(r), r, cols[i], cols[j]))

            pairs.sort(key=lambda t: t[0], reverse=True)

            for absr, r, a, b in pairs[:500]:
                if use_abs:
                    self.corr_pairs_lb.insert(tk.END, f"{a} ↔ {b}   r={r:+.3f}  |r|={absr:.3f}")
                else:
                    self.corr_pairs_lb.insert(tk.END, f"{a} ↔ {b}   r={r:+.3f}")

        self._apply_styles_all()

    def plot_correlation_matrix(self):
        mode = self.corr_mode.get()
        agg = self.corr_agg.get()
        method = self.corr_method.get()

        params = self.corr_param_tree.get_selected_params() if hasattr(self, "corr_param_tree") else []
        if not params:
            messagebox.showwarning("No parameters", "Select one or more parameters.")
            return

        stems = self._get_filtered_experiment_stems()
        if len(stems) < 1:
            messagebox.showwarning("No experiments", "Load and select at least 1 experiment (after filtering).")
            return

        if mode == "Across experiments (aggregate cycles)":
            if len(stems) < 2:
                messagebox.showwarning(
                    "Need more experiments",
                    "This mode correlates across experiments (one row per experiment), so you need at least 2 experiments."
                )
                return

            rows = []
            row_labels = []
            for stem in stems:
                exp = self.experiments[stem]
                df = exp.df
                row = {}
                for p in params:
                    if p not in df.columns:
                        row[p] = np.nan
                        continue
                    s = pd.to_numeric(df[p], errors="coerce")
                    if agg == "median":
                        row[p] = s.median(skipna=True)
                    elif agg == "last":
                        s2 = s.dropna()
                        row[p] = float(s2.iloc[-1]) if len(s2) else np.nan
                    else:
                        row[p] = s.mean(skipna=True)
                rows.append(row)
                row_labels.append(exp.stem)

            mat = pd.DataFrame(rows, index=row_labels)

        else:
            stack = []
            for stem in stems:
                exp = self.experiments[stem]
                df = exp.df
                cols = [p for p in params if p in df.columns]
                if not cols:
                    continue
                sub = df[cols].copy()
                stack.append(sub)

            if not stack:
                messagebox.showwarning("No data", "None of the selected parameters exist in the filtered experiments.")
                return

            mat = pd.concat(stack, axis=0, ignore_index=True)

        mat = mat.apply(pd.to_numeric, errors="coerce")
        mat = mat.dropna(axis=0, how="all").dropna(axis=1, how="all")

        if mat.shape[0] < 2:
            messagebox.showwarning("Not enough data", "Need at least 2 observations (rows/cycles) to compute correlation.")
            return

        corr = mat.corr(method=method, min_periods=2)

        if corr.size == 0 or np.isnan(corr.values).all():
            messagebox.showwarning(
                "Correlation undefined",
                "Correlation is all NaN (not enough valid paired data, or parameters are constant / missing)."
            )
            return

        self._corr_last_raw = corr
        self._corr_last_title = f"Correlation ({method}) | {mode}"
        self._corr_on_filter_change()

    # ========== SECTION 14: TARGET CORRELATIONS ==========
    def plot_target_correlations(self):
        target = self.tcorr_target_cb.get().strip() if hasattr(self, "tcorr_target_cb") else ""
        if not target:
            messagebox.showwarning("No target", "Choose a target parameter.")
            return

        stems = self._get_filtered_experiment_stems()
        if len(stems) < 1:
            messagebox.showwarning("No experiments", "Load and select at least 1 experiment (after filtering).")
            return

        mode = self.tcorr_mode.get()
        agg = self.tcorr_agg.get()
        method = self.tcorr_method.get()

        all_params = self._all_numeric_params(stems)
        if target not in all_params:
            messagebox.showwarning("Target not available", f"Target parameter '{target}' is not available in the current filtered set.")
            return

        if mode == "Across experiments (aggregate cycles)":
            if len(stems) < 2:
                messagebox.showwarning(
                    "Need more experiments",
                    "This mode correlates across experiments (one row per experiment), so you need at least 2 experiments."
                )
                return

            rows = []
            for stem in stems:
                exp = self.experiments[stem]
                df = exp.df
                row = {}
                for p in exp.numeric_cols:
                    s = pd.to_numeric(df[p], errors="coerce")
                    if agg == "median":
                        row[p] = s.median(skipna=True)
                    elif agg == "last":
                        s2 = s.dropna()
                        row[p] = float(s2.iloc[-1]) if len(s2) else np.nan
                    else:
                        row[p] = s.mean(skipna=True)
                rows.append(row)

            mat = pd.DataFrame(rows)
        else:
            stack = []
            for stem in stems:
                exp = self.experiments[stem]
                df = exp.df
                cols = [c for c in exp.numeric_cols if c != "cycle"]
                if cols:
                    stack.append(df[cols].copy())
            if not stack:
                messagebox.showwarning("No data", "No numeric parameters found in the filtered experiments.")
                return
            mat = pd.concat(stack, axis=0, ignore_index=True)

        mat = mat.apply(pd.to_numeric, errors="coerce")
        mat = mat.dropna(axis=0, how="all")

        if target not in mat.columns:
            messagebox.showwarning(
                "Target missing after stacking",
                f"Target '{target}' was not found as a numeric column after building the matrix."
            )
            return

        if mat.shape[0] < 2:
            messagebox.showwarning("Not enough data", "Need at least 2 observations to compute correlations.")
            return

        tgt = mat[target]
        out = {}
        for p in mat.columns:
            if p == target:
                continue
            r = tgt.corr(mat[p], method=method)
            out[p] = r

        s = pd.Series(out).dropna()
        if s.empty:
            messagebox.showwarning("No correlations", "All correlations are NaN (not enough paired data or constant values).")
            return

        self._tcorr_last_series = s
        mode_short = "all cycles" if "all cycles" in mode.lower() else "per experiment"
        self._tcorr_last_title = f"Correlation with {target} ({method}, {mode_short})"
        self._tcorr_render_last()

    def _tcorr_render_last(self):
        if self._tcorr_last_series is None:
            return

        s = self._tcorr_last_series.copy()
        minabs = float(self.tcorr_minabs_var.get())
        topn = int(self.tcorr_topn_var.get())
        sort_abs = bool(self.tcorr_abs_var.get())

        s = s[s.abs() >= minabs]
        if s.empty:
            if hasattr(self, "tcorr_lb"):
                self.tcorr_lb.delete(0, tk.END)
            self.tcorr_ax.clear()
            self._auto["tcorr"]["title"] = "No correlations in range"
            self.tcorr_ax.set_title(self._auto["tcorr"]["title"])
            self.tcorr_canvas.draw_idle()
            self._apply_styles_all()
            return

        if sort_abs:
            s = s.reindex(s.abs().sort_values(ascending=False).index)
        else:
            s = s.sort_values(ascending=False)

        s_plot = s.head(topn)

        if hasattr(self, "tcorr_lb"):
            self.tcorr_lb.delete(0, tk.END)
            for name, r in list(s.items())[:500]:
                self.tcorr_lb.insert(tk.END, f"{name}: r={r:+.3f}  |r|={abs(r):.3f}")

        self.tcorr_ax.clear()
        y = np.arange(len(s_plot))
        self.tcorr_ax.barh(y, s_plot.values)
        self.tcorr_ax.set_yticks(y)
        self.tcorr_ax.set_yticklabels(list(s_plot.index), fontsize=8)
        self.tcorr_ax.invert_yaxis()

        self._auto["tcorr"]["title"] = self._tcorr_last_title + (f", |r| ≥ {minabs:.2f}" if minabs > 0 else "")
        self.tcorr_ax.axvline(0, color="0.4", linewidth=0.8)
        self._auto["tcorr"]["xlabel"] = "Correlation r"
        self._auto["tcorr"]["ylabel"] = ""

        
        self.tcorr_canvas.draw_idle()
        self._apply_styles_all()


# ======================================================================
# SECTION 15: APPENDIX PAGE - page composition helpers (no Tk needed)
# ======================================================================
PAGE_SIZES_CM = {
    "A4 portrait": (21.0, 29.7),
    "A4 landscape": (29.7, 21.0),
    "Letter portrait": (21.59, 27.94),
    "Letter landscape": (27.94, 21.59),
}
HEADING_CM = 1.2  # height reserved for the optional page heading

# tab text in the main app -> attribute name of its matplotlib Figure
SOURCE_FIGS = {
    "Params vs Cycles": "pvc_fig",
    "Experiments vs Cycles": "evc_fig",
    "Scatter": "sc_fig",
    "Correlation Matrix": "corr_fig",
    "Target Correlations": "tcorr_fig",
}


def _wrap_title(text, width_in, fontsize):
    """Wrap a title so it fits the panel width (approximation)."""
    n = max(12, int(width_in * 72.0 / (max(float(fontsize), 1.0) * 0.58)))
    lines = []
    for line in str(text).split("\n"):
        lines.extend(textwrap.wrap(line, n, break_long_words=False) or [""])
    return "\n".join(lines)


def _scale_fig_fonts(fig, scale):
    """Multiply every font size in the figure by `scale`."""
    if abs(scale - 1.0) < 1e-9:
        return
    for ax in fig.axes:
        texts = [ax.title, ax.xaxis.label, ax.yaxis.label,
                 ax.xaxis.get_offset_text(), ax.yaxis.get_offset_text()]
        texts += list(ax.texts)
        leg = ax.get_legend()
        if leg is not None:
            texts += list(leg.get_texts())
            try:
                texts.append(leg.get_title())
            except Exception:
                pass
        for t in texts:
            try:
                t.set_fontsize(float(t.get_fontsize()) * scale)
            except Exception:
                pass
        for name, axis in (("x", ax.xaxis), ("y", ax.yaxis)):
            ticks = axis.get_major_ticks()
            if ticks:
                size = float(ticks[0].label1.get_fontsize())
                ax.tick_params(axis=name, which="both", labelsize=size * scale)


def render_panel_image(fig_bytes, title, w_px, h_px, dpi, font_scale):
    """Un-pickle a stored figure, resize it to the panel size, return RGB array."""
    fig = pickle.loads(fig_bytes)
    fig.set_dpi(dpi)
    fig.set_size_inches(w_px / dpi, h_px / dpi)
    fig.set_facecolor("white")

    _scale_fig_fonts(fig, font_scale)

    if fig.axes:
        main = fig.axes[0]
        txt = (title or "").strip("\n")
        if txt.strip():
            fs = main.title.get_fontsize()
            main.title.set_text(_wrap_title(txt, w_px / dpi, fs))
        else:
            main.title.set_text("")

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    return np.asarray(canvas.buffer_rgba())[:, :, :3].copy()


def _render_heading(text, w_px, h_px, dpi):
    fig = Figure(figsize=(w_px / dpi, h_px / dpi), dpi=dpi)
    fig.set_facecolor("white")
    fig.text(0.5, 0.45, text, ha="center", va="center", fontsize=13, fontweight="bold")
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    return np.asarray(canvas.buffer_rgba())[:, :, :3].copy()


def _paste(page, img, x, y, max_w, max_h):
    h = min(img.shape[0], max_h, page.shape[0] - y)
    w = min(img.shape[1], max_w, page.shape[1] - x)
    if h > 0 and w > 0:
        page[y:y + h, x:x + w] = img[:h, :w]


def compose_page(panels, cfg, dpi):
    """
    panels: list of dicts {"fig": bytes, "title": str}  (at most rows*cols are used)
    cfg: dict(page, rows, cols, margin_cm, gap_cm, font_scale, heading)
    Returns a PIL RGB image of the full page at `dpi`.
    """
    pw_cm, ph_cm = PAGE_SIZES_CM[cfg["page"]]

    def px(cm):
        return int(round(cm / 2.54 * dpi))

    W, H = px(pw_cm), px(ph_cm)
    page = np.full((H, W, 3), 255, dtype=np.uint8)

    m, g = px(cfg["margin_cm"]), px(cfg["gap_cm"])
    rows, cols = cfg["rows"], cfg["cols"]
    header = px(HEADING_CM) if cfg["heading"] else 0

    area_w = W - 2 * m
    area_h = H - 2 * m - header
    cw = (area_w - g * (cols - 1)) // cols
    ch = (area_h - g * (rows - 1)) // rows
    if cw < 40 or ch < 40:
        raise ValueError("The graphs do not fit on the page - reduce the margin / gap "
                         "or the number of rows / columns.")

    if header:
        strip = _render_heading(cfg["heading"], area_w, header, dpi)
        _paste(page, strip, m, m, area_w, header)

    for k, p in enumerate(panels[:rows * cols]):
        r, c = divmod(k, cols)
        x = m + c * (cw + g)
        y = m + header + r * (ch + g)
        img = render_panel_image(p["fig"], p["title"], cw, ch, dpi, cfg["font_scale"])
        _paste(page, img, x, y, cw, ch)

    return Image.fromarray(page)


def save_page_images(pages, path, dpi):
    """Save PIL pages. PDF -> one multi-page file; other formats -> one file per page."""
    base, ext = os.path.splitext(path)
    if not ext:
        ext = ".png"
        path = base + ext
    ext_l = ext.lower()
    if ext_l == ".pdf":
        pages[0].save(path, "PDF", resolution=float(dpi), save_all=True,
                      append_images=list(pages[1:]))
        return [path]
    written = []
    for i, pg in enumerate(pages):
        p = path if len(pages) == 1 else f"{base}_page{i + 1}{ext}"
        kw = {"dpi": (dpi, dpi)}
        if ext_l in (".jpg", ".jpeg"):
            kw["quality"] = 95
        pg.save(p, **kw)
        written.append(p)
    return written


# ======================================================================
# The Tk tab
# ======================================================================
class AppendixPage:
    TAB_TEXT = "Appendix Page"

    def __init__(self, app):
        self.app = app
        self.panels = []          # list of {"fig": bytes, "title": str, "orig": str}
        self.page_idx = 0
        self._loading = False
        self._after_id = None

        self.v_page = tk.StringVar(master=app, value="A4 portrait")
        self.v_rows = tk.IntVar(master=app, value=3)
        self.v_cols = tk.IntVar(master=app, value=2)
        self.v_margin = tk.DoubleVar(master=app, value=2.0)
        self.v_gap = tk.DoubleVar(master=app, value=0.5)
        self.v_fscale = tk.DoubleVar(master=app, value=0.75)
        self.v_heading = tk.StringVar(master=app, value="")
        self.v_dpi = tk.IntVar(master=app, value=300)
        self.v_title = tk.StringVar(master=app, value="")

        self._build_tab()
        self._install_topbar_button()

        for v in (self.v_page, self.v_rows, self.v_cols, self.v_margin, self.v_gap,
                  self.v_fscale, self.v_heading):
            v.trace_add("write", self._schedule_preview)
        self.v_title.trace_add("write", self._on_title_edit)

        app.after(150, self._render_preview)

    # ---------- settings ----------
    def _cfg(self):
        def num(var, default, lo, hi):
            try:
                v = float(var.get())
            except Exception:
                v = default
            return max(lo, min(hi, v))

        page = self.v_page.get()
        if page not in PAGE_SIZES_CM:
            page = "A4 portrait"
        return {
            "page": page,
            "rows": int(num(self.v_rows, 3, 1, 6)),
            "cols": int(num(self.v_cols, 2, 1, 4)),
            "margin_cm": num(self.v_margin, 2.0, 0.0, 6.0),
            "gap_cm": num(self.v_gap, 0.5, 0.0, 3.0),
            "font_scale": num(self.v_fscale, 0.75, 0.2, 2.5),
            "heading": self.v_heading.get().strip(),
            "dpi": int(num(self.v_dpi, 300, 72, 1200)),
        }

    @staticmethod
    def _capacity(cfg):
        return cfg["rows"] * cfg["cols"]

    def _num_pages(self, cfg=None):
        cfg = cfg or self._cfg()
        cap = self._capacity(cfg)
        return max(1, -(-len(self.panels) // cap))

    # ---------- UI ----------
    def _build_tab(self):
        app = self.app
        root = ttk.Frame(app.nb)
        self.tab = root
        left_container, left = app._make_scrollable_left(root)
        right = ttk.Frame(root)
        left_container.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        ttk.Label(
            left, wraplength=300, foreground="#555",
            text="Show a graph in any tab, then press “→ Appendix page” in the top bar "
                 "to send it here."
        ).pack(anchor="w", pady=(0, 6))

        # --- graph list ---
        box = ttk.LabelFrame(left, text="Graphs (page order)")
        box.pack(fill=tk.X, pady=(0, 6))
        row = ttk.Frame(box)
        row.pack(fill=tk.X, padx=6, pady=(4, 2))
        self.lb = tk.Listbox(row, width=38, height=9, exportselection=False)
        self.lb.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb = ttk.Scrollbar(row, orient="vertical", command=self.lb.yview)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.lb.config(yscrollcommand=sb.set)
        self.lb.bind("<<ListboxSelect>>", self._on_select)

        b1 = ttk.Frame(box)
        b1.pack(fill=tk.X, padx=6, pady=(2, 6))
        ttk.Button(b1, text="▲", width=3, command=lambda: self._move(-1)).pack(side=tk.LEFT)
        ttk.Button(b1, text="▼", width=3, command=lambda: self._move(+1)).pack(side=tk.LEFT, padx=(4, 0))
        ttk.Button(b1, text="Remove", command=self._remove).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(b1, text="Clear all", command=self._clear).pack(side=tk.LEFT, padx=(4, 0))

        # --- title of selected graph ---
        tb = ttk.LabelFrame(left, text="Title of selected graph (empty = no title)")
        tb.pack(fill=tk.X, pady=(0, 6))
        ttk.Entry(tb, textvariable=self.v_title).pack(fill=tk.X, padx=6, pady=(4, 2))
        ttk.Button(tb, text="↺ Original title", command=self._reset_title).pack(anchor="w", padx=6, pady=(0, 6))

        # --- page layout ---
        sbx = ttk.LabelFrame(left, text="Page layout")
        sbx.pack(fill=tk.X, pady=(0, 6))

        def add_row(r, label, widget):
            ttk.Label(sbx, text=label).grid(row=r, column=0, sticky="w", padx=6, pady=2)
            widget.grid(row=r, column=1, sticky="w", padx=6, pady=2)

        add_row(0, "Page", ttk.Combobox(sbx, state="readonly", width=16,
                                        values=list(PAGE_SIZES_CM), textvariable=self.v_page))
        add_row(1, "Rows", ttk.Spinbox(sbx, from_=1, to=6, width=7, textvariable=self.v_rows))
        add_row(2, "Columns", ttk.Spinbox(sbx, from_=1, to=4, width=7, textvariable=self.v_cols))
        add_row(3, "Margin (cm)", ttk.Spinbox(sbx, from_=0, to=6, increment=0.25, format="%.2f",
                                              width=7, textvariable=self.v_margin))
        add_row(4, "Gap (cm)", ttk.Spinbox(sbx, from_=0, to=3, increment=0.1, format="%.2f",
                                           width=7, textvariable=self.v_gap))
        add_row(5, "Font scale", ttk.Spinbox(sbx, from_=0.3, to=2.5, increment=0.05, format="%.2f",
                                             width=7, textvariable=self.v_fscale))
        add_row(6, "Page heading", ttk.Entry(sbx, width=22, textvariable=self.v_heading))
        add_row(7, "Save DPI", ttk.Spinbox(sbx, from_=100, to=1200, increment=50,
                                           width=7, textvariable=self.v_dpi))

        ttk.Label(
            left, wraplength=300, foreground="#555",
            text="Font scale multiplies the font sizes of the graphs "
                 "(set in the other tabs) when they are shrunk onto the page."
        ).pack(anchor="w", pady=(0, 6))

        ttk.Button(left, text="Save page(s) as image…", command=self.save_pages).pack(anchor="w", pady=(2, 6))

        # --- preview ---
        nav = ttk.Frame(right)
        nav.pack(side=tk.TOP, fill=tk.X, pady=(0, 4))
        ttk.Button(nav, text="◀", width=3, command=lambda: self._goto(-1)).pack(side=tk.LEFT)
        self.page_lbl = ttk.Label(nav, text="Page 1 / 1", width=14, anchor="center")
        self.page_lbl.pack(side=tk.LEFT, padx=6)
        ttk.Button(nav, text="▶", width=3, command=lambda: self._goto(+1)).pack(side=tk.LEFT)
        ttk.Label(nav, text="(preview - the saved file uses the full DPI)",
                  foreground="#666").pack(side=tk.LEFT, padx=12)

        self.fig = Figure(figsize=(6, 8), dpi=100)
        self.fig.patch.set_facecolor("#dcdcdc")
        self.ax = self.fig.add_axes([0.01, 0.01, 0.98, 0.98])
        self.ax.axis("off")
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        toolbar = NavigationToolbar2Tk(self.canvas, right)
        toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        app.nb.add(root, text=self.TAB_TEXT)

    def _install_topbar_button(self):
        app = self.app
        save_btn = None
        for child in app.winfo_children():
            for w in child.winfo_children():
                try:
                    if isinstance(w, ttk.Button) and str(w.cget("text")).startswith("Save (tight)"):
                        save_btn = w
                        break
                except tk.TclError:
                    pass
            if save_btn is not None:
                break

        orig_save = app.save_current_figure

        def save_dispatch():
            if app.nb.select() == str(self.tab):
                self.save_pages()
            else:
                orig_save()

        if save_btn is not None:
            save_btn.configure(command=save_dispatch)
            ttk.Button(save_btn.master, text="→ Appendix page",
                       command=self.send_current).pack(side=tk.LEFT, padx=(8, 0), after=save_btn)
        else:  # fallback: own thin bar above the notebook
            bar = ttk.Frame(app)
            bar.pack(side=tk.TOP, fill=tk.X, padx=10, before=app.nb)
            ttk.Button(bar, text="→ Appendix page", command=self.send_current).pack(side=tk.LEFT)

    # ---------- list handling ----------
    def _refresh_list(self, select=None):
        cfg = self._cfg()
        cap = self._capacity(cfg)
        self.lb.delete(0, tk.END)
        for i, p in enumerate(self.panels):
            pg, pos = divmod(i, cap)
            t = (p["title"] or "").replace("\n", " ").strip() or "(no title)"
            self.lb.insert(tk.END, f"p{pg + 1}·{pos + 1}   {t}")
        if select is not None and 0 <= select < len(self.panels):
            self.lb.selection_clear(0, tk.END)
            self.lb.selection_set(select)
            self.lb.see(select)
        n = len(self.panels)
        try:
            self.app.nb.tab(self.tab, text=f"{self.TAB_TEXT} ({n})" if n else self.TAB_TEXT)
        except tk.TclError:
            pass

    def _selected(self):
        sel = self.lb.curselection()
        return sel[0] if sel else None

    def _on_select(self, _event=None):
        i = self._selected()
        self._loading = True
        try:
            self.v_title.set(self.panels[i]["title"] if i is not None else "")
        finally:
            self._loading = False
        if i is not None:
            self.page_idx = i // self._capacity(self._cfg())
            self._schedule_preview()

    def _on_title_edit(self, *_):
        if self._loading:
            return
        i = self._selected()
        if i is None:
            return
        self.panels[i]["title"] = self.v_title.get()
        self._refresh_list(select=i)
        self._schedule_preview()

    def _reset_title(self):
        i = self._selected()
        if i is None:
            return
        self.v_title.set(self.panels[i]["orig"])

    def _move(self, delta):
        i = self._selected()
        if i is None:
            return
        j = i + delta
        if not (0 <= j < len(self.panels)):
            return
        self.panels[i], self.panels[j] = self.panels[j], self.panels[i]
        self.page_idx = j // self._capacity(self._cfg())
        self._refresh_list(select=j)
        self._schedule_preview()

    def _remove(self):
        i = self._selected()
        if i is None:
            return
        del self.panels[i]
        self._loading = True
        self.v_title.set("")
        self._loading = False
        self._refresh_list(select=min(i, len(self.panels) - 1) if self.panels else None)
        self._on_select()
        self._schedule_preview()

    def _clear(self):
        if not self.panels:
            return
        if not messagebox.askyesno("Clear all", "Remove all graphs from the appendix page?"):
            return
        self.panels.clear()
        self.page_idx = 0
        self._loading = True
        self.v_title.set("")
        self._loading = False
        self._refresh_list()
        self._schedule_preview()

    # ---------- sending ----------
    def send_current(self):
        app = self.app
        if app.nb.select() == str(self.tab):
            messagebox.showinfo("Appendix page",
                                "Switch to a plotting tab and send the graph from there.")
            return
        tab_text = app.nb.tab(app.nb.select(), "text")
        fig = getattr(app, SOURCE_FIGS.get(tab_text, ""), None)
        ax = fig.axes[0] if (fig is not None and fig.axes) else None
        if ax is None or not (ax.lines or ax.collections or ax.images or ax.patches):
            messagebox.showwarning("Nothing to send", "Plot a graph in this tab first.")
            return
        try:
            data = pickle.dumps(fig)
        except Exception as e:
            messagebox.showerror("Could not copy the graph", f"{type(e).__name__}: {e}")
            return

        title = ax.get_title() or ""
        self.panels.append({"fig": data, "title": title, "orig": title})
        i = len(self.panels) - 1
        self.page_idx = i // self._capacity(self._cfg())
        self._refresh_list(select=i)
        self._on_select()
        self._schedule_preview()
        app.status_var.set(f"Sent to Appendix page ({len(self.panels)} graph(s) so far).")

    # ---------- preview ----------
    def _goto(self, delta):
        n = self._num_pages()
        self.page_idx = max(0, min(self.page_idx + delta, n - 1))
        self._render_preview()

    def _schedule_preview(self, *_):
        if self._after_id is not None:
            try:
                self.app.after_cancel(self._after_id)
            except Exception:
                pass
        self._after_id = self.app.after(300, self._render_preview)

    def _render_preview(self):
        self._after_id = None
        ax = self.ax
        ax.clear()
        ax.axis("off")
        cfg = self._cfg()
        cap = self._capacity(cfg)
        npages = self._num_pages(cfg)
        self.page_idx = max(0, min(self.page_idx, npages - 1))
        self.page_lbl.config(text=f"Page {self.page_idx + 1} / {npages}")

        if not self.panels:
            ax.text(0.5, 0.5, "No graphs yet.\nPlot something in another tab,\n"
                              "then press “→ Appendix page”.",
                    ha="center", va="center", fontsize=12, color="#555")
            self.canvas.draw_idle()
            return

        chunk = self.panels[self.page_idx * cap:(self.page_idx + 1) * cap]
        try:
            img = compose_page(chunk, cfg, dpi=100)
        except Exception as e:
            ax.text(0.5, 0.5, f"{e}", ha="center", va="center", fontsize=11,
                    color="#a00", wrap=True)
            self.canvas.draw_idle()
            return
        ax.imshow(np.asarray(img))
        ax.axis("off")
        self.canvas.draw_idle()
        self._refresh_list(select=self._selected())

    # ---------- saving ----------
    def save_pages(self):
        if not self.panels:
            messagebox.showinfo("Nothing to save", "Send at least one graph to the appendix page first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".png", initialfile="appendix_page.png",
            filetypes=[("PNG", "*.png"), ("PDF (all pages in one file)", "*.pdf"),
                       ("JPEG", "*.jpg"), ("TIFF", "*.tif"), ("All files", "*.*")]
        )
        if not path:
            return
        cfg = self._cfg()
        cap = self._capacity(cfg)
        dpi = cfg["dpi"]
        self.app.config(cursor="watch")
        self.app.update_idletasks()
        try:
            pages = [compose_page(self.panels[i:i + cap], cfg, dpi)
                     for i in range(0, len(self.panels), cap)]
            written = save_page_images(pages, path, dpi)
        except Exception as e:
            messagebox.showerror("Save failed", f"{type(e).__name__}: {e}")
            return
        finally:
            self.app.config(cursor="")
        self.app.status_var.set("Saved: " + ", ".join(os.path.basename(p) for p in written))


def install_appendix_tab(app):
    """Attach the Appendix Page tab + '→ Appendix page' button to a CVExplorer window."""
    inst = AppendixPage(app)
    app.appendix = inst
    return inst


if __name__ == "__main__":
    app = CVExplorer()
    install_appendix_tab(app)  # adds the "Appendix Page" tab + "→ Appendix page" button
    app.mainloop()