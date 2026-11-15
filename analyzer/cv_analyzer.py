"""
CV Analyzer - desktop tool for analyzing cyclic voltammetry (CV) data.

Reads EC-Lab ASCII (.mpt) files and the CSV / TXT output of MPR Converter
(https://mpr-converter.vercel.app), browses files and cycles in a tree,
converts potentials to RHE, finds anodic/cathodic redox peaks, integrates
charge, optionally subtracts the OER background, computes derivatives with
smoothing, and exports per-cycle metrics to CSV. Includes an overlay window
for comparing the same cycle across many files, a printable appendix-page
builder, and a guided first-run tutorial (English / Hebrew) on synthetic data.

Author: Aviv Bryger
Run:    python analyzer/cv_analyzer.py
"""

import tkinter as tk
from tkinter import filedialog, ttk, messagebox, simpledialog
from tkinter.colorchooser import askcolor
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.legend import Legend
from matplotlib.backends.backend_agg import FigureCanvasAgg
from PIL import Image
import numpy as np
import os
import re
import pickle
import textwrap
import json
import time

# Optional SciPy for Savitzky–Golay
try:
    from scipy.signal import savgol_filter
    HAS_SAVGOL = True
except Exception:
    HAS_SAVGOL = False

# -------- Default folder path --------
# ---- Export dialog option lists ----
ANIONS = ['Iodide', 'Chloride', 'Bromide', 'Carbonate', 'Acetate',
          'Cyanate', 'Hydroxide', 'Nitrate', 'Nitrite', 'Sulphate']
STATIONS = ['A', 'B', 'C', 'D']
MEAS_TYPES = ['UAL-ACT-night', 'UOS-night']
LABORANTS = ['User']   # starting list; edited in the export window ("Edit list...")


def get_laborants():
    """The saved laborant list (falls back to LABORANTS until a list is saved)."""
    names = overlay_load_settings().get('laborants')
    if isinstance(names, list):
        return [str(n) for n in names if str(n).strip()]
    return list(LABORANTS)


def save_laborants(names):
    st = overlay_load_settings()
    st['laborants'] = [str(n).strip() for n in names if str(n).strip()]
    overlay_save_settings(st)


def remember_laborant(name):
    """Add a newly typed name to the list and make it the default next time."""
    name = str(name).strip()
    if not name:
        return
    st = overlay_load_settings()
    names = get_laborants()
    if name.lower() not in [n.lower() for n in names]:
        names.append(name)
    st['laborants'] = names
    st['last_laborant'] = name
    overlay_save_settings(st)


def default_laborant():
    names = get_laborants()
    last = overlay_load_settings().get('last_laborant')
    if last in names:
        return last
    return names[0] if names else ''


# Chemical shorthand that may appear in a .mpt file name instead of the full
# anion word. Matched against whole tokens only - a bare "cl" or "i" as a
# substring would fire on far too many words.
ANION_ALIASES = {
    'no3': 'Nitrate', 'no2': 'Nitrite', 'co3': 'Carbonate', 'so4': 'Sulphate',
    'oh': 'Hydroxide', 'cl': 'Chloride', 'br': 'Bromide', 'i': 'Iodide',
    'ac': 'Acetate', 'oac': 'Acetate', 'ocn': 'Cyanate', 'cno': 'Cyanate',
}

# Short measurement codes as they appear in file names.
MEAS_PREFIX = {'ual': 'UAL-ACT-night', 'uos': 'UOS-night'}


def guess_name_parts(mpt_name):
    """Best-effort read of the export name parts out of a .mpt file name.

    Example:  UAL_A_35710_nitrite_9_6.mpt
          ->  anion=Nitrite, station=A, electrode=35710,
              meas=UAL-ACT-night, laborant=<default>

    Every value is only a suggestion: it is pre-filled into the export
    dialog for the user to confirm or correct, never saved automatically.
    Fields that cannot be determined come back as '' (empty).
    """
    out = {'anion': '', 'station': '', 'electrode': '', 'meas': '', 'laborant': ''}
    if not mpt_name:
        return out

    stem = os.path.splitext(str(mpt_name))[0]
    low = stem.lower()
    tokens = [t for t in re.split(r'[^0-9A-Za-z]+', stem) if t]
    lowtok = [t.lower() for t in tokens]

    # --- Anion: full word anywhere, else a shorthand token ---
    for a in ANIONS:
        if a.lower() in low:
            out['anion'] = a
            break
    if not out['anion']:
        for t in lowtok:
            if t in ANION_ALIASES:
                out['anion'] = ANION_ALIASES[t]
                break

    # --- Station: a standalone single letter A-D ---
    for t in tokens:
        if len(t) == 1 and t.upper() in STATIONS:
            out['station'] = t.upper()
            break

    # --- Measurement type: full name first, then the short code ---
    for m in MEAS_TYPES:
        if m.lower() in low:
            out['meas'] = m
            break
    if not out['meas']:
        for t in lowtok:
            if t in MEAS_PREFIX:
                out['meas'] = MEAS_PREFIX[t]
                break
    if not out['meas']:
        # e.g. "UAL-ACT-day..." where the token split keeps 'ual' attached
        for code, full in MEAS_PREFIX.items():
            if low.startswith(code):
                out['meas'] = full
                break

    # --- Laborant: named in the file, else the default ---
    for l in get_laborants():
        if l.lower() in low:
            out['laborant'] = l
            break
    if not out['laborant']:
        out['laborant'] = default_laborant()

    # --- Electrode: the longest numeric token of 3+ digits.
    # The 3-digit floor keeps short date fragments such as the "9_6" in
    # "..._nitrite_9_6" from being mistaken for an electrode number.
    nums = [t for t in tokens if t.isdigit() and len(t) >= 3]
    if nums:
        out['electrode'] = max(nums, key=len)

    return out

# Default folders. Both can be changed in the program (File > Set default
# folders...) and the choice is remembered in ~/.cv_overlay_settings.json.
# Folder the file tree ("Open Folder") starts in and loads at start-up
# ('' = none until the user chooses one)
DEFAULT_FOLDER = ''

# Folder the "Export Analysis CSV" dialog writes to by default
EXPORT_FOLDER = os.path.join(os.path.expanduser('~'), 'CV_Analyzer_exports')

# =================== Robust EC-Lab ASCII reader (CV + OCV) ===================
# Data files the library lists: EC-Lab ASCII exports (.mpt), and the CSV / TXT
# files MPR Converter writes from binary .mpr files. Converter output keeps the
# original EC-Lab column names (Ewe/V, <I>/mA, cycle number, time/s), so the
# same reader handles both.
DATA_EXTS = ('.mpt', '.txt', '.csv')
DATA_FILETYPES = [("CV data (.mpt, MPR Converter .csv/.txt)", "*.mpt *.csv *.txt"),
                  ("EC-Lab ASCII", "*.mpt"),
                  ("MPR Converter output", "*.csv *.txt")]


def _header_delimiter(line):
    """Column separator of a header line: tab (EC-Lab, converter .txt),
    ';' or ',' (converter .csv)."""
    if '\t' in line:
        return '\t'
    if ';' in line:
        return ';'
    if ',' in line:
        return ','
    return '\t'


def is_cv_data_file(path):
    """True for a .mpt, or for a .csv / .txt whose header has a potential and
    a current column (MPR Converter output). Other CSV / TXT files in the
    folder, including this program's own analysis exports, are left out.
    Reads only the first lines, so listing a folder stays fast."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.mpt':
        return True
    if ext not in DATA_EXTS:
        return False
    try:
        with open(path, 'r', encoding='latin1') as f:
            for _ in range(5):
                s = f.readline()
                if not s:
                    break
                if 'Ewe' in s and ('<I>' in s or 'I/mA' in s or 'I/A' in s):
                    return True
    except OSError:
        pass
    return False


def read_mpt(path):
    """
    Returns:
        header, data, E (np.ndarray), I (np.ndarray), t (np.ndarray), cyc (np.ndarray), meta(dict)
        meta['ocv'] == True if file detected as OCV (no current column)
    """
    with open(path, 'r', encoding='latin1') as f:
        lines = f.readlines()

    header_index = None
    mode_hint = None  # "CV" or "OCV"
    for i, line in enumerate(lines):
        s = line.strip()
        if ('Ewe/V' in s or 'Ewe' in s) and ('<I>' in s or 'I/mA' in s or '</mA>' in s or '<I>/A' in s):
            header_index = i
            mode_hint = 'CV'
            break
        if s.lower().startswith('record') and '<ewe>' in s.lower():
            for j in range(i + 1, min(i + 8, len(lines))):
                if ('time' in lines[j].lower()) and ('ewe' in lines[j].lower()):
                    header_index = j
                    mode_hint = 'OCV'
                    break
        if header_index is not None:
            break

    if header_index is None:
        for i, line in enumerate(lines):
            s = line.lower()
            if ('time' in s) and ('ewe' in s):
                header_index = i
                mode_hint = 'OCV'
                break

    if header_index is None:
        raise ValueError("Could not find data header (time/Ewe).")

    delim = _header_delimiter(lines[header_index])
    if delim == '\t':
        # EC-Lab .mpt / MPR Converter .txt. A comma here can only be a
        # decimal comma (European locale exports).
        raw_header = lines[header_index].strip().replace(',', '.')
        header = [h.strip() for h in raw_header.split('\t') if h.strip() != '']
    else:
        # MPR Converter .csv (',') or a ';'-separated export. Cells keep
        # their position, so an empty cell becomes NaN instead of shifting
        # the columns after it.
        header = [h.strip().strip('"') for h in lines[header_index].strip().split(delim)]

    rows = []
    for line in lines[header_index + 1:]:
        if delim == '\t':
            parts = [p for p in line.strip().split('\t') if p != '']
        else:
            parts = [p.strip().strip('"') for p in line.strip().split(delim)]
            if parts == ['']:
                parts = []
        if len(parts) < 1:
            continue
        try:
            if delim == ',':
                row = [float(x) if x != '' else np.nan for x in parts]
            else:
                row = [float(x.replace(',', '.')) if x != '' else np.nan for x in parts]
            rows.append(row)
        except Exception:
            if line.strip().lower().startswith('mode'):
                break
            continue

    if not rows:
        data = np.zeros((0, len(header)), dtype=float)
    else:
        max_cols = max(len(r) for r in rows)
        if max_cols < len(header):
            header = header[:max_cols]
        padded = [r + [np.nan] * (len(header) - len(r)) if len(r) < len(header) else r[:len(header)] for r in rows]
        data = np.array(padded, dtype=float)

    def find(*names):
        for name in names:
            for i, h in enumerate(header):
                if name.lower() in h.lower():
                    return i
        return None

    idx_E = find('Ewe/V', 'Ewe')
    idx_I = find('<I>/mA', '<I>/A', 'I/mA', '<I>')
    idx_t = find('time/s', 'time')
    idx_cyc = find('cycle number', 'cycle', 'Cycle', 'Ns')

    if idx_E is None:
        raise ValueError("Could not find Ewe column in data.")

    E = data[:, idx_E] if data.shape[0] else np.array([], dtype=float)
    if idx_I is None:
        I = np.zeros_like(E)
    else:
        I = data[:, idx_I]

    if idx_t is None:
        t = np.arange(len(E), dtype=float)
    else:
        t = data[:, idx_t]

    if idx_cyc is None or len(data) == 0:
        cyc = np.zeros(len(E), dtype=int)
    else:
        cyc = data[:, idx_cyc].astype(int, copy=False)

    meta = {'ocv': idx_I is None, 'mode_hint': mode_hint, 'header_line': header_index}
    return header, data, E, I, t, cyc, meta

# =================== Reference conversion ===================
def convert_reference(E, ref_type, pH):
    if ref_type == 'RHE':
        return E
    REF = {
        'Ag/AgCl 3M KCl': 0.205,
        'Ag/AgCl sat. KCl': 0.197,
        'Hg/HgO 1M KOH': 0.098,
        'Hg/HgO 0.1M KOH': 0.165,
        'SCE': 0.241
    }
    return E + REF.get(ref_type, 0) + 0.059 * pH

# =================== Q integration (length-safe) ===================
def integrate_Q(I, t):
    n = len(I)
    if n == 0:
        return np.array([], dtype=float)
    if n == 1:
        return np.zeros(1, dtype=float)
    incremental = (I[1:] + I[:-1]) * 0.5 * np.diff(t)
    out = np.zeros(n, dtype=float)
    out[1:] = np.cumsum(incremental)
    return out

# =================== OER subtraction (Mode A) ===================
def oer_subtract_modeA(E, I):
    debug = {'valid': False, 'fw_pos': 0, 'bw_pos': 0, 'applied': 0}
    n = len(E)
    if n < 3:
        return I.copy(), np.zeros_like(I, dtype=bool), debug

    idx_sw = int(np.argmax(E))
    fw_mask = np.arange(n) <= idx_sw
    bw_mask = ~fw_mask

    Ef, If = E[fw_mask], I[fw_mask]
    Eb, Ib = E[bw_mask], I[bw_mask]

    pos_f = If > 0
    pos_b = Ib > 0
    debug['fw_pos'] = int(np.sum(pos_f))
    debug['bw_pos'] = int(np.sum(pos_b))
    if np.sum(pos_f) < 3 or np.sum(pos_b) < 3:
        return I.copy(), np.zeros_like(I, dtype=bool), debug

    Ef_pos, If_pos = Ef[pos_f], If[pos_f]
    Eb_pos, Ib_pos = Eb[pos_b], Ib[pos_b]

    of = np.argsort(Ef_pos)
    ob = np.argsort(Eb_pos)
    Ef_pos = Ef_pos[of]; If_pos = If_pos[of]
    Eb_pos = Eb_pos[ob]; Ib_pos = Ib_pos[ob]

    E_low = max(Ef_pos.min(), Eb_pos.min())
    E_high = min(Ef_pos.max(), Eb_pos.max())
    if not (E_high > E_low):
        return I.copy(), np.zeros_like(I, dtype=bool), debug

    in_fw_win = (Ef_pos >= E_low) & (Ef_pos <= E_high)
    if np.sum(in_fw_win) < 3:
        return I.copy(), np.zeros_like(I, dtype=bool), debug

    Ib_interp = np.interp(Ef_pos[in_fw_win], Eb_pos, Ib_pos)

    I_blend = I.copy()
    overlap_mask = np.zeros_like(I, dtype=bool)

    fw_idx_all = np.where(fw_mask)[0]
    fw_pos_idx = fw_idx_all[pos_f][of]
    apply_idx = fw_pos_idx[in_fw_win]

    I_blend[apply_idx] = If_pos[in_fw_win] - Ib_interp
    overlap_mask[apply_idx] = True

    debug['valid'] = True
    debug['applied'] = int(np.sum(overlap_mask))
    return I_blend, overlap_mask, debug

def _split_forward_backward_by_turnaround(self, E, I):
    """
    Split a CV trace into forward and backward segments by the turnaround point.
    Uses the index of maximum E as the junction (typical CV).
    Returns: (E_for, I_for, E_back, I_back, idx_turn)
    """
    E = np.asarray(E, dtype=float)
    I = np.asarray(I, dtype=float)
    if len(E) < 5:
        return E, I, np.array([]), np.array([]), None

    idx_turn = int(np.nanargmax(E))  # turnaround at max potential
    # include the turnaround point in both (helps interpolation at the junction)
    E_for  = E[:idx_turn+1]; I_for  = I[:idx_turn+1]
    E_back = E[idx_turn:];  I_back = I[idx_turn:]
    return E_for, I_for, E_back, I_back, idx_turn


def _interp_I_at_E(self, E_seg, I_seg, E_targets):
    """
    Interpolate I(E) for a segment. Handles non-monotonic segments by sorting by E.
    Returns dict: {Etarget: I_interp_or_nan}
    """
    E_seg = np.asarray(E_seg, dtype=float)
    I_seg = np.asarray(I_seg, dtype=float)

    out = {}
    if len(E_seg) < 3:
        for Et in E_targets:
            out[Et] = np.nan
        return out

    # sort by E for safe interpolation
    order = np.argsort(E_seg)
    Es = E_seg[order]
    Is = I_seg[order]

    # remove duplicate E (keep last)
    _, uniq_idx = np.unique(Es, return_index=True)
    Es = Es[uniq_idx]
    Is = Is[uniq_idx]

    Emin, Emax = float(np.nanmin(Es)), float(np.nanmax(Es))
    for Et in E_targets:
        Etf = float(Et)
        if Etf < Emin or Etf > Emax:
            out[Etf] = np.nan
        else:
            out[Etf] = float(np.interp(Etf, Es, Is))
    return out


def _fixed_currents_summary(self, E_use, I_use):
    """
    Compute fixed-voltage currents for forward/back segments + junction point.
    Returns dict with keys like:
      for_I@1.00, back_I@1.00, for_I@1.10, back_I@1.10, ...
      junction_E, junction_I
    """
    targets = [1.00, 1.10, 1.20, 1.375, 1.40, 1.50, 1.55, 1.60]

    E_for, I_for, E_back, I_back, idx_turn = self._split_forward_backward_by_turnaround(E_use, I_use)

    d_for  = self._interp_I_at_E(E_for,  I_for,  targets)
    d_back = self._interp_I_at_E(E_back, I_back, targets)

    out = {}
    for Et in targets:
        out[f"for_I@{Et:.3f}"]  = d_for.get(float(Et), np.nan)
        out[f"back_I@{Et:.3f}"] = d_back.get(float(Et), np.nan)

    # junction around turnaround (your “around 1.3–1.37” wording sounds like this)
    if idx_turn is not None and 0 <= idx_turn < len(E_use):
        out["junction_E"] = float(E_use[idx_turn])
        out["junction_I"] = float(I_use[idx_turn])
    else:
        out["junction_E"] = np.nan
        out["junction_I"] = np.nan

    return out

# =================== Cleaning helpers ===================
def clean_xy(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = min(len(x), len(y))
    x0 = x[:m]
    y0 = y[:m]
    ok = np.isfinite(x0) & np.isfinite(y0)
    removed = int(np.sum(~ok))
    return x0[ok], y0[ok], ok, removed


# =================== Scan-rate estimator ===================
def estimate_scan_rate(E, t):
    if len(E) < 3 or len(t) < 3:
        return np.nan
    dE = np.diff(E)
    dt = np.diff(t)
    mask = dt > 0
    if not np.any(mask):
        return np.nan
    sr = np.nanmedian(np.abs(dE[mask] / dt[mask]))  # V/s
    return sr

# =================== Smoothing helpers ===================
def moving_average(y, w):
    y = np.asarray(y, dtype=float)
    w = int(max(1, w))
    if w == 1 or len(y) < w:
        return y.copy()
    kernel = np.ones(w, dtype=float) / float(w)
    pad = w // 2
    ypad = np.pad(y, (pad, pad), mode='edge')
    ys = np.convolve(ypad, kernel, mode='valid')
    return ys

def safe_gradient(y, x):
    """
    Robust dy/dx with same-length output.
    Handles length mismatch + duplicate x steps.
    """
    import numpy as np

    y = np.asarray(y, dtype=float).ravel()
    x = np.asarray(x, dtype=float).ravel()

    # Align lengths defensively
    n = min(len(y), len(x))
    y = y[:n]
    x = x[:n]

    if n < 3:
        # not enough points for a stable gradient
        return np.full(n, np.nan)

    # Fix non-finite values (optional but helps avoid rare crashes)
    bad = ~(np.isfinite(y) & np.isfinite(x))
    if np.any(bad):
        good = ~bad
        if good.sum() < 3:
            return np.full(n, np.nan)
        idx = np.arange(n)
        y = np.interp(idx, idx[good], y[good])
        x = np.interp(idx, idx[good], x[good])

    # Avoid zero spacing (duplicate potentials) -> jitter very slightly
    dx = np.diff(x)
    if np.any(dx == 0):
        eps = np.finfo(float).eps * max(1.0, np.nanmax(np.abs(x)))
        x = x + eps * np.arange(n)

    # Proper dy/dx with correct shape
    return np.gradient(y, x, edge_order=1)

# ===================== Target currents @ fixed potentials =====================

# Potentials (V vs RHE) at which forward / backward / (fwd-bwd) currents are
# sampled. Defined once here - both the Analysis table columns and the
# computation read this list, so they can never drift apart.
TARGET_E = [1.0, 1.1, 1.2, 1.375, 1.4, 1.45, 1.5, 1.55, 1.6]
JUNCTION_RANGE = (1.30, 1.37)


def _ekey(E):
    # 1.1 -> "1p1", 1.375 -> "1p375"
    s = f"{float(E):g}"
    return s.replace('.', 'p').replace('-', 'm')

def _safe_interp_1d(x, y, xq):
    """
    Interpolate y(x) at xq safely.
    Returns np.nan if xq is out of range or if data is too short.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2:
        return np.nan

    # remove non-finite
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if x.size < 2:
        return np.nan

    # sort by x, and remove duplicates (keep last)
    order = np.argsort(x)
    x = x[order]; y = y[order]
    # unique x
    ux, idx = np.unique(x, return_index=True)
    # np.unique keeps first; we want last -> reverse trick
    ux2, idx2 = np.unique(x[::-1], return_index=True)
    idx_last = (x.size - 1) - idx2
    idx_last = np.sort(idx_last)
    x = x[idx_last]; y = y[idx_last]
    if x.size < 2:
        return np.nan

    if xq < x[0] or xq > x[-1]:
        return np.nan

    return float(np.interp(xq, x, y))


def split_forward_backward(E_seq, I_seq):
    """
    Split a CV cycle into forward/backward sweeps based on the max E point.
    Assumes the sequence is time-ordered.
    """
    E_seq = np.asarray(E_seq, dtype=float)
    I_seq = np.asarray(I_seq, dtype=float)
    if E_seq.size < 4:
        return None, None, None, None

    # peak E index (turning point)
    imax = int(np.nanargmax(E_seq))
    if imax < 1 or imax > E_seq.size - 2:
        return None, None, None, None

    Ef = E_seq[:imax+1]
    If = I_seq[:imax+1]
    Eb = E_seq[imax:]
    Ib = I_seq[imax:]
    return Ef, If, Eb, Ib


def sample_targets_forward_backward(E_seq, I_seq, targets):
    """
    For each target potential, return forward I, backward I, and (forward - backward).
    """
    Ef, If, Eb, Ib = split_forward_backward(E_seq, I_seq)
    out = {}
    if Ef is None:
        # fill NaNs
        for Et in targets:
            k = _ekey(Et)
            out[f"I_fwd_{k}"] = np.nan
            out[f"I_bwd_{k}"] = np.nan
            out[f"I_for_back_{k}"] = np.nan
        return out

    # For interpolation we want x ascending.
    # Forward should be increasing; backward decreasing, so sort each sweep by E for interp.
    for Et in targets:
        k = _ekey(Et)
        I_f = _safe_interp_1d(Ef, If, Et)
        I_b = _safe_interp_1d(Eb, Ib, Et)  # safe_interp sorts internally
        out[f"I_fwd_{k}"] = I_f
        out[f"I_bwd_{k}"] = I_b
        out[f"I_for_back_{k}"] = (I_f - I_b) if (np.isfinite(I_f) and np.isfinite(I_b)) else np.nan

    return out


def find_junction_in_range(E_seq, I_seq, junction_range=(1.30, 1.37), grid_n=250):
    """
    Find approximate forward/backward junction (where Ifwd ~= Ibwd) within a given E range.
    Returns dict: junction_E, junction_I_fwd, junction_I_bwd, junction_for_back
    """
    lo, hi = float(junction_range[0]), float(junction_range[1])
    Ef, If, Eb, Ib = split_forward_backward(E_seq, I_seq)
    if Ef is None:
        return {"junction_E": np.nan, "junction_I_fwd": np.nan, "junction_I_bwd": np.nan, "junction_for_back": np.nan}

    grid = np.linspace(lo, hi, int(grid_n))
    diffs = []
    Ifs = []
    Ibs = []

    for Eg in grid:
        v1 = _safe_interp_1d(Ef, If, Eg)
        v2 = _safe_interp_1d(Eb, Ib, Eg)
        Ifs.append(v1)
        Ibs.append(v2)
        diffs.append(v1 - v2 if (np.isfinite(v1) and np.isfinite(v2)) else np.nan)

    diffs = np.asarray(diffs, dtype=float)
    Ifs = np.asarray(Ifs, dtype=float)
    Ibs = np.asarray(Ibs, dtype=float)
    m = np.isfinite(diffs)
    if not np.any(m):
        return {"junction_E": np.nan, "junction_I_fwd": np.nan, "junction_I_bwd": np.nan, "junction_for_back": np.nan}

    idx = int(np.nanargmin(np.abs(diffs[m])))
    # map back to grid indices
    grid_idx = np.where(m)[0][idx]
    E_j = float(grid[grid_idx])
    I_f = float(Ifs[grid_idx])
    I_b = float(Ibs[grid_idx])
    return {"junction_E": E_j, "junction_I_fwd": I_f, "junction_I_bwd": I_b, "junction_for_back": (I_f - I_b)}


def compute_fixedE_metrics(E_seq, I_seq, targets, junction_range):
    """
    Convenience wrapper: targets + junction.
    """
    out = sample_targets_forward_backward(E_seq, I_seq, targets)
    out.update(find_junction_in_range(E_seq, I_seq, junction_range=junction_range))
    return out

# =================== Laborant list editor ===================
class LaborantListDialog(tk.Toplevel):
    """Add / remove the names offered in the Laborant field. Saved on Close."""

    def __init__(self, parent):
        super().__init__(parent)
        self.title("Laborants")
        self.resizable(False, False)
        tk.Label(self, text="Names offered in the Laborant field:",
                 font=('TkDefaultFont', 9, 'bold')).pack(anchor='w', padx=10, pady=(10, 2))
        self.lb = tk.Listbox(self, height=8, width=30, exportselection=False)
        self.lb.pack(padx=10)
        for n in get_laborants():
            self.lb.insert('end', n)

        add = tk.Frame(self)
        add.pack(fill='x', padx=10, pady=(6, 0))
        self.new_var = tk.StringVar()
        ent = tk.Entry(add, textvariable=self.new_var, width=20)
        ent.pack(side='left')
        ent.bind('<Return>', lambda e: self._add())
        tk.Button(add, text="Add", command=self._add).pack(side='left', padx=4)

        btns = tk.Frame(self)
        btns.pack(fill='x', padx=10, pady=10)
        tk.Button(btns, text="Remove selected", command=self._remove).pack(side='left')
        tk.Button(btns, text="Close", width=10, command=self._close).pack(side='right')
        tk.Label(self, text="Removing a name does not change files already saved.",
                 fg='#777', font=('TkDefaultFont', 8)).pack(anchor='w', padx=10, pady=(0, 8))

        self.protocol("WM_DELETE_WINDOW", self._close)
        self.transient(parent)
        self.grab_set()
        ent.focus_set()
        parent.wait_window(self)

    def _add(self):
        n = self.new_var.get().strip()
        if n and n.lower() not in [x.lower() for x in self.lb.get(0, 'end')]:
            self.lb.insert('end', n)
        self.new_var.set('')

    def _remove(self):
        for i in reversed(self.lb.curselection()):
            self.lb.delete(i)

    def _close(self):
        save_laborants(list(self.lb.get(0, 'end')))
        self.destroy()


# =================== Export-name dialog ===================
class ExportNameDialog(tk.Toplevel):
    """Modal dialog that collects the parts of the exported CSV file name.

    Result (self.result) is the assembled base name, or None if cancelled.
    Name format: Anion_Station_Electrode_MeasurementType_Laborant
    """

    def __init__(self, parent, mpt_names, mode='export', folder=None,
                 progress=None, report=None):
        """mode='export' -> name an exported CSV (result = base name).
        mode='rename'   -> name a data file (.mpt / .csv / .txt).
        mode='batch'    -> name + approve one file of a batch analysis; the
                           analysis report is shown alongside the fields.
        `action` tells the caller whether to proceed, skip, or stop the run."""
        super().__init__(parent)
        self.mode = mode
        self.action = 'cancel'   # 'ok' | 'skip' | 'cancel'
        # Extension kept on rename: .mpt, or .csv / .txt from MPR Converter
        self.src_ext = (os.path.splitext(mpt_names[0])[1] if mpt_names else '') or '.mpt'
        if mode == 'rename':
            title = "Rename data file"
        elif mode == 'batch':
            title = "Batch analysis - review and save"
        else:
            title = "Export Analysis CSV"
        if progress:
            title += f"   ({progress})"
        self.title(title)
        self.resizable(mode == 'batch', mode == 'batch')
        self.result = None
        self.folder = folder if folder is not None else EXPORT_FOLDER

        pad = dict(padx=10, pady=4)

        # ---- header: the .mpt file(s) concerned ----
        head = tk.Frame(self)
        head.pack(fill='x', padx=10, pady=(10, 2))
        tk.Label(head, text=("Current file name:" if mode == 'rename' else "Analyzed file:"),
                 font=('TkDefaultFont', 9, 'bold')).pack(anchor='w')
        if not mpt_names:
            shown = "(unknown)"
        elif len(mpt_names) == 1:
            shown = mpt_names[0]
        else:
            shown = mpt_names[0] + f"   (+{len(mpt_names) - 1} more)"
        tk.Label(head, text=shown, fg='#0a4d8c', wraplength=460, justify='left').pack(anchor='w')

        # Read what we can out of the file name. These are suggestions only -
        # every field stays editable and nothing is saved automatically.
        g = guess_name_parts(mpt_names[0] if mpt_names else '')
        guessed = [k for k in ('anion', 'station', 'electrode', 'meas') if g.get(k)]
        if guessed:
            note = "Auto-filled from the file name: " + ", ".join(guessed) + "  -  please verify"
        else:
            note = "Could not read any field from the file name - please fill in manually"
        tk.Label(head, text=note, fg='#777', wraplength=380,
                 justify='left', font=('TkDefaultFont', 8)).pack(anchor='w', pady=(4, 0))

        # ---- destination folder ----
        tk.Label(head, text=("In folder:" if mode == 'rename' else "Save to folder:"),
                 font=('TkDefaultFont', 9, 'bold')).pack(anchor='w', pady=(8, 0))
        folder_ok = os.path.isdir(self.folder)
        tk.Label(head, text=self.folder,
                 fg=('#0a4d8c' if folder_ok else '#b00020'),
                 wraplength=460, justify='left',
                 font=('TkDefaultFont', 8)).pack(anchor='w')
        if not folder_ok and mode != 'rename':
            tk.Label(head, text="Folder not reachable - it will be created on save.",
                     fg='#b00020', font=('TkDefaultFont', 8)).pack(anchor='w')

        ttk.Separator(self, orient='horizontal').pack(fill='x', padx=10, pady=(6, 2))

        body = tk.Frame(self)
        body.pack(fill='both', expand=True, padx=10, pady=4)

        # ---- Anion: free text + autocomplete ----
        tk.Label(body, text="Anion:").grid(row=0, column=0, sticky='w', pady=4)
        self.anion_var = tk.StringVar(value=g['anion'])
        self.anion_cb = ttk.Combobox(body, textvariable=self.anion_var,
                                     values=ANIONS, width=28)
        self.anion_cb.grid(row=0, column=1, sticky='w', pady=4)
        self.anion_cb.bind('<KeyRelease>', self._anion_autocomplete)
        self.anion_cb.focus_set()

        # ---- Station: exactly one of A/B/C/D ----
        tk.Label(body, text="Station:").grid(row=1, column=0, sticky='w', pady=4)
        self.station_var = tk.StringVar(value=g['station'])
        st = tk.Frame(body)
        st.grid(row=1, column=1, sticky='w', pady=4)
        for v in STATIONS:
            tk.Radiobutton(st, text=v, value=v, variable=self.station_var).pack(side='left')

        # ---- Electrode: free text ----
        tk.Label(body, text="Electrode:").grid(row=2, column=0, sticky='w', pady=4)
        self.electrode_var = tk.StringVar(value=g['electrode'])
        tk.Entry(body, textvariable=self.electrode_var, width=30).grid(row=2, column=1, sticky='w', pady=4)

        # ---- Measurement type: exactly one ----
        tk.Label(body, text="Measurement:").grid(row=3, column=0, sticky='nw', pady=4)
        mt = tk.Frame(body)
        mt.grid(row=3, column=1, sticky='w', pady=4)
        self.meas_var = tk.StringVar(value=g['meas'])
        for v in MEAS_TYPES:
            tk.Radiobutton(mt, text=v, value=v, variable=self.meas_var).pack(anchor='w')

        # ---- Laborant: pick from the saved list, or type a new name ----
        tk.Label(body, text="Laborant:").grid(row=4, column=0, sticky='w', pady=4)
        lb = tk.Frame(body)
        lb.grid(row=4, column=1, sticky='w', pady=4)
        self.lab_var = tk.StringVar(value=g['laborant'])
        self.lab_cb = ttk.Combobox(lb, textvariable=self.lab_var,
                                   values=get_laborants(), width=20)
        self.lab_cb.pack(side='left')
        tk.Button(lb, text="Edit list...", command=self._edit_laborants).pack(side='left', padx=(6, 0))
        tk.Label(body, text="A new name typed here is added to the list when you save.",
                 fg='#777', font=('TkDefaultFont', 8)).grid(row=5, column=1, sticky='w')

        # ---- live preview of the resulting name ----
        ttk.Separator(self, orient='horizontal').pack(fill='x', padx=10, pady=(6, 2))
        self.preview = tk.Label(self, text="", fg='#555', wraplength=460, justify='left')
        self.preview.pack(anchor='w', padx=10, pady=(0, 4))
        for var in (self.anion_var, self.station_var, self.electrode_var,
                    self.meas_var, self.lab_var):
            var.trace_add('write', lambda *_: self._update_preview())
        self._update_preview()

        # ---- analysis report (batch mode only) ----
        if report:
            ttk.Separator(self, orient='horizontal').pack(fill='x', padx=10, pady=(2, 2))
            rframe = tk.Frame(self)
            rframe.pack(fill='both', expand=True, padx=10, pady=(0, 4))
            tk.Label(rframe, text="Analysis report",
                     font=('TkDefaultFont', 9, 'bold')).pack(anchor='w')
            inner = tk.Frame(rframe)
            inner.pack(fill='both', expand=True)
            rsb = tk.Scrollbar(inner, orient='vertical')
            rsb.pack(side='right', fill='y')
            rsb_x = tk.Scrollbar(inner, orient='horizontal')
            rsb_x.pack(side='bottom', fill='x')
            rtxt = tk.Text(inner, height=18, width=110, wrap='none',
                           font=('Consolas', 8),
                           yscrollcommand=rsb.set, xscrollcommand=rsb_x.set)
            rtxt.pack(side='left', fill='both', expand=True)
            rsb.config(command=rtxt.yview)
            rsb_x.config(command=rtxt.xview)
            rtxt.insert('1.0', report)
            rtxt.config(state='disabled')   # read-only, still scrollable

        # ---- buttons ----
        btns = tk.Frame(self)
        btns.pack(fill='x', padx=10, pady=(4, 10))
        if mode == 'rename':
            tk.Button(btns, text="Rename", width=12, command=self._on_ok).pack(side='right')
            tk.Button(btns, text="Skip", width=10, command=self._on_skip).pack(side='right', padx=(0, 6))
            tk.Button(btns, text="Stop", width=10, command=self._on_cancel).pack(side='right', padx=(0, 6))
        elif mode == 'batch':
            tk.Button(btns, text="Approve & Save", width=15, command=self._on_ok).pack(side='right')
            tk.Button(btns, text="Skip", width=10, command=self._on_skip).pack(side='right', padx=(0, 6))
            tk.Button(btns, text="Stop", width=10, command=self._on_cancel).pack(side='right', padx=(0, 6))
        else:
            tk.Button(btns, text="Save file", width=12, command=self._on_ok).pack(side='right')
            tk.Button(btns, text="Cancel", width=10, command=self._on_cancel).pack(side='right', padx=(0, 6))

        self.bind('<Return>', lambda e: self._on_ok())
        self.bind('<Escape>', lambda e: self._on_cancel())

        # modal
        self.transient(parent)
        self.grab_set()
        self.protocol("WM_DELETE_WINDOW", self._on_cancel)
        self.update_idletasks()
        try:
            x = parent.winfo_rootx() + (parent.winfo_width() - self.winfo_width()) // 2
            y = parent.winfo_rooty() + (parent.winfo_height() - self.winfo_height()) // 3
            self.geometry(f"+{max(0, x)}+{max(0, y)}")
        except Exception:
            pass
        parent.wait_window(self)

    # ---------- helpers ----------
    def _anion_autocomplete(self, event):
        """Inline completion: type 'car' -> 'Carbonate' with the added part
        selected, so the next keystroke overwrites it."""
        if event.keysym in ('BackSpace', 'Delete', 'Left', 'Right', 'Up',
                            'Down', 'Escape', 'Return', 'Tab', 'Shift_L', 'Shift_R'):
            return
        typed = self.anion_var.get()
        if not typed:
            return
        for opt in ANIONS:
            if opt.lower().startswith(typed.lower()) and opt.lower() != typed.lower():
                self.anion_cb.delete(0, 'end')
                self.anion_cb.insert(0, opt)
                self.anion_cb.selection_range(len(typed), 'end')
                self.anion_cb.icursor(len(typed))
                break

    @staticmethod
    def _clean(text):
        """Make one name part safe for a Windows file name."""
        t = str(text).strip()
        for ch in '<>:"/\\|?*':
            t = t.replace(ch, '_')
        return t.replace(' ', '-')

    def _build_name(self):
        return "_".join(self._clean(v) for v in (
            self.anion_var.get(), self.station_var.get(), self.electrode_var.get(),
            self.meas_var.get(), self.lab_var.get()))

    def full_path(self):
        return os.path.join(self.folder, self._build_name() + '.csv')

    def new_mpt_path(self):
        return os.path.join(self.folder, self._build_name() + self.src_ext)

    def _update_preview(self):
        if self.mode == 'rename':
            self.preview.config(text="Will be renamed to:  " + self._build_name() + self.src_ext)
        else:
            self.preview.config(text="Will be saved as:  " + self._build_name() + ".csv")

    def _on_ok(self):
        missing = []
        if not self.anion_var.get().strip():
            missing.append("Anion")
        if not self.station_var.get():
            missing.append("Station")
        if not self.electrode_var.get().strip():
            missing.append("Electrode")
        if not self.meas_var.get():
            missing.append("Measurement type")
        if not self.lab_var.get().strip():
            missing.append("Laborant")
        if missing:
            messagebox.showwarning("Missing fields",
                                   "Please fill in:\n  - " + "\n  - ".join(missing),
                                   parent=self)
            return
        remember_laborant(self.lab_var.get())
        self.result = self._build_name()
        self.action = 'ok'
        self.destroy()

    def _edit_laborants(self):
        LaborantListDialog(self)
        self.lab_cb.config(values=get_laborants())
        try:
            self.grab_set()   # the list window took the grab; take it back
        except Exception:
            pass

    def _on_skip(self):
        self.result = None
        self.action = 'skip'
        self.destroy()

    def _on_cancel(self):
        self.result = None
        self.action = 'cancel'
        self.destroy()


# =================== Overlay window ===================
# =================== APPENDIX PAGE (for the CV Overlay window) ===================
# Graphs are "sent" from the overlay plot to an appendix page: a frozen copy of
# the figure is kept, so you can change anions / cycles and send more graphs.
# The page lays them out (default A4 portrait, 3 rows x 2 columns = 6 per page),
# every graph with an explanatory title and NO numbering, and is saved as an
# image (PNG / JPG / TIFF, one file per page) or as one multi-page PDF.

APPENDIX_PAGE_SIZES_CM = {
    "A4 portrait": (21.0, 29.7),
    "A4 landscape": (29.7, 21.0),
    "Letter portrait": (21.59, 27.94),
    "Letter landscape": (27.94, 21.59),
}
APPENDIX_HEADING_CM = 1.2  # height reserved for the optional page heading

# Measurement type (as in the file name) -> wording used in graph titles.
MEAS_DESCRIPTIONS = {
    'UAL-ACT-night': 'activation without urea',
    'UOS-night': 'activation with urea',
}


def _ap_wrap_title(text, width_in, fontsize):
    """Wrap a title so it fits the panel width (approximation)."""
    n = max(12, int(width_in * 72.0 / (max(float(fontsize), 1.0) * 0.58)))
    lines = []
    for line in str(text).split("\n"):
        lines.extend(textwrap.wrap(line, n, break_long_words=False) or [""])
    return "\n".join(lines)


def _ap_scale_fig_fonts(fig, scale):
    """Multiply every font size in the figure by `scale` (both overlay legends included)."""
    if abs(scale - 1.0) < 1e-9:
        return
    for ax in fig.axes:
        texts = [ax.title, ax.xaxis.label, ax.yaxis.label,
                 ax.xaxis.get_offset_text(), ax.yaxis.get_offset_text()]
        texts += list(ax.texts)
        for leg in [c for c in ax.get_children() if isinstance(c, Legend)]:
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


def appendix_render_panel(fig_bytes, title, w_px, h_px, dpi, font_scale):
    """Un-pickle a stored figure, resize it to the panel size, return an RGB array."""
    fig = pickle.loads(fig_bytes)
    fig.set_dpi(dpi)
    fig.set_size_inches(w_px / dpi, h_px / dpi)
    fig.set_facecolor("white")
    _ap_scale_fig_fonts(fig, font_scale)

    if fig.axes:
        main = fig.axes[0]
        txt = (title or "").strip("\n")
        if txt.strip():
            if not main.get_title():
                main.set_title(" ", fontsize=max(6.0, 14 * font_scale))
            fs = main.title.get_fontsize()
            main.title.set_text(_ap_wrap_title(txt, w_px / dpi, fs))
        else:
            main.title.set_text("")

    canvas = FigureCanvasAgg(fig)
    # the overlay uses tight_layout() once - redo it for the new (smaller) size
    eng = fig.get_layout_engine()
    if eng is None or type(eng).__name__ != 'ConstrainedLayoutEngine':
        try:
            fig.tight_layout(pad=0.4)
        except Exception:
            pass
    canvas.draw()
    return np.asarray(canvas.buffer_rgba())[:, :, :3].copy()


def _ap_render_heading(text, w_px, h_px, dpi):
    fig = Figure(figsize=(w_px / dpi, h_px / dpi), dpi=dpi)
    fig.set_facecolor("white")
    fig.text(0.5, 0.45, text, ha="center", va="center", fontsize=13, fontweight="bold")
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    return np.asarray(canvas.buffer_rgba())[:, :, :3].copy()


def _ap_paste(page, img, x, y, max_w, max_h):
    h = min(img.shape[0], max_h, page.shape[0] - y)
    w = min(img.shape[1], max_w, page.shape[1] - x)
    if h > 0 and w > 0:
        page[y:y + h, x:x + w] = img[:h, :w]


def appendix_compose_page(panels, cfg, dpi):
    """
    panels: list of dicts {"fig": bytes, "title": str} (at most rows*cols are used)
    cfg: dict(page, rows, cols, margin_cm, gap_cm, font_scale, heading)
    Returns a PIL RGB image of the full page at `dpi`.
    """
    pw_cm, ph_cm = APPENDIX_PAGE_SIZES_CM[cfg["page"]]

    def px(cm):
        return int(round(cm / 2.54 * dpi))

    W, H = px(pw_cm), px(ph_cm)
    page = np.full((H, W, 3), 255, dtype=np.uint8)

    m, g = px(cfg["margin_cm"]), px(cfg["gap_cm"])
    rows, cols = cfg["rows"], cfg["cols"]
    header = px(APPENDIX_HEADING_CM) if cfg["heading"] else 0

    area_w = W - 2 * m
    area_h = H - 2 * m - header
    cw = (area_w - g * (cols - 1)) // cols
    ch = (area_h - g * (rows - 1)) // rows
    if cw < 40 or ch < 40:
        raise ValueError("The graphs do not fit on the page - reduce the margin / gap "
                         "or the number of rows / columns.")

    if header:
        strip = _ap_render_heading(cfg["heading"], area_w, header, dpi)
        _ap_paste(page, strip, m, m, area_w, header)

    for k, p in enumerate(panels[:rows * cols]):
        r, c = divmod(k, cols)
        x = m + c * (cw + g)
        y = m + header + r * (ch + g)
        img = appendix_render_panel(p["fig"], p["title"], cw, ch, dpi, cfg["font_scale"])
        _ap_paste(page, img, x, y, cw, ch)

    return Image.fromarray(page)


def appendix_save_pages(pages, path, dpi):
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


# =================== GLOBAL Y RANGE (per measurement type) ===================
# Scanning heavy .mpt files is slow, so the min/max current of every cycle of
# every scanned file is remembered in a small JSON file. A file is only read
# again if it changed (size or modification time differ).
YRANGE_CACHE_FILE = os.path.join(os.path.expanduser('~'), '.cv_overlay_yrange_cache.json')


def _yr_load_cache():
    try:
        with open(YRANGE_CACHE_FILE, 'r', encoding='utf-8') as fh:
            d = json.load(fh)
        return d if isinstance(d, dict) else {}
    except Exception:
        return {}


def _yr_save_cache(cache):
    try:
        tmp = YRANGE_CACHE_FILE + '.tmp'
        with open(tmp, 'w', encoding='utf-8') as fh:
            json.dump(cache, fh)
        os.replace(tmp, YRANGE_CACHE_FILE)
    except Exception:
        pass  # the cache is only a speed-up


OVERLAY_SETTINGS_FILE = os.path.join(os.path.expanduser('~'), '.cv_overlay_settings.json')


def overlay_load_settings():
    try:
        with open(OVERLAY_SETTINGS_FILE, 'r', encoding='utf-8') as fh:
            d = json.load(fh)
        return d if isinstance(d, dict) else {}
    except Exception:
        return {}


def overlay_save_settings(d):
    try:
        tmp = OVERLAY_SETTINGS_FILE + '.tmp'
        with open(tmp, 'w', encoding='utf-8') as fh:
            json.dump(d, fh, indent=1, ensure_ascii=False)
        os.replace(tmp, OVERLAY_SETTINGS_FILE)
    except Exception:
        pass


_SCROLLERS = []      # [(outer frame, canvas)] of every scrollable panel
_SCROLL_BOUND = []   # non-empty once the shared wheel handler is installed


def _scroll_wheel(event):
    try:
        w = event.widget.winfo_containing(event.x_root, event.y_root)
    except Exception:
        return
    while w is not None:
        if isinstance(w, (tk.Listbox, tk.Text)):
            return            # lists scroll themselves
        for outer, canvas in list(_SCROLLERS):
            try:
                if not outer.winfo_exists():
                    _SCROLLERS.remove((outer, canvas))
                    continue
            except Exception:
                continue
            if w is outer:
                if canvas.yview() == (0.0, 1.0):
                    return
                num = getattr(event, 'num', None)
                if num == 4:
                    step = -1
                elif num == 5:
                    step = 1
                else:
                    step = -1 if event.delta > 0 else 1
                canvas.yview_scroll(step * 3, 'units')
                return
        w = getattr(w, 'master', None)


def make_scrollable(parent):
    """Return an inner Frame inside `parent` that scrolls vertically
    (scrollbar + mouse wheel while the pointer is over it)."""
    outer = tk.Frame(parent)
    outer.pack(fill='both', expand=True)
    canvas = tk.Canvas(outer, highlightthickness=0, borderwidth=0)
    vsb = tk.Scrollbar(outer, orient='vertical', command=canvas.yview)
    canvas.configure(yscrollcommand=vsb.set)
    vsb.pack(side='right', fill='y')
    canvas.pack(side='left', fill='both', expand=True)
    inner = tk.Frame(canvas)
    win = canvas.create_window((0, 0), window=inner, anchor='nw')
    inner.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
    canvas.bind('<Configure>', lambda e: canvas.itemconfigure(win, width=e.width))

    _SCROLLERS.append((outer, canvas))
    if not _SCROLL_BOUND:
        # one shared handler (added, not replacing other wheel bindings); it scrolls
        # the scrollable panel that is under the mouse pointer, if any
        for seq in ('<MouseWheel>', '<Button-4>', '<Button-5>'):
            outer.bind_all(seq, _scroll_wheel, add='+')
        _SCROLL_BOUND.append(True)
    return inner


def nice_step(span):
    """A round step of about 1/20 of `span` (1, 2 or 5 x 10^n)."""
    span = abs(float(span))
    if not np.isfinite(span) or span <= 0:
        return 1.0
    raw = span / 20.0
    mag = 10 ** np.floor(np.log10(raw))
    for m in (1, 2, 5, 10):
        if raw <= m * mag:
            return float(m * mag)
    return float(10 * mag)


def yr_cycle_extrema(I, cyc):
    """{cycle: (min I, max I)} over the finite points of every cycle."""
    I = np.asarray(I, dtype=float)
    cyc = np.asarray(cyc)
    out = {}
    for c in np.unique(cyc):
        v = I[cyc == c]
        v = v[np.isfinite(v)]
        if v.size:
            out[int(c)] = (float(v.min()), float(v.max()))
    return out


def yr_nice_limits(lo, hi, pad_frac=0.03):
    """Pad the range a little and round it outward to 'nice' tick values."""
    from matplotlib.ticker import MaxNLocator
    span = hi - lo
    pad = span * pad_frac if span > 0 else (abs(hi) * 0.05 or 1.0)
    ticks = MaxNLocator(nbins=12, steps=[1, 2, 2.5, 5, 10]).tick_values(lo - pad, hi + pad)
    return float(ticks[0]), float(ticks[-1])


class AppendixPanel:
    """The 'Appendix page' tab inside the CV Overlay window.

    The list of graphs lives on the main app (app.appendix_panels), so it is kept
    when the overlay window is closed and opened again.
    """

    def __init__(self, owner, parent, panels):
        self.owner = owner          # the OverlayWindow (for dialogs / after())
        self.panels = panels        # shared list of {"fig": bytes, "title": str, "orig": str}
        self.page_idx = 0
        self._loading = False
        self._after_id = None

        self.v_page = tk.StringVar(master=owner, value="A4 portrait")
        self.v_rows = tk.StringVar(master=owner, value="3")
        self.v_cols = tk.StringVar(master=owner, value="2")
        self.v_margin = tk.StringVar(master=owner, value="2.0")
        self.v_gap = tk.StringVar(master=owner, value="0.5")
        self.v_fscale = tk.StringVar(master=owner, value="0.75")
        self.v_heading = tk.StringVar(master=owner, value="")
        self.v_dpi = tk.StringVar(master=owner, value="300")
        self.v_title = tk.StringVar(master=owner, value="")

        self._build(parent)

        for v in (self.v_page, self.v_rows, self.v_cols, self.v_margin, self.v_gap,
                  self.v_fscale, self.v_heading):
            v.trace_add("write", self._schedule_preview)
        self.v_title.trace_add("write", self._on_title_edit)

        self._refresh_list()
        owner.after(150, self._render_preview)

    # ---------- settings ----------
    def _cfg(self):
        def num(var, default, lo, hi):
            try:
                v = float(var.get())
            except Exception:
                v = default
            return max(lo, min(hi, v))

        page = self.v_page.get()
        if page not in APPENDIX_PAGE_SIZES_CM:
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

    def _capacity(self, cfg=None):
        cfg = cfg or self._cfg()
        return cfg["rows"] * cfg["cols"]

    def _num_pages(self, cfg=None):
        return max(1, -(-len(self.panels) // self._capacity(cfg)))

    # ---------- UI ----------
    def _build(self, parent):
        pw = tk.PanedWindow(parent, orient='horizontal', sashrelief='raised', sashwidth=5)
        pw.pack(fill='both', expand=True)
        left = tk.Frame(pw)
        right = tk.Frame(pw)
        pw.add(left, minsize=250)
        pw.add(right, minsize=300)

        tk.Label(left, wraplength=240, justify='left', fg='#555',
                 text="Plot an overlay, then press “→ Appendix page” under the plot "
                      "to send it here.").pack(anchor='w', padx=6, pady=(6, 4))

        tk.Label(left, text="Graphs (page order):",
                 font=('TkDefaultFont', 9, 'bold')).pack(anchor='w', padx=6)
        lbox = tk.Frame(left)
        lbox.pack(fill='x', padx=6)
        sb = tk.Scrollbar(lbox, orient='vertical')
        sb.pack(side='right', fill='y')
        self.lb = tk.Listbox(lbox, height=8, exportselection=False, yscrollcommand=sb.set)
        self.lb.pack(side='left', fill='x', expand=True)
        sb.config(command=self.lb.yview)
        self.lb.bind("<<ListboxSelect>>", self._on_select)

        b1 = tk.Frame(left)
        b1.pack(fill='x', padx=6, pady=(2, 6))
        tk.Button(b1, text="▲", width=3, command=lambda: self._move(-1)).pack(side='left')
        tk.Button(b1, text="▼", width=3, command=lambda: self._move(+1)).pack(side='left', padx=(3, 0))
        tk.Button(b1, text="Remove", command=self._remove).pack(side='left', padx=(8, 0))
        tk.Button(b1, text="Clear all", command=self._clear).pack(side='left', padx=(3, 0))

        tk.Label(left, text="Title of selected graph (empty = no title):",
                 font=('TkDefaultFont', 9, 'bold')).pack(anchor='w', padx=6)
        tk.Entry(left, textvariable=self.v_title).pack(fill='x', padx=6)
        tk.Button(left, text="↺ Original title", command=self._reset_title).pack(anchor='w', padx=6, pady=(2, 6))

        g = tk.LabelFrame(left, text="Page layout")
        g.pack(fill='x', padx=6, pady=(0, 6))

        def row(r, label, widget):
            tk.Label(g, text=label).grid(row=r, column=0, sticky='w', padx=4, pady=1)
            widget.grid(row=r, column=1, sticky='w', padx=4, pady=1)

        row(0, "Page", ttk.Combobox(g, state='readonly', width=15,
                                    values=list(APPENDIX_PAGE_SIZES_CM), textvariable=self.v_page))
        row(1, "Rows", tk.Spinbox(g, from_=1, to=6, width=6, textvariable=self.v_rows))
        row(2, "Columns", tk.Spinbox(g, from_=1, to=4, width=6, textvariable=self.v_cols))
        row(3, "Margin (cm)", tk.Spinbox(g, from_=0, to=6, increment=0.25, format="%.2f",
                                         width=6, textvariable=self.v_margin))
        row(4, "Gap (cm)", tk.Spinbox(g, from_=0, to=3, increment=0.1, format="%.2f",
                                      width=6, textvariable=self.v_gap))
        row(5, "Font scale", tk.Spinbox(g, from_=0.3, to=2.5, increment=0.05, format="%.2f",
                                        width=6, textvariable=self.v_fscale))
        row(6, "Page heading", tk.Entry(g, width=18, textvariable=self.v_heading))
        row(7, "Save DPI", tk.Spinbox(g, from_=100, to=1200, increment=50,
                                      width=6, textvariable=self.v_dpi))
        # the Spinboxes reformat their text on creation -> put the defaults back
        self.v_rows.set("3"); self.v_cols.set("2")
        self.v_margin.set("2.00"); self.v_gap.set("0.50")
        self.v_fscale.set("0.75"); self.v_dpi.set("300")

        tk.Label(left, wraplength=240, justify='left', fg='#666', font=('TkDefaultFont', 8),
                 text="Font scale multiplies the overlay font sizes when the graph "
                      "is shrunk onto the page.").pack(anchor='w', padx=6)

        tk.Button(left, text="Save page(s) as image…",
                  command=self.save_pages).pack(fill='x', padx=6, pady=8)

        nav = tk.Frame(right)
        nav.pack(side='top', fill='x', pady=(4, 2))
        tk.Button(nav, text="◀", width=3, command=lambda: self._goto(-1)).pack(side='left', padx=(4, 0))
        self.page_lbl = tk.Label(nav, text="Page 1 / 1", width=12)
        self.page_lbl.pack(side='left', padx=4)
        tk.Button(nav, text="▶", width=3, command=lambda: self._goto(+1)).pack(side='left')
        tk.Label(nav, text="(preview - the saved file uses the full DPI)",
                 fg='#666').pack(side='left', padx=8)

        self.fig = Figure(figsize=(5, 7), dpi=100)
        self.fig.patch.set_facecolor("#dcdcdc")
        self.ax = self.fig.add_axes([0.01, 0.01, 0.98, 0.98])
        self.ax.axis("off")
        self.canvas = FigureCanvasTkAgg(self.fig, right)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        NavigationToolbar2Tk(self.canvas, right).update()

    # ---------- list handling ----------
    def _refresh_list(self, select=None):
        cap = self._capacity()
        self.lb.delete(0, 'end')
        for i, p in enumerate(self.panels):
            pg, pos = divmod(i, cap)
            t = (p["title"] or "").replace("\n", " ").strip() or "(no title)"
            self.lb.insert('end', f"p{pg + 1}·{pos + 1}   {t}")
        if select is not None and 0 <= select < len(self.panels):
            self.lb.selection_clear(0, 'end')
            self.lb.selection_set(select)
            self.lb.see(select)
        self.owner._update_appendix_tab_text()

    def _selected(self):
        sel = self.lb.curselection()
        return sel[0] if sel else None

    def _set_title_box(self, text):
        self._loading = True
        try:
            self.v_title.set(text)
        finally:
            self._loading = False

    def _on_select(self, _event=None):
        i = self._selected()
        self._set_title_box(self.panels[i]["title"] if i is not None else "")
        if i is not None:
            self.page_idx = i // self._capacity()
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
        if i is not None:
            self.v_title.set(self.panels[i]["orig"])

    def _move(self, delta):
        i = self._selected()
        if i is None:
            return
        j = i + delta
        if not (0 <= j < len(self.panels)):
            return
        self.panels[i], self.panels[j] = self.panels[j], self.panels[i]
        self.page_idx = j // self._capacity()
        self._refresh_list(select=j)
        self._schedule_preview()

    def _remove(self):
        i = self._selected()
        if i is None:
            return
        del self.panels[i]
        self._refresh_list(select=min(i, len(self.panels) - 1) if self.panels else None)
        self._on_select()
        self._schedule_preview()

    def _clear(self):
        if not self.panels:
            return
        if not messagebox.askyesno("Clear all", "Remove all graphs from the appendix page?",
                                   parent=self.owner):
            return
        self.panels.clear()
        self.page_idx = 0
        self._set_title_box("")
        self._refresh_list()
        self._schedule_preview()

    # ---------- sending ----------
    def add_figure(self, fig, title):
        try:
            data = pickle.dumps(fig)
        except Exception as e:
            messagebox.showerror("Could not copy the graph", f"{type(e).__name__}: {e}",
                                 parent=self.owner)
            return False
        self.panels.append({"fig": data, "title": title, "orig": title})
        i = len(self.panels) - 1
        self.page_idx = i // self._capacity()
        self._refresh_list(select=i)
        self._on_select()
        self._schedule_preview()
        return True

    # ---------- preview ----------
    def _goto(self, delta):
        n = self._num_pages()
        self.page_idx = max(0, min(self.page_idx + delta, n - 1))
        self._render_preview()

    def _schedule_preview(self, *_):
        if self._after_id is not None:
            try:
                self.owner.after_cancel(self._after_id)
            except Exception:
                pass
        self._after_id = self.owner.after(300, self._render_preview)

    def _render_preview(self):
        self._after_id = None
        try:
            if not self.owner.winfo_exists():
                return
        except Exception:
            return
        ax = self.ax
        ax.clear()
        ax.axis("off")
        cfg = self._cfg()
        cap = self._capacity(cfg)
        npages = self._num_pages(cfg)
        self.page_idx = max(0, min(self.page_idx, npages - 1))
        self.page_lbl.config(text=f"Page {self.page_idx + 1} / {npages}")

        if not self.panels:
            ax.text(0.5, 0.5, "No graphs yet.\nPlot an overlay, then press\n“→ Appendix page”.",
                    ha="center", va="center", fontsize=12, color="#555")
            self.canvas.draw_idle()
            return

        chunk = self.panels[self.page_idx * cap:(self.page_idx + 1) * cap]
        try:
            img = appendix_compose_page(chunk, cfg, dpi=100)
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
            messagebox.showinfo("Nothing to save",
                                "Send at least one graph to the appendix page first.",
                                parent=self.owner)
            return
        path = filedialog.asksaveasfilename(
            parent=self.owner, title="Save appendix page",
            defaultextension=".png", initialfile="appendix_page.png",
            initialdir=(EXPORT_FOLDER if os.path.isdir(EXPORT_FOLDER) else None),
            filetypes=[("PNG", "*.png"), ("PDF (all pages in one file)", "*.pdf"),
                       ("JPEG", "*.jpg"), ("TIFF", "*.tif"), ("All files", "*.*")])
        if not path:
            return
        cfg = self._cfg()
        cap = self._capacity(cfg)
        dpi = cfg["dpi"]
        self.owner.config(cursor="watch")
        self.owner.update_idletasks()
        try:
            pages = [appendix_compose_page(self.panels[i:i + cap], cfg, dpi)
                     for i in range(0, len(self.panels), cap)]
            written = appendix_save_pages(pages, path, dpi)
        except Exception as e:
            messagebox.showerror("Save failed", f"{type(e).__name__}: {e}", parent=self.owner)
            return
        finally:
            self.owner.config(cursor="")
        messagebox.showinfo("Appendix page", "Saved:\n" + "\n".join(written), parent=self.owner)


class OverlayWindow(tk.Toplevel):
    """Separate window for overlaying CV curves.

    Files are selected by anion (one or many) and measurement type (exactly
    one), both read from the .mpt file names, so no file has to be loaded
    just to appear in the lists. Only the files that actually match are
    parsed, which keeps this fast on a large library.
    """

    # Per-cycle styles, in the order cycles first appear. Markers are drawn
    # sparsely (markevery) because a CV sweep has thousands of points and a
    # marker on each one would be a solid blob.
    CYCLE_STYLES = [
        ('-',  None),
        ('--', None),
        (':',  None),
        ('-.', None),
        ('-',  'o'),
        ('--', 's'),
        (':',  '^'),
        ('-.', 'D'),
        ((0, (5, 1)), None),
        ((0, (3, 1, 1, 1)), None),
        ('-',  'v'),
        ('--', '*'),
    ]

    LEGEND_LOCS = ['best', 'upper right', 'upper left', 'lower left',
                   'lower right', 'right', 'center left', 'center right',
                   'lower center', 'upper center', 'center']

    def __init__(self, app):
        super().__init__(app.root)
        self.app = app
        self.title("CV Overlay")
        self.geometry("1150x720")

        # path -> parsed name parts, for every .mpt listed in the tree
        self.candidates = app._overlay_candidates()
        self.matched = []      # [(node, path, meta)] after "Find files"
        self.series = []       # [(label, E, I)] currently plotted

        main = tk.PanedWindow(self, orient='horizontal', sashrelief='raised', sashwidth=6)
        main.pack(fill='both', expand=True)

        left = tk.Frame(main)
        main.add(left, minsize=260)
        right = tk.Frame(main)
        main.add(right, minsize=520)

        # ---------------- selection ----------------
        nb = ttk.Notebook(left)
        nb.pack(fill='both', expand=True, padx=4, pady=4)
        sel_outer = tk.Frame(nb)
        nb.add(sel_outer, text="Select")
        fmt_outer = tk.Frame(nb)
        nb.add(fmt_outer, text="Appearance")
        sel_tab = make_scrollable(sel_outer)   # both tabs scroll when they are taller than the window
        fmt_tab = make_scrollable(fmt_outer)

        tk.Label(sel_tab, text="Anion (one or more):",
                 font=('TkDefaultFont', 9, 'bold')).pack(anchor='w', padx=6, pady=(6, 0))
        anion_box = tk.Frame(sel_tab)
        anion_box.pack(fill='x', padx=6)
        asb = tk.Scrollbar(anion_box, orient='vertical')
        asb.pack(side='right', fill='y')
        self.anion_list = tk.Listbox(anion_box, selectmode='extended', height=7,
                                     exportselection=False, yscrollcommand=asb.set)
        self.anion_list.pack(side='left', fill='x', expand=True)
        asb.config(command=self.anion_list.yview)

        avail_anions = sorted({m['anion'] for _n, _p, m in self.candidates if m['anion']})
        for a in avail_anions:
            n = sum(1 for _n, _p, m in self.candidates if m['anion'] == a)
            self.anion_list.insert('end', f"{a}  ({n})")
        self._anion_names = avail_anions

        tk.Label(sel_tab, text="Measurement type (one):",
                 font=('TkDefaultFont', 9, 'bold')).pack(anchor='w', padx=6, pady=(8, 0))
        self.meas_var = tk.StringVar(value='')
        avail_meas = sorted({m['meas'] for _n, _p, m in self.candidates if m['meas']})
        for mt in (avail_meas or MEAS_TYPES):
            n = sum(1 for _n, _p, m in self.candidates if m['meas'] == mt)
            tk.Radiobutton(sel_tab, text=f"{mt}  ({n})", value=mt,
                           variable=self.meas_var).pack(anchor='w', padx=12)
        if avail_meas:
            self.meas_var.set(avail_meas[0])

        tk.Button(sel_tab, text="Find matching files",
                  command=self.find_files).pack(fill='x', padx=6, pady=(10, 2))
        self.match_lbl = tk.Label(sel_tab, text="", fg='#555', wraplength=240, justify='left')
        self.match_lbl.pack(anchor='w', padx=6)

        tk.Label(sel_tab, text="Cycle number(s):",
                 font=('TkDefaultFont', 9, 'bold')).pack(anchor='w', padx=6, pady=(8, 0))
        self.cycles_var = tk.StringVar(value='')
        tk.Entry(sel_tab, textvariable=self.cycles_var).pack(fill='x', padx=6)
        tk.Label(sel_tab, text="e.g.  3    or   1,3,5    or   1-5,8",
                 fg='#777', font=('TkDefaultFont', 8)).pack(anchor='w', padx=6)
        self.cyc_hint = tk.Label(sel_tab, text="", fg='#555',
                                 wraplength=240, justify='left',
                                 font=('TkDefaultFont', 8))
        self.cyc_hint.pack(anchor='w', padx=6)

        # ---- common Y axis for all graphs of one measurement type ----
        yfr = tk.LabelFrame(sel_tab, text="Common Y axis (this measurement type)")
        yfr.pack(fill='x', padx=6, pady=(8, 0))
        ybtn = tk.Frame(yfr)
        ybtn.pack(fill='x', padx=4, pady=(4, 2))
        tk.Button(ybtn, text="Scan for global Y range",
                  command=self.scan_global_y).pack(side='left')
        tk.Button(ybtn, text="Stop", command=self._stop_y_scan).pack(side='left', padx=4)
        tk.Button(ybtn, text="Auto Y", command=self._clear_y_limits).pack(side='left')
        tk.Button(yfr, text="Save current Y min / max as default for this type",
                  command=self._save_y_default_manual).pack(fill='x', padx=4, pady=(0, 2))
        self.yscan_lbl = tk.Label(yfr, text="Scans ALL files of the chosen measurement type "
                                             "(every anion), for the cycle(s) entered above "
                                             "(empty = all cycles).",
                                  fg='#555', wraplength=240, justify='left',
                                  font=('TkDefaultFont', 8))
        self.yscan_lbl.pack(anchor='w', padx=4, pady=(0, 2))
        tk.Label(yfr, text="The result fills Y min / Y max under Appearance - "
                           "you can also type them there yourself.",
                 fg='#777', wraplength=240, justify='left',
                 font=('TkDefaultFont', 8)).pack(anchor='w', padx=4, pady=(0, 4))
        self._y_scan_stop = False

        # ---- files going into the plot ----
        tk.Label(sel_tab, text="Files in the plot:",
                 font=('TkDefaultFont', 9, 'bold')).pack(anchor='w', padx=6, pady=(8, 0))
        flist_box = tk.Frame(sel_tab)
        flist_box.pack(fill='both', expand=True, padx=6)
        fsb = tk.Scrollbar(flist_box, orient='vertical')
        fsb.pack(side='right', fill='y')
        fsb_x = tk.Scrollbar(flist_box, orient='horizontal')
        fsb_x.pack(side='bottom', fill='x')
        self.file_list = tk.Listbox(flist_box, height=7, activestyle='none',
                                    yscrollcommand=fsb.set, xscrollcommand=fsb_x.set)
        self.file_list.pack(side='left', fill='both', expand=True)
        fsb.config(command=self.file_list.yview)
        fsb_x.config(command=self.file_list.xview)

        tk.Button(sel_tab, text="Plot overlay",
                  command=self.plot_overlay).pack(fill='x', padx=6, pady=6)

        # ---------------- appearance ----------------
        self._build_format_tab(fmt_tab)

        # ---------------- right side: Overlay | Appendix page ----------------
        self.right_nb = ttk.Notebook(right)
        self.right_nb.pack(fill='both', expand=True)
        plot_tab = tk.Frame(self.right_nb)
        self.right_nb.add(plot_tab, text="Overlay")
        self.appendix_tab = tk.Frame(self.right_nb)
        self.right_nb.add(self.appendix_tab, text="Appendix page")
        right = plot_tab  # the plot and its save bar go into the 'Overlay' tab

        # ---------------- figure ----------------
        self.fig = Figure(figsize=(7, 5))
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, right)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        NavigationToolbar2Tk(self.canvas, right).update()

        save_bar = tk.Frame(right)
        save_bar.pack(fill='x', pady=(0, 4))
        tk.Button(save_bar, text="Save figure...", command=self.save_figure).pack(side='left', padx=4)
        tk.Label(save_bar, text="DPI:").pack(side='left')
        self.dpi_var = tk.StringVar(value='300')
        tk.Entry(save_bar, textvariable=self.dpi_var, width=5).pack(side='left')
        tk.Button(save_bar, text="Copy data to clipboard",
                  command=self.copy_data).pack(side='left', padx=8)
        self.btn_appendix = tk.Button(save_bar, text="→ Appendix page",
                                      command=self.send_to_appendix)
        self.btn_appendix.pack(side='left', padx=8)

        # The appendix list is kept on the main app, so it survives closing this window.
        if not hasattr(app, 'appendix_panels'):
            app.appendix_panels = []
        self._plotted_meas = ''
        self.appendix = AppendixPanel(self, self.appendix_tab, app.appendix_panels)

        # ---- automatic updates: every change of the selection re-plots ----
        self._last_auto_title = ''
        self._fmt_after = None
        self._replot_after = None
        self.anion_list.bind('<<ListboxSelect>>', lambda e: self._schedule_replot(150), add='+')
        self.meas_var.trace_add('write', lambda *_: (self._load_y_default(), self._schedule_replot(150)))
        self.cycles_var.trace_add('write', lambda *_: (self._load_y_default(), self._schedule_replot(700)))
        self.after(50, self._load_y_default)

        self._redraw_empty()

    # ---------------- appearance tab ----------------
    def _build_format_tab(self, tab):
        r = 0

        def entry(label, var, width=22):
            nonlocal r
            tk.Label(tab, text=label).grid(row=r, column=0, sticky='w', padx=6, pady=2)
            tk.Entry(tab, textvariable=var, width=width).grid(row=r, column=1, sticky='w', padx=4)
            r += 1

        def num_entry(label, var, step, width=7, lo=None):
            """Entry with − / + buttons. `step` is a number or a function returning one."""
            nonlocal r
            tk.Label(tab, text=label).grid(row=r, column=0, sticky='w', padx=6, pady=2)
            fr = tk.Frame(tab)
            fr.grid(row=r, column=1, sticky='w', padx=4)
            tk.Button(fr, text="−", width=2, pady=0,
                      command=lambda: self._bump(var, step, -1, lo)).pack(side='left')
            tk.Entry(fr, textvariable=var, width=width, justify='center').pack(side='left', padx=2)
            tk.Button(fr, text="+", width=2, pady=0,
                      command=lambda: self._bump(var, step, +1, lo)).pack(side='left')
            r += 1

        self.title_var = tk.StringVar(value='')
        self.xlabel_var = tk.StringVar(value='E vs RHE (V)')
        self.ylabel_var = tk.StringVar(value='Current (mA)')
        entry("Title:", self.title_var, 26)
        tk.Label(tab, text="(auto: anion + measurement type; type your own to override)",
                 fg='#777', font=('TkDefaultFont', 8)).grid(
            row=r, column=0, columnspan=2, sticky='w', padx=6); r += 1
        tk.Button(tab, text="↺ Auto title", command=self._reset_auto_title).grid(
            row=r, column=1, sticky='w', padx=4, pady=(0, 4)); r += 1
        entry("X label:", self.xlabel_var)
        entry("Y label:", self.ylabel_var)

        tk.Label(tab, text="Font sizes", font=('TkDefaultFont', 9, 'bold')).grid(
            row=r, column=0, sticky='w', padx=6, pady=(8, 2)); r += 1
        self.fs_title = tk.StringVar(value='14')
        self.fs_label = tk.StringVar(value='12')
        self.fs_tick = tk.StringVar(value='10')
        self.fs_legend = tk.StringVar(value='9')
        num_entry("  Title:", self.fs_title, 1, 5, lo=1)
        num_entry("  Axis labels:", self.fs_label, 1, 5, lo=1)
        num_entry("  Ticks:", self.fs_tick, 1, 5, lo=1)
        num_entry("  Legend:", self.fs_legend, 1, 5, lo=1)

        tk.Label(tab, text="Lines & legend", font=('TkDefaultFont', 9, 'bold')).grid(
            row=r, column=0, sticky='w', padx=6, pady=(8, 2)); r += 1
        self.lw_var = tk.StringVar(value='1.8')
        num_entry("  Line width:", self.lw_var, 0.2, 5, lo=0.1)
        self.ms_var = tk.StringVar(value='5')
        num_entry("  Marker size:", self.ms_var, 1, 5, lo=0)

        self.show_legend = tk.BooleanVar(value=True)
        tk.Checkbutton(tab, text="Show legend", variable=self.show_legend).grid(
            row=r, column=0, columnspan=2, sticky='w', padx=6); r += 1

        tk.Label(tab, text="  Legend position:").grid(row=r, column=0, sticky='w', padx=6)
        self.legend_loc = tk.StringVar(value='best')
        ttk.Combobox(tab, textvariable=self.legend_loc, values=self.LEGEND_LOCS,
                     width=15, state='readonly').grid(row=r, column=1, sticky='w', padx=4); r += 1

        self.show_cycle_legend = tk.BooleanVar(value=True)
        tk.Checkbutton(tab, text="Show cycle legend (line styles)",
                       variable=self.show_cycle_legend).grid(
            row=r, column=0, columnspan=2, sticky='w', padx=6); r += 1

        tk.Label(tab, text="  Cycle legend position:").grid(row=r, column=0, sticky='w', padx=6)
        self.cycle_legend_loc = tk.StringVar(value='lower right')
        ttk.Combobox(tab, textvariable=self.cycle_legend_loc, values=self.LEGEND_LOCS,
                     width=15, state='readonly').grid(row=r, column=1, sticky='w', padx=4); r += 1

        self.show_grid = tk.BooleanVar(value=True)
        tk.Checkbutton(tab, text="Show grid", variable=self.show_grid).grid(
            row=r, column=0, columnspan=2, sticky='w', padx=6); r += 1

        tk.Label(tab, text="Axis limits (blank = auto)",
                 font=('TkDefaultFont', 9, 'bold')).grid(
            row=r, column=0, sticky='w', padx=6, pady=(8, 2)); r += 1
        self.xmin_var = tk.StringVar(value='')
        self.xmax_var = tk.StringVar(value='')
        self.ymin_var = tk.StringVar(value='')
        self.ymax_var = tk.StringVar(value='')
        # step follows the current axis span (about 1/20 of it)
        xstep = lambda: nice_step(np.diff(self.ax.get_xlim())[0])
        ystep = lambda: nice_step(np.diff(self.ax.get_ylim())[0])
        num_entry("  X min:", self.xmin_var, xstep, 8)
        num_entry("  X max:", self.xmax_var, xstep, 8)
        num_entry("  Y min:", self.ymin_var, ystep, 8)
        num_entry("  Y max:", self.ymax_var, ystep, 8)

        # The legend now shows one entry per anion, so these labels apply to
        # the clipboard export rather than to the legend.
        tk.Label(tab, text="Data-export labels", font=('TkDefaultFont', 9, 'bold')).grid(
            row=r, column=0, sticky='w', padx=6, pady=(8, 2)); r += 1
        # Default includes the electrode number: several files can share an
        # anion, and 'Nitrite c3' three times over is a useless legend.
        self.label_mode = tk.StringVar(value='anion+electrode+cycle')
        for val, txt in [('anion+cycle', 'Anion + cycle'),
                         ('anion+electrode+cycle', 'Anion + electrode + cycle'),
                         ('file+cycle', 'File name + cycle')]:
            tk.Radiobutton(tab, text=txt, value=val, variable=self.label_mode).grid(
                row=r, column=0, columnspan=2, sticky='w', padx=12); r += 1

        tk.Label(tab, text="Changes are applied automatically.", fg='#777',
                 font=('TkDefaultFont', 8)).grid(
            row=r, column=0, columnspan=2, sticky='w', padx=6, pady=(10, 0)); r += 1
        tk.Button(tab, text="Apply", command=self.apply_format).grid(
            row=r, column=0, columnspan=2, sticky='we', padx=6, pady=(2, 10))

        # ---- every appearance change redraws (no file loading - just restyling) ----
        for v in (self.title_var, self.xlabel_var, self.ylabel_var,
                  self.fs_title, self.fs_label, self.fs_tick, self.fs_legend,
                  self.lw_var, self.ms_var, self.show_legend, self.legend_loc,
                  self.show_cycle_legend, self.cycle_legend_loc, self.show_grid,
                  self.xmin_var, self.xmax_var, self.ymin_var, self.ymax_var):
            v.trace_add('write', self._schedule_format)

    def _bump(self, var, step, sign, lo=None):
        """− / + button: change a numeric field by one step. An empty axis-limit
        field starts from the limit currently shown on the plot."""
        try:
            st = float(step() if callable(step) else step)
        except Exception:
            st = 1.0
        cur = self._num(var)
        if cur is None:
            axis_now = {id(self.xmin_var): lambda: self.ax.get_xlim()[0],
                        id(self.xmax_var): lambda: self.ax.get_xlim()[1],
                        id(self.ymin_var): lambda: self.ax.get_ylim()[0],
                        id(self.ymax_var): lambda: self.ax.get_ylim()[1]}.get(id(var))
            if axis_now is None:
                return
            cur = float(axis_now())
            cur = round(cur / st) * st          # snap to the step grid
        new = cur + sign * st
        if lo is not None:
            new = max(lo, new)
        new = round(new, 10)
        var.set(f"{new:g}")

    # ---------------- data ----------------
    def find_files(self, silent=False):
        picked = [self._anion_names[i] for i in self.anion_list.curselection()]
        meas = self.meas_var.get()
        if not picked:
            if not silent:
                messagebox.showwarning("Overlay", "Select at least one anion.", parent=self)
            self.matched = []
            return
        if not meas:
            if not silent:
                messagebox.showwarning("Overlay", "Select a measurement type.", parent=self)
            self.matched = []
            return

        self.matched = [(n, p, m) for n, p, m in self.candidates
                        if m['anion'] in picked and m['meas'] == meas]
        if not self.matched:
            self.match_lbl.config(text="No files match that combination.", fg='#b00020')
            self.cyc_hint.config(text="")
            self._show_files([], "matched")
            return

        # Matching is done on file NAMES only - nothing is read from disk
        # here. The files are opened in plot_overlay, once the cycle numbers
        # are known, so this stays instant however large the folder is.
        self.match_lbl.config(text=f"{len(self.matched)} file(s) matched.", fg='#555')
        self.cyc_hint.config(text="Enter cycle number(s), then press Plot overlay.")
        self._show_files([p for _n, p, _m in self.matched], "matched")

    def _show_files(self, paths, kind):
        """Fill the 'Files in the plot' box."""
        self.file_list.delete(0, 'end')
        if not paths:
            self.file_list.insert('end', "(none)")
            return
        for p in paths:
            self.file_list.insert('end', os.path.basename(p))
        self.file_list.insert('end', f"--- {len(paths)} file(s) {kind} ---")

    def _parse_cycles(self):
        """Read the cycle box: '1,3,5' or '1-5,8' or '3'. Returns a sorted
        list, or None if something could not be understood."""
        txt = self.cycles_var.get().replace(';', ',')
        out = set()
        for part in txt.split(','):
            part = part.strip()
            if not part:
                continue
            if '-' in part:
                a, _sep, b = part.partition('-')
                try:
                    a, b = int(a.strip()), int(b.strip())
                except ValueError:
                    return None
                if a > b:
                    a, b = b, a
                out.update(range(a, b + 1))
            else:
                try:
                    out.add(int(part))
                except ValueError:
                    return None
        return sorted(out)

    def _label_for(self, meta, path, cyc):
        mode = self.label_mode.get()
        if mode == 'file+cycle':
            return f"{os.path.splitext(os.path.basename(path))[0]}  c{cyc}"
        if mode == 'anion+electrode+cycle':
            return f"{meta['anion'] or '?'} {meta['electrode'] or ''} c{cyc}".replace('  ', ' ')
        return f"{meta['anion'] or '?'}  c{cyc}"

    def plot_overlay(self, silent=False):
        if not self.matched:
            if not silent:
                messagebox.showwarning("Overlay", "Press 'Find matching files' first.", parent=self)
            return

        want = self._parse_cycles()
        if want is None:
            if silent:
                return
            messagebox.showwarning(
                "Overlay",
                "Could not read the cycle numbers.\n\n"
                "Use a number, a comma-separated list, or a range:\n"
                "   3        1,3,5        1-5,8",
                parent=self)
            return
        if not want:
            if not silent:
                messagebox.showwarning("Overlay", "Enter at least one cycle number.", parent=self)
            return

        try:
            pH = float(self.app.ph_var.get())
        except Exception:
            pH = 14.0
        ref = self.app.ref_var.get()

        # ---- files are opened HERE, now that the cycles are known ----
        self.series = []
        used_paths = []
        missing = []
        unreadable = []
        total = len(self.matched)

        for i, (node, path, meta) in enumerate(self.matched, start=1):
            # Progress, because this is the slow step when many files are
            # selected and the library sits on a synced network drive.
            self.match_lbl.config(text=f"Loading {i} of {total}:  {os.path.basename(path)}")
            try:
                self.update_idletasks()
            except Exception:
                pass

            if not self.app._ensure_file_loaded(path, node):
                unreadable.append(os.path.basename(path))
                continue
            d = self.app.cache.get(path)
            if not d:
                unreadable.append(os.path.basename(path))
                continue

            cyc_arr = d['cyc']
            hit = False
            for c in want:
                sel = (cyc_arr == c)
                if not np.any(sel):
                    continue
                E = convert_reference(d['E'][sel], ref, pH)
                I = d['I'][sel]
                # anion and cycle carried alongside: colour is keyed on the
                # anion, line style on the cycle.
                self.series.append((self._label_for(meta, path, c),
                                    meta.get('anion') or '?', c, E, I))
                hit = True
            if hit:
                used_paths.append(path)
            else:
                missing.append(os.path.basename(path))

        wanted_txt = ', '.join(str(c) for c in want)

        if not self.series:
            self.match_lbl.config(text=f"{total} file(s) matched.", fg='#555')
            if not silent:
                messagebox.showwarning(
                    "Overlay",
                    f"None of the {total} matched file(s) contain cycle(s) {wanted_txt}."
                    + (f"\n\n{len(unreadable)} file(s) could not be read." if unreadable else ""),
                    parent=self)
            self._show_files([], "plotted")
            self._redraw_empty()
            return

        status = f"{len(used_paths)} of {total} file(s) plotted, cycle(s) {wanted_txt}."
        if unreadable:
            status += f"  {len(unreadable)} unreadable."
        self.match_lbl.config(text=status, fg='#555')

        self._plotted_meas = self.meas_var.get()
        self._update_auto_title()
        self._show_files(used_paths, "plotted")
        if missing:
            self.file_list.insert('end', f"(skipped, no cycle {wanted_txt}: {len(missing)})")
        if unreadable:
            self.file_list.insert('end', f"(could not read: {len(unreadable)})")

        self.apply_format()

    # ---------------- drawing ----------------
    @staticmethod
    def _num(var, default=None):
        try:
            return float(var.get())
        except Exception:
            return default

    def _redraw_empty(self):
        self.ax.clear()
        self.ax.grid(True)
        self.ax.set_xlabel(self.xlabel_var.get())
        self.ax.set_ylabel(self.ylabel_var.get())
        self.ax.text(0.5, 0.5, "Select anion + measurement type,\nthen Find matching files",
                     ha='center', va='center', transform=self.ax.transAxes, color='#888')
        self.canvas.draw()

    def apply_format(self):
        """Redraw with the current appearance settings. Safe to call with no
        data - it just restyles the empty axes."""
        self.ax.clear()
        if not self.series:
            self._redraw_empty()
            return

        lw = self._num(self.lw_var, 1.8) or 1.8
        ms = self._num(self.ms_var, 5) or 5

        # Colour identifies the anion; line style identifies the cycle.
        anions_in_order = []
        cycles_in_order = []
        for _lab, anion, cyc, _E, _I in self.series:
            if anion not in anions_in_order:
                anions_in_order.append(anion)
            if cyc not in cycles_in_order:
                cycles_in_order.append(cyc)
        cycles_in_order.sort()

        palette = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']
        colour_of = {a: palette[i % len(palette)] for i, a in enumerate(anions_in_order)}
        style_of = {c: self.CYCLE_STYLES[i % len(self.CYCLE_STYLES)]
                    for i, c in enumerate(cycles_in_order)}

        # The anion is named once in the legend; later curves of the same
        # anion are drawn with label '_nolegend_'.
        legended = set()
        for label, anion, cyc, E, I in self.series:
            ls, mk = style_of[cyc]
            kw = dict(lw=lw, color=colour_of[anion], linestyle=ls)
            if mk:
                kw.update(marker=mk, markersize=ms,
                          markevery=max(1, len(E) // 25),
                          markerfacecolor='none')
            if anion in legended:
                self.ax.plot(E, I, label='_nolegend_', **kw)
            else:
                self.ax.plot(E, I, label=anion, **kw)
                legended.add(anion)

        fs_t = int(self._num(self.fs_title, 14) or 14)
        fs_l = int(self._num(self.fs_label, 12) or 12)
        fs_k = int(self._num(self.fs_tick, 10) or 10)
        fs_g = int(self._num(self.fs_legend, 9) or 9)

        if self.title_var.get().strip():
            self.ax.set_title(self.title_var.get(), fontsize=fs_t)
        self.ax.set_xlabel(self.xlabel_var.get(), fontsize=fs_l)
        self.ax.set_ylabel(self.ylabel_var.get(), fontsize=fs_l)
        self.ax.tick_params(axis='both', labelsize=fs_k)
        self.ax.grid(self.show_grid.get())

        xmin, xmax = self._num(self.xmin_var), self._num(self.xmax_var)
        ymin, ymax = self._num(self.ymin_var), self._num(self.ymax_var)
        if xmin is not None or xmax is not None:
            self.ax.set_xlim(left=xmin, right=xmax)
        if ymin is not None or ymax is not None:
            self.ax.set_ylim(bottom=ymin, top=ymax)

        leg1 = None
        if self.show_legend.get():
            leg1 = self.ax.legend(fontsize=fs_g, loc=self.legend_loc.get(),
                                  title="Anion")

        # A second legend keyed on line style, so the cycle numbers are
        # readable. Without it the styles are visible but unexplained.
        if self.show_cycle_legend.get() and len(cycles_in_order) > 1:
            handles = []
            for c in cycles_in_order:
                ls, mk = style_of[c]
                handles.append(Line2D([], [], color='0.25', lw=lw, linestyle=ls,
                                      marker=mk, markersize=ms,
                                      markerfacecolor='none',
                                      label=f"cycle {c}"))
            leg2 = self.ax.legend(handles=handles, fontsize=fs_g,
                                  loc=self.cycle_legend_loc.get(), title="Cycle")
            # Adding a second legend removes the first, so re-attach it.
            if leg1 is not None:
                self.ax.add_artist(leg1)

        self.fig.tight_layout()
        self.canvas.draw()

    # ---------------- automatic updates ----------------
    def _schedule_format(self, *_):
        if getattr(self, '_fmt_after', None) is not None:
            try:
                self.after_cancel(self._fmt_after)
            except Exception:
                pass
        self._fmt_after = self.after(250, self._auto_format)

    def _auto_format(self):
        self._fmt_after = None
        if self.series:
            self.apply_format()

    def _schedule_replot(self, delay=300):
        """Debounced: typing '1-5' or clicking several anions re-plots once."""
        if getattr(self, '_replot_after', None) is not None:
            try:
                self.after_cancel(self._replot_after)
            except Exception:
                pass
        self._replot_after = self.after(delay, self._auto_replot)

    def _auto_replot(self):
        self._replot_after = None
        self.find_files(silent=True)
        want = self._parse_cycles()
        if not self.matched:
            self.series = []
            self._redraw_empty()
            return
        if not want:   # empty or not readable yet (still typing)
            return
        self.plot_overlay(silent=True)

    # ---------------- auto title ----------------
    def _auto_title_text(self):
        anions = []
        for _lab, anion, _cyc, _E, _I in self.series:
            if anion and anion != '?' and anion not in anions:
                anions.append(anion)
        meas = self._plotted_meas or self.meas_var.get()
        desc = MEAS_DESCRIPTIONS.get(meas, meas)
        return " - ".join(p for p in (", ".join(anions), desc) if p)

    def _update_auto_title(self):
        """Put the auto title into the Title box, unless the user typed their own."""
        auto = self._auto_title_text()
        cur = self.title_var.get().strip()
        if cur == '' or cur == self._last_auto_title:
            self.title_var.set(auto)
        self._last_auto_title = auto

    def _reset_auto_title(self):
        self._last_auto_title = self._auto_title_text()
        self.title_var.set(self._last_auto_title)

    # ---------------- saved Y range per measurement type ----------------
    def _cycles_key(self):
        want = self._parse_cycles()
        if want is None:
            return None
        return ','.join(str(c) for c in want) if want else 'all'

    def _store_y_default(self, meas, info):
        st = overlay_load_settings()
        st.setdefault('y_defaults', {}).setdefault(meas, {})[info['cycles']] = info
        overlay_save_settings(st)

    def _save_y_default_manual(self):
        meas = self.meas_var.get()
        key = self._cycles_key()
        ymin, ymax = self._num(self.ymin_var), self._num(self.ymax_var)
        if not meas or key is None or ymin is None or ymax is None:
            messagebox.showwarning("Overlay", "Choose a measurement type and fill both Y min and Y max.",
                                   parent=self)
            return
        self._store_y_default(meas, {'cycles': key, 'ymin': ymin, 'ymax': ymax, 'source': 'manual',
                                     'saved': time.strftime('%Y-%m-%d %H:%M')})
        self.yscan_lbl.config(text=f"Saved as default for {meas} ({key if key != 'all' else 'all cycles'}): "
                                   f"{ymin:g} to {ymax:g}", fg='#1a5c1a')

    def _load_y_default(self):
        """When the measurement type or cycles change: fill Y min / max from the saved
        default for that type (same cycles if saved, otherwise the latest one saved for
        the type). No saved default -> automatic Y."""
        meas = self.meas_var.get()
        key = self._cycles_key()
        if not meas or key is None:
            return
        saved = overlay_load_settings().get('y_defaults', {}).get(meas, {})
        info = saved.get(key)
        if info is None and saved:
            info = max(saved.values(), key=lambda d: d.get('saved', ''))
        if info is None:
            self.ymin_var.set('')
            self.ymax_var.set('')
            self.yscan_lbl.config(text=f"No saved Y range for {meas} yet - press Scan, "
                                       f"or type Y min / max under Appearance.", fg='#555')
            return
        self.ymin_var.set(f"{info['ymin']:g}")
        self.ymax_var.set(f"{info['ymax']:g}")
        cyc_txt = 'all cycles' if info['cycles'] == 'all' else f"cycle(s) {info['cycles']}"
        txt = f"Saved Y range for {meas}, {cyc_txt}:  {info['ymin']:g} to {info['ymax']:g}"
        if info.get('source') == 'scan':
            txt += (f"\nMax {info['max_val']:.4g} mA in {info['max_file']} (cycle {info['max_cycle']})"
                    f"\nMin {info['min_val']:.4g} mA in {info['min_file']} (cycle {info['min_cycle']})")
        else:
            txt += "\n(set manually)"
        if info['cycles'] != key:
            txt += "\n(saved for other cycles - scan again for these)"
        txt += f"\nsaved {info.get('saved', '')}"
        self.yscan_lbl.config(text=txt, fg='#1a5c1a' if info['cycles'] == key else '#b36b00')

    # ---------------- global Y range ----------------
    def _stop_y_scan(self):
        self._y_scan_stop = True

    def _clear_y_limits(self):
        self.ymin_var.set('')
        self.ymax_var.set('')
        self.yscan_lbl.config(text="Y axis back to automatic.", fg='#555')
        if self.series:
            self.apply_format()

    def scan_global_y(self):
        meas = self.meas_var.get()
        if not meas:
            messagebox.showwarning("Overlay", "Select a measurement type.", parent=self)
            return
        want = self._parse_cycles()
        if want is None:
            messagebox.showwarning("Overlay", "Could not read the cycle numbers.", parent=self)
            return
        files = [(n, p) for n, p, m in self.candidates if m['meas'] == meas]
        if not files:
            messagebox.showwarning("Overlay", f"No files of type {meas}.", parent=self)
            return

        disk = _yr_load_cache()
        self._y_scan_stop = False
        best_hi = (-np.inf, None, None)   # (value, path, cycle)
        best_lo = (np.inf, None, None)
        used = read_new = 0
        unreadable = []
        total = len(files)

        for k, (node, path) in enumerate(files, start=1):
            if self._y_scan_stop:
                break
            name = os.path.basename(path)
            try:
                st = os.stat(path)
                sig = [st.st_size, int(st.st_mtime)]
            except OSError:
                unreadable.append(name)
                continue

            ent = disk.get(path)
            if ent and ent.get('sig') == sig:
                ext = {int(c): tuple(v) for c, v in ent['cycles'].items()}
            else:
                self.yscan_lbl.config(text=f"Reading {k} of {total}:  {name}", fg='#555')
                try:
                    self.update()          # keep the window (and Stop) responsive
                except Exception:
                    pass
                d = self.app.cache.get(path)
                try:
                    if d is not None:
                        I, cyc = d['I'], d['cyc']
                    else:
                        # read directly WITHOUT keeping the data in memory
                        _h, _d, _E, I, _t, cyc, _m = read_mpt(path)
                    ext = yr_cycle_extrema(I, cyc)
                except Exception:
                    unreadable.append(name)
                    continue
                disk[path] = {'sig': sig, 'cycles': {str(c): list(v) for c, v in ext.items()}}
                read_new += 1

            cycles = want if want else list(ext.keys())
            hit = False
            for c in cycles:
                if c not in ext:
                    continue
                lo, hi = ext[c]
                hit = True
                if hi > best_hi[0]:
                    best_hi = (hi, path, c)
                if lo < best_lo[0]:
                    best_lo = (lo, path, c)
            used += hit

        if read_new:
            _yr_save_cache(disk)

        stopped = self._y_scan_stop
        cyc_txt = ('cycle(s) ' + ', '.join(str(c) for c in want)) if want else 'all cycles'
        if best_hi[1] is None:
            self.yscan_lbl.config(
                text=f"No data found for {meas}, {cyc_txt}."
                     + (" (stopped)" if stopped else ""), fg='#b00020')
            return

        ymin, ymax = yr_nice_limits(best_lo[0], best_hi[0])
        self.ymin_var.set(f"{ymin:g}")
        self.ymax_var.set(f"{ymax:g}")
        if not stopped:   # only a complete scan becomes the default for this type
            self._store_y_default(meas, {
                'cycles': ','.join(str(c) for c in want) if want else 'all',
                'ymin': ymin, 'ymax': ymax, 'source': 'scan',
                'max_val': best_hi[0], 'max_file': os.path.basename(best_hi[1]), 'max_cycle': best_hi[2],
                'min_val': best_lo[0], 'min_file': os.path.basename(best_lo[1]), 'min_cycle': best_lo[2],
                'saved': time.strftime('%Y-%m-%d %H:%M')})

        txt = (f"{meas}, {cyc_txt}: {used} of {total} file(s)"
               + (" - STOPPED, partial result" if stopped else "") + "\n"
               f"Max: {best_hi[0]:.4g} mA  in  {os.path.basename(best_hi[1])}  (cycle {best_hi[2]})\n"
               f"Min: {best_lo[0]:.4g} mA  in  {os.path.basename(best_lo[1])}  (cycle {best_lo[2]})\n"
               f"Y axis for all graphs:  {ymin:g}  to  {ymax:g}")
        if unreadable:
            txt += f"\n({len(unreadable)} file(s) could not be read)"
        txt += "\n(saved as default for this type)" if not stopped else "\n(not saved - scan was stopped)"
        self.yscan_lbl.config(text=txt, fg='#b36b00' if stopped else '#1a5c1a')
        if self.series:
            self.apply_format()

    # ---------------- appendix ----------------
    def _appendix_default_title(self):
        """The overlay Title if one was typed; otherwise anion(s) + measurement wording,
        e.g. 'Carbonate, Nitrate - activation with urea'."""
        typed = self.title_var.get().strip()
        if typed:
            return typed
        return self._auto_title_text()

    def send_to_appendix(self):
        if not self.series:
            messagebox.showwarning("Overlay", "Plot an overlay first.", parent=self)
            return
        self.apply_format()  # make sure the copy matches the current appearance settings
        if self.appendix.add_figure(self.fig, self._appendix_default_title()):
            self.match_lbl.config(
                text=f"Sent to Appendix page ({len(self.appendix.panels)} graph(s) so far).",
                fg='#555')

    def _update_appendix_tab_text(self):
        n = len(self.app.appendix_panels) if hasattr(self.app, 'appendix_panels') else 0
        try:
            self.right_nb.tab(self.appendix_tab,
                              text=f"Appendix page ({n})" if n else "Appendix page")
        except Exception:
            pass

    # ---------------- output ----------------
    def save_figure(self):
        if not self.series:
            messagebox.showwarning("Overlay", "Nothing to save yet.", parent=self)
            return
        fn = filedialog.asksaveasfilename(
            parent=self, title="Save overlay figure",
            defaultextension='.png',
            initialdir=(EXPORT_FOLDER if os.path.isdir(EXPORT_FOLDER) else None),
            filetypes=[("PNG", "*.png"), ("SVG (vector)", "*.svg"),
                       ("PDF (vector)", "*.pdf"), ("TIFF", "*.tif")])
        if not fn:
            return
        try:
            dpi = int(float(self.dpi_var.get()))
        except Exception:
            dpi = 300
        try:
            self.fig.savefig(fn, dpi=dpi, bbox_inches='tight')
        except Exception as e:
            messagebox.showerror("Overlay", f"Could not save:\n{e}", parent=self)
            return
        messagebox.showinfo("Overlay", f"Saved:\n{fn}", parent=self)

    def copy_data(self):
        """Put the plotted curves on the clipboard as TSV, for replotting
        elsewhere (Origin, Excel)."""
        if not self.series:
            messagebox.showwarning("Overlay", "Nothing to copy yet.", parent=self)
            return
        cols = []
        for label, _anion, _cyc, E, I in self.series:
            cols.append((f"E_{label}", E))
            cols.append((f"I_{label}", I))
        n = max(len(v) for _h, v in cols)
        lines = ["\t".join(h for h, _v in cols)]
        for i in range(n):
            lines.append("\t".join(
                (f"{v[i]:.6g}" if i < len(v) else "") for _h, v in cols))
        try:
            self.clipboard_clear()
            self.clipboard_append("\n".join(lines))
            messagebox.showinfo("Overlay",
                                f"{len(self.series)} curve(s) copied as TSV.", parent=self)
        except Exception as e:
            messagebox.showwarning("Overlay", f"Could not copy: {e}", parent=self)


# =================== DEFAULT FOLDERS ===================
# After the tutorial the user picks the data folder once (it can be a local
# folder or a synced cloud folder, e.g. Google Drive for desktop); the choice
# is saved in the settings.


def find_shared_drives_root():
    """'X:\\Shared drives' of Google Drive for desktop, if one exists."""
    for letter in 'GHIJKLMNOPQRSTUVWXYZDEF':
        for name in ('Shared drives', 'Shared Drives', 'אחסון שיתופי', 'תיקיות אחסון שיתופי'):
            p = f"{letter}:\\{name}"
            if os.path.isdir(p):
                return p
    return None


FOLDER_TEXT = {
    'en': ("Choose your data folder",
           "Choose the folder that holds your .mpt files, and then the folder for "
           "the exported CSV files.\n\n"
           "A synced cloud folder works too, for example a shared drive in "
           "Google Drive for desktop:\nG:\\Shared drives\\<drive>\\...\n\n"
           "You can change both later under File > Set default folders...",
           "Choose folders...", "Later"),
    'he': ("בחירת תיקיית הנתונים",
           "בחר את התיקייה שבה נמצאים קובצי ה-mpt, ואחר כך את התיקייה לקובצי ה-CSV "
           "המיוצאים.\n\n"
           "אפשר גם תיקייה מסונכרנת מהענן, למשל דרייב שיתופי דרך Google Drive for desktop:\n"
           "G:\\Shared drives\\<drive>\\...\n\n"
           "אפשר לשנות את שתיהן בהמשך דרך File > Set default folders...",
           "בחירת תיקיות...", "אחר כך"),
}


class FolderSetupDialog(tk.Toplevel):
    """Explains where the data folder is set. choose=True if the user wants to pick folders now."""

    def __init__(self, parent, lang='en'):
        super().__init__(parent)
        self.choose = False
        title, body, b_choose, b_later = FOLDER_TEXT.get(lang, FOLDER_TEXT['en'])
        rtl = lang == 'he'
        self.title(title)
        self.resizable(False, False)
        tk.Label(self, text=title, font=('Arial', 12, 'bold')).pack(
            anchor=('e' if rtl else 'w'), padx=14, pady=(12, 4))
        tk.Label(self, text=body, font=('Arial', 10),
                 wraplength=440, justify=('right' if rtl else 'left')).pack(padx=14)
        row = tk.Frame(self)
        row.pack(fill='x', padx=10, pady=12)
        tk.Button(row, text=b_choose, font=('Arial', 10), command=self._choose).pack(
            side=('left' if rtl else 'right'), padx=3)
        tk.Button(row, text=b_later, relief='flat', fg='#888', font=('Arial', 9),
                  command=self.destroy).pack(side=('right' if rtl else 'left'))
        self.transient(parent)
        self.update_idletasks()
        try:
            x = parent.winfo_rootx() + (parent.winfo_width() - self.winfo_width()) // 2
            y = parent.winfo_rooty() + (parent.winfo_height() - self.winfo_height()) // 3
            self.geometry(f"+{max(0, x)}+{max(0, y)}")
        except Exception:
            pass
        self.grab_set()
        parent.wait_window(self)

    def _choose(self):
        self.choose = True
        self.destroy()

def _set_folders(mpt=None, export=None):
    """Replace the default .mpt / export folders (from saved settings, or the tutorial)."""
    global DEFAULT_FOLDER, EXPORT_FOLDER
    if mpt:
        DEFAULT_FOLDER = mpt
    if export:
        EXPORT_FOLDER = export


# =================== TUTORIAL (runs on first start) ===================
# A guided tour on two synthetic sample files. Every action step waits until
# the user has actually done it before "Next" unlocks. The tour starts by
# itself on every launch until it has been finished once; it can be replayed
# from Help > Tutorial.

TUTORIAL_DIR = os.path.join(os.path.expanduser('~'), 'CV_Analyzer_tutorial')
TUTORIAL_EXPORT_DIR = os.path.join(TUTORIAL_DIR, 'Processed')
# (file name, peak shift in V, peak scale) - names follow the lab naming scheme,
# so the export dialog and the overlay window can read them.
TUTORIAL_SAMPLES = [
    ('UAL_A_10001_Nitrate_tutorial.mpt', 0.000, 1.00),
    ('UAL_A_10002_Carbonate_tutorial.mpt', 0.015, 0.80),
]


def make_tutorial_samples(folder=TUTORIAL_DIR):
    """Write the sample .mpt files (EC-Lab ASCII layout) if they are missing.

    Synthetic Ni(OH)2-like CV vs RHE, 10 mV/s, 3 cycles: an anodic peak near
    1.42 V, a cathodic peak near 1.33 V, and an OER tail above 1.5 V.
    """
    os.makedirs(folder, exist_ok=True)
    rng = np.random.default_rng(7)
    step, rate = 0.002, 0.010                     # V per point, V/s
    up = np.arange(1.00, 1.60 + 1e-9, step)
    sweep = np.concatenate([up, up[-2::-1]])
    for name, shift, scale in TUTORIAL_SAMPLES:
        path = os.path.join(folder, name)
        if os.path.exists(path):
            continue
        lines = ["EC-Lab ASCII FILE",
                 "Nb header lines : 4",
                 "Synthetic sample data written by CV Analyzer for the tutorial",
                 "mode\ttime/s\tEwe/V\t<I>/mA\tcycle number"]
        t = 0.0
        for cyc in (1, 2, 3):
            grow = 1.0 + 0.08 * (cyc - 1)
            for k, E in enumerate(sweep):
                if k < len(up):   # forward sweep
                    I = 0.03 + 1.2 * scale * grow * np.exp(-((E - (1.42 + shift)) / 0.025) ** 2)
                else:             # backward sweep
                    I = -0.03 - 0.9 * scale * grow * np.exp(-((E - (1.33 + shift)) / 0.022) ** 2)
                I += 0.01 * np.exp((E - 1.50) / 0.03) + rng.normal(0, 0.004)
                lines.append(f"2\t{t:.3f}\t{E:.5f}\t{I:.6f}\t{cyc}")
                t += step / rate
        with open(path, 'w', encoding='latin1') as fh:
            fh.write("\n".join(lines) + "\n")


TUT_UI = {
    'en': {
        'title': "CV Analyzer - Tutorial",
        'step': "Step {i} of {n}",
        'next': "Next", 'back': "Back", 'finish': "Finish", 'exit': "Exit tutorial",
        'waiting': "Waiting for you to do this step...",
        'done': "✓ Done - press Next",
        'exit_q': "Exit the tutorial?\n\nIt will start again the next time you open "
                  "the program, until you finish it once.",
    },
    'he': {
        'title': "CV Analyzer - הדרכה",
        'step': "שלב {i} מתוך {n}",
        'next': "הבא", 'back': "הקודם", 'finish': "סיום", 'exit': "יציאה מההדרכה",
        'waiting': "ממתין שתבצע את השלב...",
        'done': "✓ בוצע - לחץ 'הבא'",
        'exit_q': "לצאת מההדרכה?\n\nהיא תופעל שוב בפעם הבאה שתפתח את התוכנה, "
                  "עד שתסיים אותה פעם אחת.",
    },
}

# key -> {lang: (title, body)}.  {folder} / {export} are filled in at run time.
TUT_STEPS = {
    'welcome': {
        'en': ("Welcome to CV Analyzer",
               "This short tour runs a complete analysis on two sample files that were "
               "just created for you in:\n{folder}\n\n"
               "Each step asks you to do one thing in the program. The Next button "
               "unlocks once you have done it.\n\n"
               "The file list was cleared for the tour - none of your files were changed."),
        'he': ("ברוכים הבאים ל-CV Analyzer",
               "הסיור הקצר הזה מריץ ניתוח מלא על שני קבצים לדוגמה שנוצרו עבורך בתיקייה:\n"
               "{folder}\n\n"
               "בכל שלב תתבקש לבצע פעולה אחת בתוכנה. כפתור 'הבא' ייפתח אחרי שתבצע אותה.\n\n"
               "רשימת הקבצים נוקתה לצורך הסיור - אף קובץ שלך לא שונה."),
    },
    'open_folder': {
        'en': ("Open a folder",
               "Click 'Open Folder' (top left) and choose the tutorial folder. "
               "The dialog already opens there - just confirm.\n\n"
               "Every .mpt file in the folder and its subfolders is listed in the tree, "
               "along with CSV / TXT files exported by MPR Converter."),
        'he': ("פתיחת תיקייה",
               "לחץ על 'Open Folder' (למעלה משמאל) ובחר את תיקיית ההדרכה. "
               "החלון כבר נפתח בה - רק אשר.\n\n"
               "כל קובצי ה-mpt בתיקייה ובתתי-התיקיות יופיעו בעץ, "
               "וגם קובצי CSV ו-TXT שיוצאו מ-MPR Converter."),
    },
    'tree': {
        'en': ("Select cycles",
               "Tick the ✔ column next to a file name. This selects all of its cycles "
               "and plots them.\n\n"
               "Tip: the ▸ arrow shows each cycle with its scan rate, and single cycles "
               "can be ticked. Clicking a cycle's 'Color' cell changes its color."),
        'he': ("בחירת מחזורים",
               "סמן את עמודת ה-✔ ליד שם של קובץ. כך נבחרים כל המחזורים שלו ומוצגים בגרף.\n\n"
               "טיפ: החץ ▸ מציג כל מחזור עם קצב הסריקה שלו, ואפשר לסמן מחזורים בודדים. "
               "לחיצה על תא ה-Color של מחזור משנה את הצבע שלו."),
    },
    'reference': {
        'en': ("Reference electrode and pH",
               "All potentials are converted to RHE using the reference electrode and "
               "the pH set here.\n\n"
               "The sample files are already vs RHE, so leave 'RHE'. With your own data, "
               "choose the reference electrode and the electrolyte pH before analyzing."),
        'he': ("אלקטרודת ייחוס ו-pH",
               "כל הפוטנציאלים מומרים ל-RHE לפי אלקטרודת הייחוס וה-pH שנבחרו כאן.\n\n"
               "קובצי הדוגמה כבר מול RHE, אז השאר 'RHE'. בנתונים שלך, בחר את אלקטרודת "
               "הייחוס ואת ה-pH של האלקטרוליט לפני הניתוח."),
    },
    'plot_q': {
        'en': ("Plot modes",
               "'Plot Mode' changes what is drawn. Choose 'QvsE' (charge vs potential) - "
               "the plot updates by itself, as it does for every setting in the toolbar."),
        'he': ("מצבי תצוגה",
               "'Plot Mode' קובע מה מוצג בגרף. בחר 'QvsE' (מטען מול פוטנציאל) - "
               "הגרף מתעדכן לבד, כמו בכל הגדרה בסרגל הכלים."),
    },
    'plot_i': {
        'en': ("Back to current",
               "Now choose 'IvsE' again.\n\n"
               "The other modes: IvsT - current vs time, EvsT - potential vs time."),
        'he': ("חזרה לזרם",
               "עכשיו בחר שוב 'IvsE'.\n\n"
               "המצבים האחרים: IvsT - זרם מול זמן, EvsT - פוטנציאל מול זמן."),
    },
    'run': {
        'en': ("Run the analysis",
               "With the cycles still ticked, click 'Run Analysis on Checked Cycles' "
               "in the Analysis panel on the right."),
        'he': ("הרצת הניתוח",
               "כשהמחזורים עדיין מסומנים, לחץ על 'Run Analysis on Checked Cycles' "
               "בחלונית Analysis שמימין."),
    },
    'report': {
        'en': ("Reading the results",
               "Top: one row per cycle - anodic and cathodic peaks, charge, OER-corrected "
               "values and the currents at fixed potentials. Scroll sideways for all columns.\n\n"
               "Bottom: the report - key parameters first, then anything missing, then the "
               "full statistics. 'Copy report' copies it to the clipboard."),
        'he': ("קריאת התוצאות",
               "למעלה: שורה לכל מחזור - פיקים אנודיים וקתודיים, מטען, ערכים אחרי תיקון OER "
               "והזרמים בפוטנציאלים קבועים. גלול הצידה כדי לראות את כל העמודות.\n\n"
               "למטה: הדוח - קודם הפרמטרים המרכזיים, אחר כך מה שחסר, ובסוף הסטטיסטיקה "
               "המלאה. 'Copy report' מעתיק אותו ללוח."),
    },
    'export': {
        'en': ("Export a CSV",
               "Click 'Export Analysis CSV'. The name fields are filled in from the .mpt "
               "file name - check them and press 'Save file'. Type your name in "
               "'Laborant' - a new name is kept for next time ('Edit list...' manages the list).\n\n"
               "During the tour, files are saved to:\n{export}"),
        'he': ("ייצוא ל-CSV",
               "לחץ על 'Export Analysis CSV'. שדות השם ממולאים אוטומטית משם קובץ ה-mpt - "
               "בדוק אותם ולחץ על 'Save file'. הקלד את שמך ב-'Laborant' - שם חדש נשמר "
               "לפעם הבאה ('Edit list...' לעריכת הרשימה).\n\n"
               "במהלך הסיור הקבצים נשמרים ב:\n{export}"),
    },
    'batch': {
        'en': ("Process a whole folder",
               "The second sample file has no export yet. Click 'Find unanalyzed...' above "
               "the tree and confirm. For each file you see its report - press "
               "'Approve & Save'."),
        'he': ("עיבוד תיקייה שלמה",
               "לקובץ הדוגמה השני עוד אין ייצוא. לחץ על 'Find unanalyzed...' מעל העץ ואשר. "
               "לכל קובץ תראה את הדוח שלו - לחץ על 'Approve & Save'."),
    },
    'overlay': {
        'en': ("Overlay window",
               "Click 'Overlay...' above the tree. It compares the same cycle across "
               "many files."),
        'he': ("חלון ההשוואה",
               "לחץ על 'Overlay...' מעל העץ. החלון משווה את אותו מחזור בין קבצים רבים."),
    },
    'ov_select': {
        'en': ("Pick what to compare",
               "Select both anions in the list (Ctrl+click) and type 2 in "
               "'Cycle number(s)'. The plot updates by itself.\n\n"
               "Files are matched by name, so the anion and the measurement type must "
               "appear in the .mpt file name."),
        'he': ("בחירת מה להשוות",
               "בחר את שני האניונים ברשימה (Ctrl+קליק) והקלד 2 ב-'Cycle number(s)'. "
               "הגרף מתעדכן לבד.\n\n"
               "ההתאמה נעשית לפי שם הקובץ, ולכן שם קובץ ה-mpt חייב לכלול את האניון "
               "ואת סוג המדידה."),
    },
    'appendix': {
        'en': ("Build an appendix page",
               "Click '→ Appendix page' under the plot. The graph is added to a printable "
               "page (tab 'Appendix page'), where graphs can be arranged and the page "
               "saved as PNG or PDF."),
        'he': ("בניית דף נספח",
               "לחץ על '→ Appendix page' מתחת לגרף. הגרף נוסף לדף להדפסה (לשונית "
               "'Appendix page'), שבו אפשר לסדר גרפים ולשמור את הדף כ-PNG או PDF."),
    },
    'done': {
        'en': ("You're ready",
               "The usual workflow:\n"
               "Open Folder > tick cycles > Run Analysis > Export Analysis CSV.\n\n"
               "For many files use 'Find unanalyzed...'.\n"
               "Rarely used options (OER subtraction, onset lines, Sanity Fix) are in the "
               "⋮ menu at the right end of the toolbar.\n"
               "To choose your own folders: File > Set default folders...\n"
               "To replay this tour: Help > Tutorial."),
        'he': ("סיימת!",
               "סדר העבודה הרגיל:\n"
               "Open Folder, סימון מחזורים, Run Analysis, ואז Export Analysis CSV.\n\n"
               "להרבה קבצים השתמש ב-'Find unanalyzed...'.\n"
               "אפשרויות נדירות (חיסור OER, קווי onset, Sanity Fix) נמצאות בתפריט ⋮ "
               "בקצה הימני של סרגל הכלים.\n"
               "לבחירת התיקיות שלך: File ואז Set default folders...\n"
               "להפעלת הסיור שוב: Help ואז Tutorial."),
    },
}


class TutorialLanguageDialog(tk.Toplevel):
    """Small modal window: pick English or Hebrew. result = 'en' | 'he' | None."""

    def __init__(self, parent, first_run):
        super().__init__(parent)
        self.result = None
        self.title("CV Analyzer - Tutorial")
        self.resizable(False, False)
        txt = ("Welcome! A short guided tour will show you how to use the program.\n"
               "Choose the language of the tour:" if first_run else
               "Choose the language of the tour:")
        tk.Label(self, text=txt, justify='left', padx=16, pady=10).pack()
        tk.Label(self, text="בחר את שפת ההדרכה:", font=('Arial', 10)).pack(padx=16)
        row = tk.Frame(self)
        row.pack(pady=12)
        tk.Button(row, text="English", width=12,
                  command=lambda: self._pick('en')).pack(side='left', padx=6)
        tk.Button(row, text="עברית", width=12, font=('Arial', 10),
                  command=lambda: self._pick('he')).pack(side='left', padx=6)
        self.protocol("WM_DELETE_WINDOW", self.destroy)
        self.transient(parent)
        self.update_idletasks()
        try:
            x = parent.winfo_rootx() + (parent.winfo_width() - self.winfo_width()) // 2
            y = parent.winfo_rooty() + (parent.winfo_height() - self.winfo_height()) // 3
            self.geometry(f"+{max(0, x)}+{max(0, y)}")
        except Exception:
            pass
        self.grab_set()
        parent.wait_window(self)

    def _pick(self, lang):
        self.result = lang
        self.destroy()


class Tutorial:
    """Guided tour: a floating instruction window plus a coloured frame around
    the control the current step is about. Action steps are checked every
    0.3 s; Next unlocks only when the action has really been done."""

    HL_WAIT = '#ff7a00'   # frame colour while waiting
    HL_DONE = '#1a9e3a'   # frame colour once the step is done

    def __init__(self, app, lang):
        self.app = app
        self.root = app.root
        self.lang = lang if lang in TUT_UI else 'en'
        self.ui = TUT_UI[self.lang]
        self.i = 0
        self.done = set()
        self._after = None
        self._placed_for = None
        self._hidden = False
        self._panels0 = len(getattr(app, 'appendix_panels', []))

        # sample data + a private export folder for the tour
        make_tutorial_samples(TUTORIAL_DIR)
        os.makedirs(TUTORIAL_EXPORT_DIR, exist_ok=True)
        for fn in os.listdir(TUTORIAL_EXPORT_DIR):     # start the export steps clean
            if fn.lower().endswith('.csv'):
                try:
                    os.remove(os.path.join(TUTORIAL_EXPORT_DIR, fn))
                except Exception:
                    pass
        self._old_export = EXPORT_FOLDER
        _set_folders(export=TUTORIAL_EXPORT_DIR)
        app.tutorial_folder = TUTORIAL_DIR

        self.steps = self._make_steps()
        self._build_window()
        self._make_edges()
        self._show(0)
        self._tick()

    # ---------- steps ----------
    def _ow(self):
        ow = getattr(self.app, '_overlay_win', None)
        try:
            return ow if (ow is not None and ow.winfo_exists()) else None
        except Exception:
            return None

    def _n_csv(self):
        try:
            return sum(1 for f in os.listdir(TUTORIAL_EXPORT_DIR) if f.lower().endswith('.csv'))
        except Exception:
            return 0

    def _plotted_modes(self):
        return {v.get('mode') for v in getattr(self.app, 'last_debug', {}).values()}

    def _make_steps(self):
        a = self.app
        norm = lambda p: os.path.normcase(os.path.normpath(p))
        return [
            dict(key='welcome', target=lambda: None, check=None),
            dict(key='open_folder', target=lambda: a.btn_open_folder,
                 check=lambda: any(norm(r) == norm(TUTORIAL_DIR) for r in a.open_roots)),
            dict(key='tree', target=lambda: a.tree,
                 check=lambda: bool(a._all_checked_cycles())),
            dict(key='reference', target=lambda: a.ref_combo, check=None),
            dict(key='plot_q', target=lambda: a.plot_mode_combo,
                 check=lambda: 'QvsE' in self._plotted_modes()),
            dict(key='plot_i', target=lambda: a.plot_mode_combo,
                 check=lambda: self._plotted_modes() == {'IvsE'}),
            dict(key='run', target=lambda: a.btn_run_analysis,
                 check=lambda: bool(a.analysis_rows)),
            dict(key='report', target=lambda: a.analysis_text, check=None),
            dict(key='export', target=lambda: a.btn_export,
                 check=lambda: self._n_csv() >= 1),
            dict(key='batch', target=lambda: a.btn_find_unanalyzed,
                 check=lambda: self._n_csv() >= 2),
            dict(key='overlay', target=lambda: a.btn_overlay,
                 check=lambda: self._ow() is not None),
            dict(key='ov_select',
                 target=lambda: self._ow().anion_list if self._ow() else None,
                 check=lambda: self._ow() is not None
                 and len({s[1] for s in self._ow().series}) >= 2),
            dict(key='appendix',
                 target=lambda: self._ow().btn_appendix if self._ow() else None,
                 check=lambda: len(getattr(a, 'appendix_panels', [])) > self._panels0),
            dict(key='done', target=lambda: None, check=None),
        ]

    # ---------- windows ----------
    def _build_window(self):
        rtl = self.lang == 'he'
        anchor, justify = ('e', 'right') if rtl else ('w', 'left')
        b = tk.Toplevel(self.root)
        b.title(self.ui['title'])
        b.resizable(False, False)
        b.attributes('-topmost', True)
        b.protocol("WM_DELETE_WINDOW", self._ask_exit)
        self.bubble = b

        self.lbl_count = tk.Label(b, fg='#777', font=('Arial', 9))
        self.lbl_count.pack(anchor=anchor, padx=14, pady=(10, 0))
        self.lbl_title = tk.Label(b, font=('Arial', 12, 'bold'), wraplength=380, justify=justify)
        self.lbl_title.pack(anchor=anchor, padx=14, pady=(2, 6))
        self.lbl_body = tk.Label(b, font=('Arial', 10), wraplength=380, justify=justify)
        self.lbl_body.pack(anchor=anchor, padx=14)
        self.lbl_status = tk.Label(b, font=('Arial', 10, 'bold'))
        self.lbl_status.pack(anchor=anchor, padx=14, pady=(10, 4))

        row = tk.Frame(b)
        row.pack(fill='x', padx=10, pady=(4, 10))
        self.btn_next = tk.Button(row, width=10, font=('Arial', 10), command=self._next)
        self.btn_back = tk.Button(row, width=8, font=('Arial', 10), command=self._back)
        self.btn_exit = tk.Button(row, text=self.ui['exit'], relief='flat', fg='#888',
                                  font=('Arial', 9), command=self._ask_exit)
        far, near = ('left', 'right') if rtl else ('right', 'left')
        self.btn_next.pack(side=far, padx=3)
        self.btn_back.pack(side=far, padx=3)
        self.btn_exit.pack(side=near)
        self.btn_back.config(text=self.ui['back'])

    def _make_edges(self):
        """Four thin borderless windows that frame the target control."""
        self.edges = []
        for _ in range(4):
            e = tk.Toplevel(self.root)
            e.overrideredirect(True)
            e.attributes('-topmost', True)
            e.configure(bg=self.HL_WAIT)
            e.withdraw()
            self.edges.append(e)

    def _hide_edges(self):
        for e in self.edges:
            try:
                e.withdraw()
            except Exception:
                pass

    def _frame(self, w, color):
        try:
            ok = w is not None and w.winfo_exists() and w.winfo_viewable()
        except Exception:
            ok = False
        if not ok:
            self._hide_edges()
            return
        x, y = w.winfo_rootx(), w.winfo_rooty()
        ww, hh = w.winfo_width(), w.winfo_height()
        t, pad = 3, 2
        X, Y = x - pad - t, y - pad - t
        W, H = ww + 2 * (pad + t), hh + 2 * (pad + t)
        for e, (gw, gh, gx, gy) in zip(self.edges, [(W, t, X, Y), (W, t, X, Y + H - t),
                                                      (t, H, X, Y), (t, H, X + W - t, Y)]):
            e.configure(bg=color)
            e.geometry(f"{gw}x{gh}+{gx}+{gy}")
            if e.state() == 'withdrawn':
                e.deiconify()

    def _place(self, w):
        """Put the instruction window next to the target (once per step, so the
        user can move it freely afterwards)."""
        b = self.bubble
        b.update_idletasks()
        bw, bh = b.winfo_reqwidth(), b.winfo_reqheight()
        sw, sh = b.winfo_screenwidth(), b.winfo_screenheight()
        try:
            ok = w is not None and w.winfo_viewable()
        except Exception:
            ok = False
        if ok:
            tx, ty, tw, th = w.winfo_rootx(), w.winfo_rooty(), w.winfo_width(), w.winfo_height()
            x = tx + tw + 20
            if x + bw > sw - 10:
                x = tx - bw - 20
            if x < 10:            # no room on either side -> below the target
                x = min(max(10, tx), sw - bw - 10)
                y = ty + th + 20
            else:
                y = ty
        else:
            r = self.root
            x = r.winfo_rootx() + (r.winfo_width() - bw) // 2
            y = r.winfo_rooty() + (r.winfo_height() - bh) // 3
        x = min(max(10, x), sw - bw - 10)
        y = min(max(10, y), sh - bh - 60)
        b.geometry(f"+{int(x)}+{int(y)}")

    # ---------- flow ----------
    def _show(self, i):
        self.i = i
        st = self.steps[i]
        title, body = TUT_STEPS[st['key']][self.lang]
        body = body.format(folder=TUTORIAL_DIR, export=TUTORIAL_EXPORT_DIR)
        self.lbl_count.config(text=self.ui['step'].format(i=i + 1, n=len(self.steps)))
        self.lbl_title.config(text=title)
        self.lbl_body.config(text=body)
        last = i == len(self.steps) - 1
        self.btn_next.config(text=self.ui['finish'] if last else self.ui['next'])
        self.btn_back.config(state=('disabled' if i == 0 else 'normal'))
        self._placed_for = None
        self._refresh()

    def _refresh(self):
        st = self.steps[self.i]
        if st['check'] is None:
            self.lbl_status.config(text='')
            self.btn_next.config(state='normal')
        elif st['key'] in self.done:
            self.lbl_status.config(text=self.ui['done'], fg=self.HL_DONE)
            self.btn_next.config(state='normal')
        else:
            self.lbl_status.config(text=self.ui['waiting'], fg=self.HL_WAIT)
            self.btn_next.config(state='disabled')

    def _tick(self):
        self._after = None
        try:
            if not self.bubble.winfo_exists():
                return
        except Exception:
            return
        st = self.steps[self.i]
        if st['check'] is not None and st['key'] not in self.done:
            try:
                if st['check']():
                    self.done.add(st['key'])
            except Exception:
                pass
        self._refresh()

        # Step aside while a dialog of the program (e.g. the export-name window)
        # holds the input, or while the main window is minimised.
        try:
            busy = self.root.grab_current() is not None or self.root.state() == 'iconic'
        except Exception:
            busy = False
        if busy:
            if not self._hidden:
                self.bubble.withdraw()
                self._hide_edges()
                self._hidden = True
        else:
            if self._hidden:
                self.bubble.deiconify()
                self._hidden = False
            try:
                target = st['target']()
            except Exception:
                target = None
            done = st['check'] is None or st['key'] in self.done
            self._frame(target, self.HL_DONE if done else self.HL_WAIT)
            if self._placed_for != self.i:
                self._place(target)
                self._placed_for = self.i
        self._after = self.root.after(300, self._tick)

    def _next(self):
        if self.i == len(self.steps) - 1:
            st = overlay_load_settings()
            st['tutorial_done'] = True
            st['tutorial_lang'] = self.lang
            overlay_save_settings(st)
            app, lang = self.app, self.lang
            self._close()
            app.root.after(300, lambda: app.offer_folder_setup(lang))
            return
        self._show(self.i + 1)

    def _back(self):
        if self.i > 0:
            self._show(self.i - 1)

    def _ask_exit(self):
        if messagebox.askyesno(self.ui['title'], self.ui['exit_q'], parent=self.bubble):
            self._close()

    def _close(self):
        if self._after is not None:
            try:
                self.root.after_cancel(self._after)
            except Exception:
                pass
            self._after = None
        for w in [self.bubble] + self.edges:
            try:
                w.destroy()
            except Exception:
                pass
        _set_folders(export=self._old_export)
        self.app.tutorial_folder = None
        self.app.tutorial = None


# =================== PARAMETER HELP (the "i" button in the Analysis tab) ===================
# Each entry: ('section', en, he) for a heading, or (column name(s), en, he).
# {an}, {cath}, {junc}, {targets} are filled in from the current settings.
PARAM_HELP = [
    ('section', "General", "כללי"),
    ('', "All potentials are in V vs RHE (converted with the Reference and pH set in the "
         "toolbar). Currents are in mA, as in the .mpt file. One row = one cycle. A blank "
         "cell means the value could not be computed for that cycle.",
         "כל הפוטנציאלים ב-V מול RHE (מומרים לפי ה-Reference וה-pH שבסרגל הכלים). הזרמים "
         "ב-mA, כמו בקובץ ה-mpt. כל שורה היא מחזור אחד. תא ריק פירושו שלא ניתן היה לחשב "
         "את הערך במחזור הזה."),
    ('file, cycle', "The .mpt file and the cycle number in it.",
     "קובץ ה-mpt ומספר המחזור בתוכו."),

    ('section', "Peaks (raw curve)", "פיקים (עקומה גולמית)"),
    ('anodic_E_raw, anodic_I_raw',
     "Potential and current of the highest current inside the anodic window ({an}). "
     "This is the extreme point in the window, not a fitted peak.",
     "הפוטנציאל והזרם של הזרם הגבוה ביותר בתוך החלון האנודי ({an}). זו הנקודה הקיצונית "
     "בחלון, לא פיק מותאם."),
    ('cathodic_E_raw, cathodic_I_raw',
     "Potential and current of the lowest (most negative) current inside the cathodic "
     "window ({cath}).",
     "הפוטנציאל והזרם של הזרם הנמוך ביותר (השלילי ביותר) בתוך החלון הקתודי ({cath})."),
    ('I_ratio_raw', "anodic_I / |cathodic_I| - the ratio of the two peak currents.",
     "anodic_I / |cathodic_I| - היחס בין שני זרמי הפיק."),
    ('dEp_raw', "anodic_E - cathodic_E - peak separation, in V.",
     "anodic_E - cathodic_E - ההפרש בין פוטנציאלי הפיקים, ב-V."),

    ('section', "Charge", "מטען"),
    ('Q_end_raw',
     "Charge accumulated over the whole cycle: the integral of current over time, in "
     "mA·s (= mC). Net value - anodic charge counts positive, cathodic negative.",
     "המטען שהצטבר לאורך כל המחזור: האינטגרל של הזרם על הזמן, ב-mA·s (= mC). ערך נטו - "
     "מטען אנודי נספר חיובי וקתודי שלילי."),

    ('section', "OER-corrected values (*_corr)", "ערכים אחרי תיקון OER (*_corr)"),
    ('anodic_E_corr ... Q_end_corr, I_ratio_corr, dEp_corr',
     "The same quantities as the raw ones, measured on the OER-corrected curve. The "
     "correction (Mode A): where both sweeps carry positive current, the backward-sweep "
     "current at the same potential is subtracted from the forward-sweep current; "
     "elsewhere the curve is unchanged. These columns are always computed, whether or not "
     "'Subtract OER' is ticked in the ⋮ menu, and are blank when the correction cannot be "
     "made for that cycle.",
     "אותם גדלים כמו בערכים הגולמיים, אבל על העקומה המתוקנת ל-OER. התיקון (Mode A): "
     "בתחום שבו בשתי הסריקות הזרם חיובי, מחסרים מהזרם בסריקה קדימה את הזרם בסריקה "
     "אחורה באותו פוטנציאל; בשאר העקומה אין שינוי. העמודות האלה מחושבות תמיד, בלי קשר "
     "לסימון 'Subtract OER' בתפריט ⋮, ונשארות ריקות כשאי אפשר לבצע את התיקון במחזור."),

    ('section', "Onset potentials", "פוטנציאלי onset"),
    ('onset_raw_E_p5 / p10 / p15 / p25  (and onset_corr_E_p...)',
     "On the forward sweep, inside the anodic window: the first potential at which the "
     "current reaches 5 / 10 / 15 / 25 % of the way from the lowest to the highest "
     "current in the window.",
     "בסריקה קדימה, בתוך החלון האנודי: הפוטנציאל הראשון שבו הזרם מגיע ל-5 / 10 / 15 / "
     "25 % מהדרך בין הזרם הנמוך לגבוה בחלון."),
    ('onset_raw_E_k2 ... k20  (and onset_corr_E_k...)',
     "On the forward sweep, inside the anodic window: the positive part of the second "
     "derivative d²I/dE² is integrated cumulatively. The first 20 % of the points give a "
     "baseline (mean and standard deviation). The onset is the first potential where the "
     "integral stays above mean + k × std for 'min_run' points in a row (⋮ menu). Higher "
     "k = a more pronounced rise is required.",
     "בסריקה קדימה, בתוך החלון האנודי: מחשבים אינטגרל מצטבר של החלק החיובי של הנגזרת "
     "השנייה d²I/dE². 20 % הנקודות הראשונות משמשות כקו בסיס (ממוצע וסטיית תקן). ה-onset "
     "הוא הפוטנציאל הראשון שבו האינטגרל נשאר מעל ממוצע + k × סטיית תקן במשך 'min_run' "
     "נקודות ברצף (תפריט ⋮). k גבוה יותר = נדרשת עלייה בולטת יותר."),

    ('section', "Currents at fixed potentials", "זרמים בפוטנציאלים קבועים"),
    ('raw_I_fwd_<E>, raw_I_bwd_<E>, raw_I_for_back_<E>',
     "The cycle is split at its highest potential into a forward and a backward sweep. "
     "fwd / bwd = current on that sweep at potential E (interpolated); for_back = fwd - bwd. "
     "In the column name 'p' replaces the decimal point: 1p375 = 1.375 V. Potentials: "
     "{targets}. Blank if the sweep does not reach that potential.",
     "המחזור מחולק בפוטנציאל הגבוה ביותר שלו לסריקה קדימה ולסריקה אחורה. fwd / bwd = "
     "הזרם בסריקה הזו בפוטנציאל E (באינטרפולציה); for_back = fwd - bwd. בשם העמודה 'p' "
     "מחליף את הנקודה העשרונית: 1p375 = 1.375 V. הפוטנציאלים: {targets}. ריק אם הסריקה "
     "לא מגיעה לפוטנציאל הזה."),
    ('raw_junction_E, raw_junction_I_fwd, raw_junction_I_bwd, raw_junction_for_back',
     "Inside {junc}: the potential where the forward and backward currents are closest "
     "to each other (where the two sweeps meet), the two currents there, and their "
     "difference. If the sweeps do not cross in this range, it is the point where they "
     "come closest.",
     "בתוך {junc}: הפוטנציאל שבו הזרם בסריקה קדימה והזרם בסריקה אחורה הכי קרובים זה לזה "
     "(המקום שבו הסריקות נפגשות), שני הזרמים בנקודה הזו וההפרש ביניהם. אם הסריקות לא "
     "נחתכות בתחום, זו הנקודה שבה הן הכי קרובות."),

    ('section', "Report statistics (bottom pane)", "סטטיסטיקה בדוח (החלונית התחתונה)"),
    ('n, miss, miss%', "Number of cycles with a value, and how many / what share are blank.",
     "מספר המחזורים שיש בהם ערך, וכמה / איזה חלק ריקים."),
    ('mean, std, variance', "Average, sample standard deviation and variance over the cycles.",
     "ממוצע, סטיית תקן של המדגם ושונות על פני המחזורים."),
    ('min, max, median', "Smallest, largest and middle value.",
     "הערך הקטן, הגדול והאמצעי."),
    ('CV%', "std / |mean| × 100 - the spread relative to the value, comparable between "
            "parameters of different size.",
     "std / |mean| × 100 - הפיזור ביחס לערך, מאפשר להשוות פרמטרים בסדרי גודל שונים."),
]


class ParamHelpWindow(tk.Toplevel):
    """What every column of the analysis table means, in English or Hebrew."""

    def __init__(self, app, lang='en'):
        super().__init__(app.root)
        self.app = app
        self.title("Analysis parameters")
        self.geometry("680x640")
        top = tk.Frame(self)
        top.pack(fill='x', padx=8, pady=6)
        tk.Label(top, text="Analysis parameters",
                 font=('Arial', 12, 'bold')).pack(side='left')
        tk.Button(top, text="עברית", font=('Arial', 9),
                  command=lambda: self.render('he')).pack(side='right', padx=2)
        tk.Button(top, text="English", font=('Arial', 9),
                  command=lambda: self.render('en')).pack(side='right', padx=2)

        fr = tk.Frame(self)
        fr.pack(fill='both', expand=True, padx=8, pady=(0, 8))
        sb = tk.Scrollbar(fr, orient='vertical')
        sb.pack(side='right', fill='y')
        self.txt = tk.Text(fr, wrap='word', font=('Arial', 10), padx=10, pady=8,
                           yscrollcommand=sb.set, relief='flat', bg='#fbfbfb')
        self.txt.pack(side='left', fill='both', expand=True)
        sb.config(command=self.txt.yview)
        self.render(lang)

    def _fill_values(self):
        a = self.app
        try:
            an = f"{float(a.an_low.get()):g}-{float(a.an_high.get()):g} V"
            cath = f"{float(a.cath_low.get()):g}-{float(a.cath_high.get()):g} V"
        except Exception:
            an, cath = "1.35-1.50 V", "1.25-1.40 V"
        return dict(an=an, cath=cath,
                    junc=f"{JUNCTION_RANGE[0]:g}-{JUNCTION_RANGE[1]:g} V",
                    targets=", ".join(f"{e:g}" for e in TARGET_E) + " V")

    def render(self, lang):
        rtl = lang == 'he'
        T = self.txt
        T.config(state='normal')
        T.delete('1.0', 'end')
        side = 'right' if rtl else 'left'
        T.tag_configure('sec', font=('Arial', 12, 'bold'), foreground='#0a4d8c',
                        spacing1=14, spacing3=4, justify=side)
        T.tag_configure('name', font=('Consolas', 10, 'bold'), spacing1=8, justify=side)
        T.tag_configure('body', font=('Arial', 10), lmargin1=14, lmargin2=14,
                        rmargin=(14 if rtl else 0), spacing3=2, justify=side)
        vals = self._fill_values()
        for name, en, he in PARAM_HELP:
            text = (he if rtl else en).format(**vals)
            if name == 'section':
                T.insert('end', text + "\n", 'sec')
                continue
            if name:
                T.insert('end', name + "\n", 'name')
            T.insert('end', text + "\n", 'body')
        T.config(state='disabled')
        T.see('1.0')


# =================== GUI ===================
class CVApp:
    def __init__(self, root):
        self.root = root
        root.title("CV Analyzer v8.37")

        # Folders chosen under File > Set default folders... replace the
        # built-in defaults, so the script works on any computer.
        _st = overlay_load_settings()
        _set_folders(_st.get('mpt_folder'), _st.get('export_folder'))
        self.tutorial = None
        self.tutorial_folder = None

        # Fit the window to the actual screen (helps on 14" laptops, where
        # a hard-coded size would push controls off the bottom/right edge).
        try:
            sw = root.winfo_screenwidth()
            sh = root.winfo_screenheight()
            root.geometry(f"{int(sw * 0.95)}x{int(sh * 0.90)}+10+10")
            root.minsize(900, 560)
        except Exception:
            pass

        # Data & state
        self.files = []
        self.cache = {}       # path -> dict(E,I,t,cyc,ocv)
        self.file_nodes = {}
        self.checked = {}
        self.file_colors = {}  # per-file default, or per-cycle override with key (path, cyc)
        self.display_names = {}
        self.dir_nodes = {}
        self.scan_rates = {}   # (path, cyc) -> V/s

        # Peaks: (file_path, cycle) -> list of {'E','Y','type'}
        self.peaks = {}
        self.peak_index_map = []

        # Recent folders
        self.recents = []
        # Root folders currently listed in the tree (used to reload after a
        # batch rename, since renaming invalidates every stored path).
        self.open_roots = []

        # Modes
        self.manual_peak_mode = False
        self.delete_peak_mode = False
        self.last_debug = {}

        # ---------- Smoothing settings ----------
        self.use_smoothing_for_peaks = tk.BooleanVar(value=False)
        self.use_smoothing_for_derivatives = tk.BooleanVar(value=True)
        self.smooth_method = tk.StringVar(value='None')  # None | Moving average | Savitzky–Golay
        self.smooth_window = tk.StringVar(value='9')     # odd recommended for Savgol
        self.savgol_poly = tk.StringVar(value='2')
        
        # Live-derivatives tracking (for instant updates when smoothing changes)
        self._deriv_last_pair = None  # (path, cyc) or None

        # Sliders (range bars) mirror into IntVars for Tk Scale widgets
        try:
            _w0 = int(float(self.smooth_window.get()))
        except Exception:
            _w0 = 9
        try:
            _p0 = int(float(self.savgol_poly.get()))
        except Exception:
            _p0 = 2

        self.smooth_window_i = tk.IntVar(value=max(1, _w0))
        self.savgol_poly_i = tk.IntVar(value=max(0, _p0))
        
        # ---------- 2nd-derivative integrator detection (tunable) ----------
        self.int_k_sigma = tk.DoubleVar(value=3.0)      # significance threshold in sigma units
        self.int_min_run = tk.IntVar(value=15)          # consecutive points above threshold
        
        # Onset-by-kσ markers (raw plot)
        self.show_onset_k_lines = tk.BooleanVar(value=False)
        
        # minimum consecutive points above threshold to accept onset (tuneable)
        self.onset_min_run = tk.IntVar(value=4)




        # ---------- Analysis / Comparison state ----------
        self.analysis_rows = []   # list of dicts
        self.analysis_map = {}    # (path, cyc) -> dict of metrics

        # Backward-compatible single-param (kept), plus multi-select
        self.cmp_param = tk.StringVar(value='anodic_E_raw')

        # -------------------- MENU --------------------
        menu = tk.Menu(root)
        root.config(menu=menu)
        file_menu = tk.Menu(menu, tearoff=False)
        menu.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Open Folder...", command=self.open_folder)
        file_menu.add_command(label="Add MPT...", command=self.add_file)
        self.recent_menu = tk.Menu(file_menu, tearoff=False)
        file_menu.add_cascade(label="Open Recent", menu=self.recent_menu)
        file_menu.add_separator()
        file_menu.add_command(label="Save PNG...", command=self.save_png)
        file_menu.add_separator()
        file_menu.add_command(label="Set default folders...", command=self.set_default_folders)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=root.quit)
        help_menu = tk.Menu(menu, tearoff=False)
        menu.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="Tutorial...", command=self.start_tutorial)

        # -------------------- TOP BAR --------------------
        # Split into two rows so everything stays reachable on a 14" screen
        # instead of overflowing off the right edge of a single row.
        bar_wrap = tk.Frame(root)
        bar_wrap.pack(fill='x')
        bar = tk.Frame(bar_wrap)
        bar.pack(fill='x')
        bar2 = tk.Frame(bar_wrap)
        bar2.pack(fill='x')

        self.btn_open_folder = tk.Button(bar, text="Open Folder", command=self.open_folder)
        self.btn_open_folder.pack(side='left')
        tk.Button(bar, text="Add MPT", command=self.add_file).pack(side='left', padx=(4,0))
        tk.Button(bar, text="Clear All", command=self.clear_all).pack(side='left', padx=4)

        tk.Label(bar, text='Reference:').pack(side='left')
        self.ref_var = tk.StringVar(value='RHE')
        self.ref_combo = ttk.Combobox(bar, textvariable=self.ref_var, width=18,
                     values=['RHE','Ag/AgCl 3M KCl','Ag/AgCl sat. KCl',
                             'Hg/HgO 1M KOH','Hg/HgO 0.1M KOH','SCE'])
        self.ref_combo.pack(side='left')

        tk.Label(bar, text='pH:').pack(side='left')
        self.ph_var = tk.StringVar(value='14')
        tk.Entry(bar, textvariable=self.ph_var, width=4).pack(side='left')

        tk.Label(bar, text='Font:').pack(side='left')
        self.font_var = tk.StringVar(value='12')
        tk.Entry(bar, textvariable=self.font_var, width=4).pack(side='left')

        tk.Label(bar, text='Thick:').pack(side='left')
        self.line_var = tk.StringVar(value='2')
        tk.Entry(bar, textvariable=self.line_var, width=4).pack(side='left')

        tk.Label(bar, text="Plot Mode:").pack(side='left', padx=(10,0))
        self.plot_mode = tk.StringVar(value='IvsE')
        self.plot_mode_combo = ttk.Combobox(bar, textvariable=self.plot_mode, width=8,
                     values=['IvsE','QvsE','IvsT','EvsT'])
        self.plot_mode_combo.pack(side='left')

        self.time_zero = tk.BooleanVar(value=True)
        tk.Checkbutton(bar, text="t=0", variable=self.time_zero, command=lambda: self.plot()).pack(side='left', padx=(6,0))

        self.hide_legend = tk.BooleanVar(value=True)
        tk.Checkbutton(bar, text="Hide Legend", variable=self.hide_legend, command=self.plot).pack(side='left')

        self.subtract_oer = tk.BooleanVar(value=False)
        self.show_onset_lines = tk.BooleanVar(value=False)
        self.sanity_mode = tk.BooleanVar(value=True)

        # Rarely used options live in the "⋮" menu at the right end of the toolbar.
        self.btn_more = tk.Menubutton(bar, text="⋮", relief='flat', padx=10,
                                      font=('TkDefaultFont', 13, 'bold'), cursor='hand2')
        self.btn_more.pack(side='right', padx=4)
        more = tk.Menu(self.btn_more, tearoff=False)
        self.btn_more.config(menu=more)
        more.add_checkbutton(label="Subtract OER (Mode A)",
                             variable=self.subtract_oer, command=self.plot)
        more.add_separator()
        more.add_checkbutton(label="Show onset % (raw)",
                             variable=self.show_onset_lines, command=self.plot)
        more.add_checkbutton(label="Show onset kσ (raw)",
                             variable=self.show_onset_k_lines, command=self.plot)
        more.add_command(label="Onset min_run...", command=self._ask_onset_min_run)
        more.add_separator()
        more.add_checkbutton(label="Sanity Fix Mode (skip NaN points)",
                             variable=self.sanity_mode, command=self.plot)

        self.btn_plot = tk.Button(bar2, text="Plot", command=self.plot)
        self.btn_plot.pack(side='left', padx=6)
        tk.Button(bar2, text="Save PNG", command=self.save_png).pack(side='left')

        

        # -------------------- SIDE-BY-SIDE LAYOUT --------------------
        body = tk.Frame(root)
        body.pack(fill='both', expand=True)

        # Thin always-visible strip holding the collapse/expand toggle.
        # It stays put when the tree panel is hidden, so there is always a
        # way to bring the panel back.
        strip = tk.Frame(body, width=16)
        strip.pack(side='left', fill='y')
        strip.pack_propagate(False)
        self.collapse_btn = tk.Button(
            strip, text='◀', relief='flat', bd=0, padx=0,
            command=self.toggle_left_panel
        )
        self.collapse_btn.pack(fill='both', expand=True)

        main = tk.PanedWindow(body, orient='horizontal', sashrelief='raised', sashwidth=6)
        main.pack(fill='both', expand=True)
        self.main_paned = main

        # Left column (Tree)
        left = tk.Frame(main)
        main.add(left, minsize=220)
        self.left_panel = left
        self._left_visible = True

        style = ttk.Style()
        style.configure('Treeview', rowheight=22)
        style.map('Treeview', background=[('selected', '#cce6ff')])

        # Library toolbar
        lib_bar = tk.Frame(left)
        lib_bar.pack(fill='x', pady=(2, 2))
        tk.Button(lib_bar, text="Rename MPT files...",
                  command=self.batch_rename_mpt).pack(side='left', padx=2)
        self.btn_find_unanalyzed = tk.Button(lib_bar, text="Find unanalyzed...",
                                             command=self.find_unanalyzed_mpt)
        self.btn_find_unanalyzed.pack(side='left', padx=2)
        self.btn_overlay = tk.Button(lib_bar, text="Overlay...",
                                     command=self.open_overlay_window)
        self.btn_overlay.pack(side='left', padx=2)

        tree_frame = tk.Frame(left)
        tree_frame.pack(fill="both", expand=True)

        tree_scroll_y = tk.Scrollbar(tree_frame, orient="vertical")
        tree_scroll_x = tk.Scrollbar(tree_frame, orient="horizontal")

        self.tree = ttk.Treeview(
            tree_frame,
            columns=('chk','type','color','path'),
            displaycolumns=('chk','type','color'),  # 'path' still stored, just not shown
            selectmode='none',
            yscrollcommand=tree_scroll_y.set,
            xscrollcommand=tree_scroll_x.set
        )
        self.tree.heading('#0', text="Folder / File / Cycle")
        self.tree.heading('chk', text="✔")
        self.tree.heading('type', text="Type")
        self.tree.heading('color', text="Color")
        self.tree.heading('path', text="Path")
        self.tree.column('#0', width=190, stretch=True)
        self.tree.column('chk', width=30, anchor='center', stretch=False)
        self.tree.column('type', width=56, anchor='center', stretch=False)
        self.tree.column('color', width=46, anchor='center', stretch=False)
        self.tree.column('path', width=160, anchor='w', stretch=True)

        tree_scroll_y.config(command=self.tree.yview)
        tree_scroll_y.pack(side="right", fill="y")
        tree_scroll_x.config(command=self.tree.xview)
        tree_scroll_x.pack(side="bottom", fill="x")
        self.tree.pack(side="left", fill="both", expand=True)

        self.tree.bind('<Button-1>', self.on_tree_click)
        self.tree.bind('<Double-1>', self.rename_item)
        self.tree.bind('<<TreeviewOpen>>', self.on_tree_expand)

        # Center column (Main Plot)
        center = tk.Frame(main)
        main.add(center, minsize=420)
        self.center_panel = center

        fig = Figure(figsize=(6, 4))
        self.ax = fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(fig, center)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        toolbar = NavigationToolbar2Tk(self.canvas, center)
        toolbar.update()
        self.canvas.mpl_connect('button_press_event', self.on_click_plot)

        # Right column (Tabbed)
        right = tk.Frame(main)
        main.add(right, minsize=340)

        tabs = ttk.Notebook(right)
        tabs.pack(fill='both', expand=True)

       # ---------- Peak detection params (still used by the analysis code) ----------
        self.min_height_var = tk.StringVar(value='0.001')
        self.cath_low = tk.StringVar(value='1.25')
        self.cath_high = tk.StringVar(value='1.40')
        self.an_low = tk.StringVar(value='1.35')
        self.an_high = tk.StringVar(value='1.50')

        # ---------- Tafel params (state kept; Tafel tab removed) ----------
        self.tafel_E0 = tk.StringVar(value='1.48')
        self.tafel_N = tk.StringVar(value='20')
        self.tafel_logI_min = tk.StringVar(value='-6')
        self.tafel_logI_max = tk.StringVar(value='-2')

        # Hidden holders. Several existing methods write into these widgets
        # (peak list, debug log, tafel log). They are created but never
        # packed, so removing the Peaks/Debug/Tafel/Smoothing/Comparison
        # tabs cannot raise AttributeError anywhere in the codebase.
        self._hidden = tk.Frame(right)
        self.peak_list = tk.Listbox(self._hidden)
        self.debug_text = tk.Text(self._hidden)
        self.tafel_text = tk.Text(self._hidden)

        # -------- Analysis tab (the only tab) --------
        analysis_tab = tk.Frame(tabs)
        tabs.add(analysis_tab, text="Analysis")
        self._build_analysis_tab(analysis_tab)

        # Colors
        self.color_cycle = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']
        
    

        if os.path.isdir(DEFAULT_FOLDER):
            self.open_folder(path_override=DEFAULT_FOLDER)

        # Any change of the plot settings re-draws the main plot by itself
        # (debounced, so typing "14" in pH draws once).
        self._plot_after = None
        for v in (self.ref_var, self.ph_var, self.font_var, self.line_var, self.plot_mode):
            v.trace_add('write', self._schedule_plot)

        # First run: the tutorial starts by itself until it is finished once.
        root.after(700, self._maybe_start_tutorial)

    def _schedule_plot(self, *_):
        if self._plot_after is not None:
            try:
                self.root.after_cancel(self._plot_after)
            except Exception:
                pass
        self._plot_after = self.root.after(350, self._auto_plot)

    def _auto_plot(self):
        self._plot_after = None
        try:
            self.plot()
        except Exception:
            pass   # e.g. a half-typed value - the next change re-plots

    def _ask_onset_min_run(self):
        try:
            cur = int(self.onset_min_run.get())
        except Exception:
            cur = 4
        v = simpledialog.askinteger("Onset min_run",
                                    "Minimum consecutive points above the threshold:",
                                    initialvalue=cur, minvalue=1, maxvalue=500, parent=self.root)
        if v:
            self.onset_min_run.set(v)
            self.plot()

    def offer_folder_setup(self, lang='en'):
        """After the tutorial: point the program at the user's data folder."""
        if overlay_load_settings().get('mpt_folder'):
            return
        if FolderSetupDialog(self.root, lang).choose:
            self.set_default_folders()

    # =================== Default folders / tutorial ===================
    def set_default_folders(self):
        """Choose (and remember) the .mpt folder and the CSV export folder."""
        if self.tutorial is not None:
            messagebox.showinfo("Default folders", "Finish or exit the tutorial first.")
            return
        mpt = filedialog.askdirectory(
            title="Folder with the .mpt files",
            initialdir=DEFAULT_FOLDER if os.path.isdir(DEFAULT_FOLDER) else find_shared_drives_root())
        if not mpt:
            return
        exp = filedialog.askdirectory(
            title="Folder for the exported CSV files",
            initialdir=EXPORT_FOLDER if os.path.isdir(EXPORT_FOLDER) else mpt)
        if not exp:
            return
        st = overlay_load_settings()
        st['mpt_folder'] = mpt
        st['export_folder'] = exp
        overlay_save_settings(st)
        _set_folders(mpt, exp)
        if messagebox.askyesno("Default folders",
                               f"Saved.\n\n.mpt files:\n{mpt}\n\nExports:\n{exp}\n\n"
                               "Open the .mpt folder now?"):
            self.clear_all()
            self.open_folder(path_override=mpt)

    def _maybe_start_tutorial(self):
        if not overlay_load_settings().get('tutorial_done'):
            self.start_tutorial(first_run=True)

    def start_tutorial(self, first_run=False):
        if self.tutorial is not None:
            try:
                self.tutorial.bubble.lift()
                return
            except Exception:
                self.tutorial = None
        lang = TutorialLanguageDialog(self.root, first_run).result
        if not lang:
            return
        # Fresh view, so only the sample files are listed (no file is changed).
        ow = getattr(self, '_overlay_win', None)
        if ow is not None:
            try:
                ow.destroy()
            except Exception:
                pass
            self._overlay_win = None
        self.clear_all()
        self.tutorial = Tutorial(self, lang)
            
    def _cumtrapz(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        n = min(len(x), len(y))
        if n < 2:
            return np.zeros(n, dtype=float)
        x = x[:n]; y = y[:n]
        dx = np.diff(x)
        inc = 0.5 * (y[1:] + y[:-1]) * dx
        out = np.zeros(n, dtype=float)
        out[1:] = np.cumsum(inc)
        out[~np.isfinite(out)] = 0.0
        return out
    
    def _find_sig_positive_E(
        self, E, d2,
        baseline_frac=0.15,   # fraction of early points used for baseline stats
        k_sigma=4.0,          # significance strength
        min_run_frac=0.05     # require a short consecutive run above threshold
    ):
        E = np.asarray(E, dtype=float)
        d2 = np.asarray(d2, dtype=float)
        n = min(len(E), len(d2))
        if n < 10:
            return np.nan, np.array([], dtype=float), np.nan
    
        E = E[:n]; d2 = d2[:n]
    
        # cumulative integral of 2nd derivative
        C = self._cumtrapz(E, d2)
    
        # baseline stats on early scan
        nb = max(5, int(baseline_frac * n))
        base = C[:nb]
        mu = float(np.nanmean(base)) if np.any(np.isfinite(base)) else 0.0
        sd = float(np.nanstd(base))  if np.any(np.isfinite(base)) else 0.0
    
        thr = mu + k_sigma * sd
    
        # first index where C crosses threshold and stays above it for a run
        min_run = max(3, int(min_run_frac * n))
        above = C > thr
        if np.any(above):
            idxs = np.where(above)[0]
            for i in idxs:
                j = min(n, i + min_run)
                if np.all(above[i:j]):
                    return float(E[i]), C, thr
    
        return np.nan, C, thr
    
    def _align_zeros(self, ax_left, ax_right):
        y1min, y1max = ax_left.get_ylim()
        y2min, y2max = ax_right.get_ylim()

        def frac_zero(ymin, ymax):
            if ymax == ymin:
                return 0.5
            return (0 - ymin) / (ymax - ymin)

        f1 = frac_zero(y1min, y1max)
    
        # ensure right axis includes zero
        if not (y2min < 0 < y2max):
            m = max(abs(y2min), abs(y2max), 1e-9)
            y2min, y2max = -m, m
    
        span2 = (y2max - y2min)
        new_y2min = -f1 * span2
        new_y2max = (1 - f1) * span2
        ax_right.set_ylim(new_y2min, new_y2max)

    
    def _integrator_onset_from_d2(self, E, d2):
        """
        Finds onset E where cumulative integral of positive part of d2 becomes
        significantly positive for a sustained run.
    
        Uses GUI-tunable:
          - self.int_k_sigma (float)
          - self.int_min_run (int points)
    
        Returns: onset_E (float or np.nan), debug dict
        """
        E = np.asarray(E, dtype=float)
        d2 = np.asarray(d2, dtype=float)
    
        dbg = {"onset_E": np.nan, "k_sigma": float(self.int_k_sigma.get()),
               "min_run": int(self.int_min_run.get()), "thr": np.nan}
    
        if len(E) < 10 or len(d2) < 10:
            return np.nan, dbg
    
        # sort by E (important!)
        o = np.argsort(E)
        Eo = E[o]
        d2o = d2[o]
    
        # integrate only positive curvature contribution
        dE = np.gradient(Eo)
        dE[~np.isfinite(dE)] = 0.0
        pos = np.maximum(d2o, 0.0)
    
        integ = np.cumsum(pos * np.abs(dE))
    
        # baseline: early portion (first 20% of points, at least 10)
        n0 = max(10, int(0.2 * len(integ)))
        base = integ[:n0]
        mu = float(np.nanmean(base))
        sd = float(np.nanstd(base))
        if not np.isfinite(mu) or not np.isfinite(sd) or sd == 0:
            return np.nan, dbg
    
        k = float(self.int_k_sigma.get())
        thr = mu + k * sd
        dbg["thr"] = thr
    
        above = integ >= thr
    
        # require sustained run
        min_run = int(self.int_min_run.get())
        min_run = max(1, min_run)
    
        run = 0
        for i, a in enumerate(above):
            if a:
                run += 1
                if run >= min_run:
                    onset_E = float(Eo[i - min_run + 1])
                    dbg["onset_E"] = onset_E
                    return onset_E, dbg
            else:
                run = 0
    
        return np.nan, dbg


    def _build_peaks_tab(self, tab):
        # ---------- Controls ----------
        ctrl = tk.LabelFrame(tab, text="Peak Detection / Marking")
        ctrl.pack(fill='x', padx=8, pady=8)
    
        r0 = tk.Frame(ctrl); r0.pack(fill='x', padx=8, pady=(6,2))
        tk.Button(r0, text="Run Auto Peaks (checked cycles)", command=self.run_auto_peaks).pack(side='left')
        tk.Button(r0, text="Manual Peak Mode", command=self.enable_manual_peak_mode).pack(side='left', padx=(8,0))
        tk.Button(r0, text="Delete Peak Mode", command=self.enable_delete_peak_mode).pack(side='left', padx=(8,0))
    
        r1 = tk.Frame(ctrl); r1.pack(fill='x', padx=8, pady=2)
        tk.Label(r1, text="Min Height:").pack(side='left')
        tk.Entry(r1, textvariable=self.min_height_var, width=8).pack(side='left', padx=(6,16))
    
        r2 = tk.Frame(ctrl); r2.pack(fill='x', padx=8, pady=2)
        tk.Label(r2, text="Cathodic window (V):").pack(side='left')
        tk.Entry(r2, textvariable=self.cath_low, width=7).pack(side='left', padx=(6,4))
        tk.Label(r2, text="to").pack(side='left')
        tk.Entry(r2, textvariable=self.cath_high, width=7).pack(side='left', padx=(4,16))
    
        tk.Label(r2, text="Anodic window (V):").pack(side='left')
        tk.Entry(r2, textvariable=self.an_low, width=7).pack(side='left', padx=(6,4))
        tk.Label(r2, text="to").pack(side='left')
        tk.Entry(r2, textvariable=self.an_high, width=7).pack(side='left', padx=(4,0))
    
        hint = tk.Label(
            ctrl,
            text=("Manual mode: click on a checked cycle curve to add a peak. "
                  "Delete mode: click near a labeled peak to remove it."),
            fg='#444'
        )
        hint.pack(anchor='w', padx=8, pady=(4,6))
    
        # ---------- Peak list ----------
        body = tk.Frame(tab)
        body.pack(fill='both', expand=True, padx=8, pady=(0,8))
    
        left = tk.Frame(body)
        left.pack(side='left', fill='both', expand=True)
    
        tk.Label(left, text="Detected / Manual Peaks:").pack(anchor='w')
        self.peak_list = tk.Listbox(left, height=18, selectmode='extended')
        self.peak_list.pack(side='left', fill='both', expand=True)
    
        sb = tk.Scrollbar(left, orient='vertical', command=self.peak_list.yview)
        sb.pack(side='left', fill='y')
        self.peak_list.config(yscrollcommand=sb.set)
    
        right = tk.Frame(body)
        right.pack(side='left', fill='y', padx=(10,0))
    
        tk.Button(right, text="Delete Selected Peaks", command=self.delete_selected_peak).pack(fill='x', pady=(0,6))
        tk.Button(right, text="Refresh List", command=self.refresh_peak_panel).pack(fill='x')


    # =================== Smoothing UI ===================
    def _build_smoothing_tab(self, tab):
        top = tk.Frame(tab); top.pack(fill='x', padx=8, pady=8)
    
        tk.Checkbutton(
            top, text="Use smoothing for peak detection",
            variable=self.use_smoothing_for_peaks,
            command=lambda: self._on_smoothing_change(live=True)
        ).pack(anchor='w')
    
        tk.Checkbutton(
            top, text="Use smoothing for derivatives",
            variable=self.use_smoothing_for_derivatives,
            command=lambda: self._on_smoothing_change(live=True)
        ).pack(anchor='w')
    
        row = tk.Frame(tab); row.pack(fill='x', padx=8, pady=(6,6))
        tk.Label(row, text="Method:").pack(side='left')
        cmb = ttk.Combobox(
            row, textvariable=self.smooth_method, width=18,
            values=['None', 'Moving average', 'Savitzky-Golay'],
            state='readonly'
        )
        cmb.pack(side='left', padx=(6,12))
        cmb.bind("<<ComboboxSelected>>", lambda e: self._on_smoothing_change(live=True))
    
        # Window slider
        srow = tk.Frame(tab); srow.pack(fill='x', padx=8, pady=(4,2))
        tk.Label(srow, text="Window:").pack(side='left')
        self.window_val_lbl = tk.Label(srow, text=str(self.smooth_window_i.get()), width=4)
        self.window_val_lbl.pack(side='right')
        tk.Scale(
            srow, from_=1, to=301, orient='horizontal',
            variable=self.smooth_window_i, showvalue=False,
            command=lambda v: self._on_smoothing_change(live=True)
        ).pack(side='left', fill='x', expand=True, padx=(6,8))
    
        # Poly slider
        prow = tk.Frame(tab); prow.pack(fill='x', padx=8, pady=(4,2))
        tk.Label(prow, text="Poly (SG):").pack(side='left')
        self.poly_val_lbl = tk.Label(prow, text=str(self.savgol_poly_i.get()), width=4)
        self.poly_val_lbl.pack(side='right')
        tk.Scale(
            prow, from_=0, to=9, orient='horizontal',
            variable=self.savgol_poly_i, showvalue=False,
            command=lambda v: self._on_smoothing_change(live=True)
        ).pack(side='left', fill='x', expand=True, padx=(6,8))
    
        hint = tk.Label(
            tab,
            text=("Savitzky–Golay requires SciPy (scipy.signal). "
                  + ("SciPy detected ✅" if HAS_SAVGOL else "SciPy NOT detected ❌ (SG will act as None).")),
            fg=('green' if HAS_SAVGOL else 'red')
        )
        hint.pack(anchor='w', padx=8, pady=(6,8))
        
        # ---------- Integrator detection tuning ----------
        int_frame = tk.LabelFrame(tab, text="2nd-derivative integrator onset (tuning)")
        int_frame.pack(fill='x', padx=8, pady=(6, 8))
        
        rA = tk.Frame(int_frame); rA.pack(fill='x', padx=8, pady=(6,2))
        tk.Label(rA, text="k_sigma (higher = more significant):").pack(side='left')
        self.k_sigma_lbl = tk.Label(rA, text=f"{self.int_k_sigma.get():.2f}", width=6)
        self.k_sigma_lbl.pack(side='right')
        tk.Scale(
            rA, from_=0.5, to=15.0, resolution=0.1, orient='horizontal',
            variable=self.int_k_sigma, showvalue=False,
            command=lambda v: (self.k_sigma_lbl.config(text=f"{self.int_k_sigma.get():.2f}"),
                               self._replot_derivatives_smart())
        ).pack(side='left', fill='x', expand=True, padx=(8,8))
        
        rB = tk.Frame(int_frame); rB.pack(fill='x', padx=8, pady=(4,6))
        tk.Label(rB, text="min_run (consecutive points):").pack(side='left')
        self.min_run_lbl = tk.Label(rB, text=str(self.int_min_run.get()), width=6)
        self.min_run_lbl.pack(side='right')
        tk.Scale(
            rB, from_=3, to=200, orient='horizontal',
            variable=self.int_min_run, showvalue=False,
            command=lambda v: (self.min_run_lbl.config(text=str(self.int_min_run.get())),
                               self._replot_derivatives_smart())
        ).pack(side='left', fill='x', expand=True, padx=(8,8))

    
        btns = tk.Frame(tab); btns.pack(fill='x', padx=8, pady=(2,6))
        tk.Button(btns, text="Re-plot main", command=self.plot).pack(side='left')
        tk.Button(btns, text="Re-run Auto Peaks", command=self.run_auto_peaks).pack(side='left', padx=(8,0))
        tk.Button(btns, text="Plot derivatives (first checked)", command=self.plot_derivatives_for_first_checked).pack(side='left', padx=(8,0))
    
        # -------------------- Derivatives plot (NOW IN SMOOTHING TAB) --------------------
        # -------------------- Derivatives plots (3 stacked) --------------------
        plot_frame = tk.Frame(tab)
        plot_frame.pack(fill='both', expand=True, padx=8, pady=(0,8))
        
        fig = Figure(figsize=(6, 7.2))
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312, sharex=ax1)
        ax3 = fig.add_subplot(313, sharex=ax1)
        
        self.deriv_axes = [ax1, ax2, ax3]
        self.deriv_canvas = FigureCanvasTkAgg(fig, plot_frame)
        self.deriv_canvas.get_tk_widget().pack(fill='both', expand=True)
        
        tb = NavigationToolbar2Tk(self.deriv_canvas, plot_frame)
        tb.update()
        


    def _get_smooth_params(self):
        method = self.smooth_method.get()
        w = int(self.smooth_window_i.get())
        w = max(1, w)
        if method == 'Savitzky-Golay' and w % 2 == 0:
            w += 1
        return method, w
    

    def _maybe_smooth(self, y, for_what='peaks'):
        """
        for_what: 'peaks' or 'derivatives'
        """
        y = np.asarray(y, dtype=float)
        if len(y) == 0:
            return y

        if for_what == 'peaks' and (not self.use_smoothing_for_peaks.get()):
            return y
        if for_what == 'derivatives' and (not self.use_smoothing_for_derivatives.get()):
            return y

        method, w = self._get_smooth_params()

        if method == 'None':
            return y

        if method == 'Moving average':
            return moving_average(y, w)

        if method == 'Savitzky-Golay':
            if not HAS_SAVGOL:
                return y
            try:
                poly = int(float(self.savgol_poly.get()))
            except Exception:
                poly = 2
            poly = max(0, poly)
            if w <= poly:
                w = poly + 3
                if w % 2 == 0:
                    w += 1
            if w < 3:
                w = 3
            w = min(w, len(y) if len(y) % 2 == 1 else len(y) - 1)
            if w < 3:
                return y
            try:
                return savgol_filter(y, w, poly)
            except Exception:
                return y

        return y
    
    def _sync_smoothing_stringvars_from_sliders(self):
        # keep the original StringVars updated (used by _get_smooth_params and existing code)
        self.smooth_window.set(str(int(self.smooth_window_i.get())))
        self.savgol_poly.set(str(int(self.savgol_poly_i.get())))
    
        # update small numeric labels if exist
        if hasattr(self, 'window_val_lbl'):
            self.window_val_lbl.config(text=str(int(self.smooth_window_i.get())))
        if hasattr(self, 'poly_val_lbl'):
            self.poly_val_lbl.config(text=str(int(self.savgol_poly_i.get())))

    def _on_smoothing_change(self, live=True):
        # called by sliders / method dropdown / checkbox
        # enforce odd window in UI when SG selected
        if self.smooth_method.get() == 'Savitzky-Golay':
            w = int(self.smooth_window_i.get())
            if w % 2 == 0:
                self.smooth_window_i.set(w + 1)
    
        self._sync_smoothing_stringvars_from_sliders()
    
        if not live:
            return
    
        self._replot_derivatives_smart()


    def _replot_derivatives_smart(self):
        # if we already plotted derivatives before, replot same pair
        if self._deriv_last_pair is not None:
            p, cyc = self._deriv_last_pair
            return self._plot_derivatives_pair(p, cyc)
    
        # otherwise: use first checked cycle
        p, cyc = self._first_checked_cycle()
        if p is None:
            return
        return self._plot_derivatives_pair(p, cyc)
    
    def _plot_derivatives_pair(self, p, cyc):
        if not hasattr(self, 'deriv_axes') or not hasattr(self, 'deriv_canvas'):
            return
    
        # -------- preserve zoom/pan per subplot (left axes only) --------
        prev_limits = []
        for ax in self.deriv_axes:
            try:
                if len(ax.lines) > 0:
                    prev_limits.append((ax.get_xlim(), ax.get_ylim()))
                else:
                    prev_limits.append((None, None))
            except Exception:
                prev_limits.append((None, None))
    
        try:
            pH = float(self.ph_var.get())
        except Exception:
            pH = 14.0
    
        E_rhe, t_sel, Y_raw, Y_corr, mask_corr, I_for_peaks, I_sel, dbg, is_ocv = self.get_display_series(p, cyc, pH)
        if E_rhe is None or is_ocv:
            return
    
        E = np.asarray(E_rhe, dtype=float)
        I_raw = np.asarray(I_sel, dtype=float)
    
        # Current used for derivatives (keeps your checkbox behavior)
        if self.subtract_oer.get():
            I_base, _, _ = oer_subtract_modeA(E, I_raw)
        else:
            I_base = I_raw.copy()
    
        # Clean finite points (keep original order!)
        E_use, I_base2, ok_mask, _ = clean_xy(E, I_base)
        if E_use is None or len(E_use) < 8:
            return
    
        # ---------- helpers ----------
        def _get_k_sigma_and_min_run():
            # Use your GUI vars if you already created them; otherwise defaults
            k_sigma = 3.0
            min_run = 6
            for name in ("k_sigma_var", "k_sigma", "int_k_sigma", "int_k_sigma_var"):
                if hasattr(self, name):
                    try:
                        v = getattr(self, name)
                        k_sigma = float(v.get() if hasattr(v, "get") else v)
                        break
                    except Exception:
                        pass
            for name in ("min_run_var", "min_run", "int_min_run", "int_min_run_var"):
                if hasattr(self, name):
                    try:
                        v = getattr(self, name)
                        min_run = int(float(v.get() if hasattr(v, "get") else v))
                        break
                    except Exception:
                        pass
            min_run = max(2, min_run)
            k_sigma = max(0.0, k_sigma)
            return k_sigma, min_run
    
        def _cumulative_trapz(y, x):
            y = np.asarray(y, dtype=float)
            x = np.asarray(x, dtype=float)
            n = len(y)
            if n < 2:
                return np.zeros_like(y)
            dx = np.diff(x)
            inc = 0.5 * (y[1:] + y[:-1]) * dx
            out = np.zeros(n, dtype=float)
            out[1:] = np.cumsum(inc)
            out[~np.isfinite(out)] = 0.0
            return out
    
        def _detect_positive_rise(Eo, A, k_sigma=3.0, min_run=6):
            """
            Find first index where A becomes significantly positive for >=min_run consecutive samples.
            Baseline = first 20% of points.
            """
            Eo = np.asarray(Eo, dtype=float)
            A = np.asarray(A, dtype=float)
            n = len(A)
            if n < max(10, min_run + 3):
                return None
    
            bN = max(5, int(0.2 * n))
            base = A[:bN]
            mu = float(np.nanmedian(base))
            sd = float(np.nanstd(base))
            if not np.isfinite(sd) or sd <= 1e-12:
                return None
    
            z = (A - mu) / sd
            thr = k_sigma
    
            run = 0
            for i in range(n):
                if np.isfinite(z[i]) and (z[i] >= thr) and (A[i] > 0):
                    run += 1
                    if run >= min_run:
                        j = i - min_run + 1
                        return j
                else:
                    run = 0
            return None
    
        def _clear_twins(ax):
            # remove any twinx axes previously created for this subplot
            fig = ax.figure
            to_remove = []
            for a in list(fig.axes):
                if a is ax:
                    continue
                try:
                    same_pos = (a.get_position().bounds == ax.get_position().bounds)
                except Exception:
                    same_pos = False
                if same_pos:
                    if a.yaxis.get_ticks_position() in ('right', 'both'):
                        to_remove.append(a)
            for a in to_remove:
                try:
                    a.remove()
                except Exception:
                    pass
    
        def _compute_d1(Eseg, Iseg):
            Eseg = np.asarray(Eseg, dtype=float)
            Iseg = np.asarray(Iseg, dtype=float)
            if len(Eseg) < 5:
                return Eseg, np.zeros_like(Iseg)
            Is = self._maybe_smooth(Iseg, for_what='derivatives')
            d1 = safe_gradient(Is, Eseg)
            return Eseg, d1
    
        def _compute_d2(Eseg, Iseg):
            Eseg = np.asarray(Eseg, dtype=float)
            Iseg = np.asarray(Iseg, dtype=float)
            if len(Eseg) < 6:
                return Eseg, np.zeros_like(Iseg)
            Is = self._maybe_smooth(Iseg, for_what='derivatives')
            d1 = safe_gradient(Is, Eseg)
            d2 = safe_gradient(d1, Eseg)
            return Eseg, d2
    
        # ---------- split forward/back by switching potential index (in ORIGINAL order) ----------
        idx_sw = int(np.argmax(E_use))
        if idx_sw < 4 or idx_sw > len(E_use) - 5:
            # too close to edge -> unstable split
            idx_sw = len(E_use) // 2
    
        Ef = E_use[:idx_sw+1]
        If = I_base2[:idx_sw+1]
    
        Eb = E_use[idx_sw:]
        Ib = I_base2[idx_sw:]
    
        # Backward branch is typically decreasing E; reverse so x increases for plots/integration
        if len(Eb) > 2 and Eb[0] > Eb[-1]:
            Eb_plot = Eb[::-1]
            Ib_plot = Ib[::-1]
        else:
            Eb_plot = Eb
            Ib_plot = Ib
    
        # ---------- compute derivatives ----------
        Ef1, d1f = _compute_d1(Ef, If)
        Eb1, d1b = _compute_d1(Eb_plot, Ib_plot)
    
        Ef2, d2f = _compute_d2(Ef, If)
        Eb2, d2b = _compute_d2(Eb_plot, Ib_plot)


        if getattr(self, "var_zero_d2", None) is not None and self.var_zero_d2.get():
            if d2f is not None and len(d2f) > 0:
                d2f = d2f - d2f[0]
            if d2b is not None and len(d2b) > 0:
                d2b = d2b - d2b[0]    
        # ---------- integrate 2nd derivative (cumulative) ----------
        Af = _cumulative_trapz(d2f, Ef2) if (d2f is not None and Ef2 is not None and len(d2f)>1 and len(Ef2)>1) else None
        Ab = _cumulative_trapz(d2b, Eb2) if (d2b is not None and Eb2 is not None and len(d2b)>1 and len(Eb2)>1) else None
        
        if getattr(self, "var_zero_cum", None) is not None and self.var_zero_cum.get():
            if Af is not None and len(Af) > 0:
                Af = Af - Af[0]
            if Ab is not None and len(Ab) > 0:
                Ab = Ab - Ab[0]
    
        # ---------- significance detection on integral ----------
        k_sigma, min_run = _get_k_sigma_and_min_run()
        j_f = _detect_positive_rise(Ef2, Af, k_sigma=k_sigma, min_run=min_run) if (Af is not None) else None
        j_b = _detect_positive_rise(Eb2, Ab, k_sigma=k_sigma, min_run=min_run) if (Ab is not None) else None
    
        # ---------- draw ----------
        ax1, ax2, ax3 = self.deriv_axes
    
        for ax in (ax1, ax2, ax3):
            _clear_twins(ax)
            ax.clear()
            ax.grid(True)
    
        # TOP: 1st derivative only
        ax1.plot(Ef1, d1f, label="dI/dE (forward)")
        ax1.plot(Eb1, d1b, label="dI/dE (backward)")
        ax1.set_title("1st derivative")
        ax1.set_xlabel("E vs RHE (V)")
        ax1.set_ylabel("dI/dE (arb.)")
        ax1.legend(loc='upper left', fontsize=9)
    
        # MIDDLE: 2nd derivative + integral on twin axis
        ax2.plot(Ef2, d2f, label="d²I/dE² (forward)")
        ax2.set_title("2nd derivative — forward scan")
        ax2.set_xlabel("E vs RHE (V)")
        ax2.set_ylabel("d²I/dE² (arb.)")
        ax2.legend(loc='upper left', fontsize=9)
    
        ax2r = ax2.twinx()
        ax2r.plot(Ef2, Af, alpha=0.9, linestyle='--', label="∫ d²I/dE² dE (forward)")
        ax2r.set_ylabel("Accumulated (arb.)")
    
        # mark detected onset on integral
        # mark detected onset on integral
        if j_f is not None and 0 <= j_f < len(Ef2):
            x0 = float(Ef2[j_f])
            ax2.axvline(x0, linestyle=':', linewidth=1.6, alpha=0.9)
        
            # place text at 98% of axis height, independent of data scaling
            ax2.text(
                x0, 0.98, f" rise @ {x0:.3f} V",
                transform=ax2.get_xaxis_transform(),  # x in data, y in axes fraction
                rotation=90, va='top', ha='right',
                fontsize=9, alpha=0.9, clip_on=True
            )

    
        # align zeros (requires your helper to exist in CVApp)
        if hasattr(self, "_align_zeros"):
            self._align_zeros(ax2, ax2r)
    
        # BOTTOM: 2nd derivative + integral on twin axis
        ax3.plot(Eb2, d2b, label="d²I/dE² (backward)")
        ax3.set_title("2nd derivative — backward scan")
        ax3.set_xlabel("E vs RHE (V)")
        ax3.set_ylabel("d²I/dE² (arb.)")
        ax3.legend(loc='upper left', fontsize=9)
    
        ax3r = ax3.twinx()
        ax3r.plot(Eb2, Ab, alpha=0.9, linestyle='--', label="∫ d²I/dE² dE (backward)")
        ax3r.set_ylabel("Accumulated (arb.)")
    
        if j_b is not None and 0 <= j_b < len(Eb2):
            x0 = float(Eb2[j_b])
            ax3.axvline(x0, linestyle=':', linewidth=1.6, alpha=0.9)
            ax3.text(x0, ax3.get_ylim()[1], f" rise @ {x0:.3f} V",
                     rotation=90, va='top', ha='right', fontsize=9, alpha=0.9)
    
        if hasattr(self, "_align_zeros"):
            self._align_zeros(ax3, ax3r)
    
        # -------- restore limits (zoom) on left axes only --------
        for ax, (xyl) in zip((ax1, ax2, ax3), prev_limits):
            xlim, ylim = xyl
            if xlim is not None and ylim is not None:
                try:
                    ax.set_xlim(xlim)
                    ax.set_ylim(ylim)
                except Exception:
                    pass
    
        self.deriv_canvas.draw()
        self._deriv_last_pair = (p, cyc)
            


    # =================== Utilities ===================
    def get_display_name(self, item, default):
        return self.display_names.get(item, default)

    def _insert_dir_node(self, parent_id, dir_path):
        if dir_path in self.dir_nodes:
            return self.dir_nodes[dir_path]
        node_id = self.tree.insert(parent_id, 'end', text=os.path.basename(dir_path) or dir_path,
                                   values=('', 'dir', '', dir_path))
        self.checked[node_id] = False
        self.dir_nodes[dir_path] = node_id
        return node_id

    def _add_recent(self, folder):
        if folder in self.recents:
            self.recents.remove(folder)
        self.recents.insert(0, folder)
        self.recents = self.recents[:8]
        self.recent_menu.delete(0, 'end')
        for f in self.recents:
            self.recent_menu.add_command(label=f, command=lambda p=f: self.open_folder(path_override=p))

    # =================== Folder / File loading ===================
    def open_folder(self, path_override=None):
        initial = DEFAULT_FOLDER if os.path.isdir(DEFAULT_FOLDER) else None
        if self.tutorial_folder and os.path.isdir(self.tutorial_folder):
            initial = self.tutorial_folder
        if path_override is None:
            folder = filedialog.askdirectory(initialdir=initial) if initial else filedialog.askdirectory()
        else:
            folder = path_override
        if not folder:
            return

        self._add_recent(folder)
        if folder not in self.open_roots:
            self.open_roots.append(folder)
        self._insert_dir_node('', folder)

        for current_root, dirs, files in os.walk(folder):
            parent_id = self._insert_dir_node('', current_root) if current_root not in self.dir_nodes else self.dir_nodes[current_root]
            for d in sorted(dirs):
                self._insert_dir_node(parent_id, os.path.join(current_root, d))
            for fn in sorted(files):
                fpath = os.path.join(current_root, fn)
                if not is_cv_data_file(fpath):
                    continue

                # Lazy: only list the file name here. The file itself is
                # parsed on demand (first expand or first check) by
                # _ensure_file_loaded, so opening a folder never reads data
                # for files the user never looks at.
                file_id = self.tree.insert(parent_id, 'end', text=fn, values=('', 'file', '', fpath))
                self.file_nodes[file_id] = fpath
                self.checked[file_id] = False
                self.tree.insert(file_id, 'end', text='(loading…)', values=('', 'placeholder', '', fpath))

    def add_file(self):
        paths = filedialog.askopenfilenames(filetypes=DATA_FILETYPES)
        if not paths:
            return
        manual_root = self.dir_nodes.get('__MANUAL__')
        if manual_root is None:
            manual_root = self.tree.insert('', 'end', text='Manual Files', values=('', 'dir', '', '__MANUAL__'))
            self.dir_nodes['__MANUAL__'] = manual_root
            self.checked[manual_root] = False

        for p in paths:
            try:
                _, _, E, I, t, cyc, meta = read_mpt(p)
            except Exception as e:
                messagebox.showerror("Error", str(e))
                continue

            if p not in self.files:
                self.files.append(p)
            self.cache[p] = {'E':E,'I':I,'t':t,'cyc':cyc,'ocv':meta.get('ocv', False)}
            if p not in self.file_colors:
                self.file_colors[p] = self.color_cycle[len(self.file_colors) % len(self.color_cycle)]

            file_id = self.tree.insert(manual_root, 'end', text=os.path.basename(p), values=('', 'file', '', p))
            self.file_nodes[file_id] = p
            self.checked[file_id] = False
            for c in sorted(set(cyc)):
                node = self.tree.insert(file_id, 'end', text=f"Cycle {c}", values=('', 'cycle', '■', p))
                self.checked[node] = False
                sr = self._compute_cycle_scan_rate(p, c)
                self.scan_rates[(p, c)] = sr
                self._update_cycle_text_with_sr(node, c, sr)

    # =================== Clear ===================
    def clear_all(self):
        for it in self.tree.get_children():
            self.tree.delete(it)
        self.files.clear()
        self.cache.clear()
        self.file_nodes.clear()
        self.checked.clear()
        self.file_colors.clear()
        self.display_names.clear()
        self.dir_nodes.clear()
        self.open_roots.clear()
        self.peaks.clear()
        self.peak_index_map.clear()
        self.scan_rates.clear()
        self.analysis_rows.clear()
        self.analysis_map.clear()
        self.ax.clear()
        self.canvas.draw()

        self.peak_list.delete(0, 'end')
        self.debug_text.delete('1.0', 'end')
        self.tafel_text.delete('1.0', 'end')

        self.analysis_text.delete('1.0', 'end')
        for it in self.analysis_tree.get_children():
            self.analysis_tree.delete(it)
        for _axname, _cvname in (('analysis_ax', 'analysis_canvas'),
                                 ('cmp_ax', 'cmp_canvas')):
            if hasattr(self, _axname) and hasattr(self, _cvname):
                getattr(self, _axname).clear()
                getattr(self, _cvname).draw()
        self.last_debug = {}
        self._deriv_last_pair = None
        if hasattr(self, 'deriv_axes') and hasattr(self, 'deriv_canvas'):
            for ax in self.deriv_axes:
                ax.clear()
            self.deriv_canvas.draw()




    # =================== Helpers: cycle scan rate & text ===================
    def _compute_cycle_scan_rate(self, path, cycle):
        try:
            E = self.cache[path]['E']; t = self.cache[path]['t']; cyc_arr = self.cache[path]['cyc']
            sel = cyc_arr == cycle
            if not np.any(sel):
                return np.nan
            return estimate_scan_rate(E[sel], t[sel])
        except Exception:
            return np.nan

    def _update_cycle_text_with_sr(self, node_id, cycle, sr_Vs):
        try:
            if np.isfinite(sr_Vs):
                sr_mVs = sr_Vs * 1000.0
                self.tree.item(node_id, text=f"Cycle {cycle} ({sr_mVs:.2f} mV/s)")
            else:
                self.tree.item(node_id, text=f"Cycle {cycle}")
        except Exception:
            pass

    # =================== Left panel collapse ===================
    def toggle_left_panel(self):
        """Hide/show the file tree so the plots get the full width on a
        small screen. The toggle button itself always stays visible."""
        try:
            if self._left_visible:
                try:
                    self.main_paned.paneconfigure(self.left_panel, hide=True)
                except Exception:
                    self.main_paned.forget(self.left_panel)
                self.collapse_btn.config(text='▶')
                self._left_visible = False
            else:
                try:
                    self.main_paned.paneconfigure(self.left_panel, hide=False)
                except Exception:
                    self.main_paned.add(self.left_panel, before=self.center_panel, minsize=220)
                self.collapse_btn.config(text='◀')
                self._left_visible = True
        except Exception:
            pass

    # =================== Lazy file loading ===================
    def _ensure_file_loaded(self, fpath, file_id):
        """Parse fpath (once) and replace its '(loading…)' placeholder with
        real Cycle nodes. Safe to call repeatedly - a no-op once cached."""
        if fpath in self.cache:
            return True
        try:
            _, _, E, I, t, cyc, meta = read_mpt(fpath)
        except Exception:
            return False  # leave the placeholder; nothing usable to show

        if fpath not in self.files:
            self.files.append(fpath)
        self.cache[fpath] = {'E': E, 'I': I, 't': t, 'cyc': cyc, 'ocv': meta.get('ocv', False)}
        if fpath not in self.file_colors:
            self.file_colors[fpath] = self.color_cycle[len(self.file_colors) % len(self.color_cycle)]

        for ch in self.tree.get_children(file_id):
            self.tree.delete(ch)

        for c in sorted(set(cyc)):
            node = self.tree.insert(file_id, 'end', text=f"Cycle {c}", values=('', 'cycle', '■', fpath))
            self.checked[node] = False
            sr = self._compute_cycle_scan_rate(fpath, c)
            self.scan_rates[(fpath, c)] = sr
            self._update_cycle_text_with_sr(node, c, sr)
        return True

    def _ensure_dir_loaded(self, dir_id):
        """Recursively load every file listed under a directory node
        (used when the user checks a whole folder at once)."""
        for ch in self.tree.get_children(dir_id):
            fpath = self.file_nodes.get(ch)
            if fpath:
                self._ensure_file_loaded(fpath, ch)
            else:
                self._ensure_dir_loaded(ch)

    def on_tree_expand(self, event):
        """Lazy-load a file's cycles the first time its node is expanded."""
        item = self.tree.focus()
        fpath = self.file_nodes.get(item)
        if fpath:
            self._ensure_file_loaded(fpath, item)

    # =================== Unified Tree Click ===================
    def on_tree_click(self, event):
        col = self.tree.identify_column(event.x)
        item = self.tree.identify_row(event.y)
        if not item:
            return

        # Checkbox column
        if col == '#1':
            vals0 = list(self.tree.item(item, 'values'))
            vals0 = (vals0 + ['']*(4-len(vals0)))
            if vals0[1] == 'placeholder':
                return  # ignore clicks on the "(loading…)" row itself

            new = not self.checked.get(item, False)
            self.checked[item] = new
            vals = list(self.tree.item(item, 'values'))
            vals = (vals + ['']*(4-len(vals)))
            vals[0] = '✔' if new else ''
            self.tree.item(item, values=vals)

            # If turning a file/dir ON, make sure its data is actually
            # parsed first, so the propagation below finds real cycle
            # nodes instead of the placeholder.
            if new:
                ntype = vals[1] if len(vals) > 1 else None
                if ntype == 'file':
                    fpath = self.file_nodes.get(item)
                    if fpath:
                        self._ensure_file_loaded(fpath, item)
                elif ntype == 'dir':
                    self._ensure_dir_loaded(item)

            # Propagate to children
            for ch in self.tree.get_children(item):
                self.checked[ch] = new
                cv = list(self.tree.item(ch, 'values'))
                cv = (cv + ['']*(4-len(cv)))
                cv[0] = '✔' if new else ''
                self.tree.item(ch, values=cv)
                for gch in self.tree.get_children(ch):
                    self.checked[gch] = new
                    gcv = list(self.tree.item(gch, 'values'))
                    gcv = (gcv + ['']*(4-len(gcv)))
                    gcv[0] = '✔' if new else ''
                    self.tree.item(gch, values=gcv)

            # Update parent check state
            parent = self.tree.parent(item)
            if parent:
                all_on = all(self.checked.get(ch, False) for ch in self.tree.get_children(parent))
                pv = list(self.tree.item(parent, 'values')); pv = (pv + ['']*(4-len(pv)))
                pv[0] = '✔' if all_on else ''
                self.tree.item(parent, values=pv)

            gp = self.tree.parent(parent)
            if gp:
                all_on = all(self.checked.get(ch, False) for ch in self.tree.get_children(gp))
                ppv = list(self.tree.item(gp, 'values')); ppv = (ppv + ['']*(4-len(ppv)))
                ppv[0] = '✔' if all_on else ''
                self.tree.item(gp, values=ppv)

            self.plot()
            return

        # Color column
        if col == '#3':
            vals = self.tree.item(item, 'values')
            if len(vals) < 2 or vals[1] != 'cycle':
                return

            parent = self.tree.parent(item)
            path = self.tree.item(parent, 'values')[3] or self.file_nodes.get(parent)
            if not path:
                return

            try:
                label = self.tree.item(item, 'text')
                cyc = int(label.split()[1])
            except Exception:
                return

            chosen = askcolor(title="Choose cycle color")
            if not chosen or not chosen[1]:
                return
            chosen_hex = chosen[1]

            self.file_colors[(path, cyc)] = chosen_hex

            new_vals = list(vals)
            new_vals = (new_vals + ['']*(4-len(new_vals)))
            new_vals[2] = '■'
            tagname = f"color_{path}_{cyc}"
            self.tree.item(item, values=new_vals, tags=(tagname,))
            self.tree.tag_configure(tagname, foreground=chosen_hex)

            self.plot()
            return

    # =================== Rename ===================
    def rename_item(self, event):
        item = self.tree.focus()
        if not item:
            return
        old = self.tree.item(item, 'text')
        new = simpledialog.askstring("Rename", f"Rename '{old}' to:")
        if new:
            self.tree.item(item, text=new)
            self.display_names[item] = new
            self.plot()

    # =================== Build series (optional OER) ===================
    def get_display_series(self, p, cyc, pH):
        try:
            E = self.cache[p]['E']; I = self.cache[p]['I']; t = self.cache[p]['t']
            cyc_arr = self.cache[p]['cyc']; is_ocv = self.cache[p].get('ocv', False)
        except Exception:
            return None, None, None, None, None, None, None, {}, False

        sel = (cyc_arr == cyc)
        if not np.any(sel):
            return None, None, None, None, None, None, None, {}, is_ocv

        E_rhe = convert_reference(E[sel], self.ref_var.get(), pH)
        I_sel = I[sel]
        t_sel = t[sel]

        mask_corr = np.zeros_like(I_sel, dtype=bool)
        I_corr = I_sel.copy()
        oer_dbg = {'enabled': False}

        if (not is_ocv) and self.subtract_oer.get():
            I_blend, mask_corr, dbg = oer_subtract_modeA(E_rhe, I_sel)
            I_corr = I_blend
            oer_dbg = {'enabled': True, **dbg}

        mode = self.plot_mode.get()
        if mode in ('IvsE', 'QvsE'):
            if mode == 'QvsE':
                Y_raw = integrate_Q(I_sel, t_sel)
                Y_corr = integrate_Q(I_corr, t_sel)
            else:
                Y_raw = I_sel
                Y_corr = I_corr
        elif mode == 'IvsT':
            Y_raw = I_sel
            Y_corr = I_corr
        elif mode == 'EvsT':
            Y_raw = E_rhe
            Y_corr = E_rhe
        else:
            Y_raw = I_sel
            Y_corr = I_corr

        I_for_peaks = I_corr
        dbg_total = {
            'E_len': len(E_rhe), 'I_raw_len': len(I_sel), 'I_corr_len': len(I_corr),
            't_len': len(t_sel), 'mask_corr_applied': int(np.sum(mask_corr)),
            'oer': oer_dbg, 'is_ocv': is_ocv
        }
        return E_rhe, t_sel, Y_raw, Y_corr, mask_corr, I_for_peaks, I_sel, dbg_total, is_ocv

    # =================== Auto peaks ===================
    def run_auto_peaks(self):
        try:
            thresh = float(self.min_height_var.get())
            cath_low = float(self.cath_low.get()); cath_high = float(self.cath_high.get())
            an_low = float(self.an_low.get()); an_high = float(self.an_high.get())
        except Exception:
            messagebox.showerror("Error", "Invalid peak detection values.")
            return

        self.peaks.clear()

        for it, ch in self.checked.items():
            if not ch or self.tree.item(it, 'values')[1] != 'cycle':
                continue

            parent = self.tree.parent(it)
            p = self.tree.item(parent, 'values')[3] or self.file_nodes.get(parent)
            if not p:
                continue
            try:
                label = self.tree.item(it, 'text')
                cyc = int(label.split()[1])
            except Exception:
                continue

            try:
                pH = float(self.ph_var.get())
            except Exception:
                pH = 14.0

            E_rhe, t_sel, Y_raw, Y_corr, mask_corr, I_for_peaks, I_sel, dbg, is_ocv = self.get_display_series(p, cyc, pH)
            if E_rhe is None or is_ocv:
                continue

            x_peaks = E_rhe
            y_peaks = I_for_peaks
            if self.sanity_mode.get():
                x_peaks, y_peaks, _, _ = clean_xy(x_peaks, y_peaks)

            # OPTIONAL: smooth current for peak detection
            y_detect = self._maybe_smooth(y_peaks, for_what='peaks')

            found = []
            for i in range(1, len(y_detect) - 1):
                Ei = x_peaks[i]
                Yi = y_detect[i]
                if not np.isfinite(Yi) or abs(Yi) < thresh:
                    continue
                in_anod = (an_low <= Ei <= an_high)
                in_cath = (cath_low <= Ei <= cath_high)
                is_max = (y_detect[i] > y_detect[i-1] and y_detect[i] > y_detect[i+1])
                is_min = (y_detect[i] < y_detect[i-1] and y_detect[i] < y_detect[i+1])
                if in_anod and is_max:
                    found.append({'E': float(x_peaks[i]), 'Y': float(y_peaks[i]), 'type': 'anodic'})
                elif in_cath and is_min:
                    found.append({'E': float(x_peaks[i]), 'Y': float(y_peaks[i]), 'type': 'cathodic'})
            self.peaks[(p, cyc)] = found

        self.refresh_peak_panel()
        self.plot()

    # =================== Modes ===================
    def enable_manual_peak_mode(self):
        self.manual_peak_mode = True
        self.delete_peak_mode = False
        messagebox.showinfo("Manual Peak Mode", "Click on curve to add a peak (within defined windows). For OCV files, peaks are disabled.")

    def enable_delete_peak_mode(self):
        self.delete_peak_mode = True
        self.manual_peak_mode = False
        messagebox.showinfo("Delete Peak Mode", "Click on a labeled peak to delete it.")

    # =================== Plot click ===================
    def on_click_plot(self, event):
        if event.inaxes != self.ax:
            return
        if self.delete_peak_mode:
            self.delete_peak_by_click(event.xdata, event.ydata)
            return
        if not self.manual_peak_mode:
            return

        target = None
        for it, ch in self.checked.items():
            if ch and self.tree.item(it, 'values')[1] == 'cycle':
                parent = self.tree.parent(it)
                p = self.tree.item(parent, 'values')[3] or self.file_nodes.get(parent)
                if not p:
                    continue
                try:
                    label = self.tree.item(it, 'text')
                    cyc = int(label.split()[1])
                except Exception:
                    continue
                target = (p, cyc)
                break
        if target is None:
            messagebox.showerror("No cycle selected", "Check a cycle in the tree.")
            return

        p, cyc = target
        if self.cache.get(p, {}).get('ocv', False):
            messagebox.showwarning("OCV file", "Peaks are not applicable for OCV files (no current).")
            return

        try:
            pH = float(self.ph_var.get())
        except Exception:
            pH = 14.0

        E_rhe, t_sel, Y_raw, Y_corr, mask_corr, I_for_peaks, I_sel, dbg, is_ocv = self.get_display_series(p, cyc, pH)
        if E_rhe is None or is_ocv:
            return

        mode = self.plot_mode.get()
        if mode in ('IvsE', 'QvsE'):
            xarr = E_rhe
            yarr = I_for_peaks
        elif mode == 'IvsT':
            xarr = t_sel - (t_sel[0] if len(t_sel) > 0 and self.time_zero.get() else 0)
            yarr = I_for_peaks
        elif mode == 'EvsT':
            xarr = t_sel - (t_sel[0] if len(t_sel) > 0 and self.time_zero.get() else 0)
            yarr = E_rhe
        else:
            xarr = E_rhe
            yarr = I_for_peaks

        if self.sanity_mode.get():
            xarr, yarr, _, _ = clean_xy(xarr, yarr)
        if len(xarr) == 0:
            return

        idx = int(np.argmin(np.abs(xarr - event.xdata)))
        Ep = float(xarr[idx]); Yp = float(yarr[idx])

        try:
            cath_low = float(self.cath_low.get()); cath_high = float(self.cath_high.get())
            an_low = float(self.an_low.get()); an_high = float(self.an_high.get())
        except Exception:
            cath_low, cath_high, an_low, an_high = 1.25, 1.40, 1.35, 1.50

        in_cath = (cath_low <= Ep <= cath_high) if mode in ('IvsE', 'QvsE') else False
        in_anod = (an_low <= Ep <= an_high) if mode in ('IvsE', 'QvsE') else False

        typ = 'cathodic' if in_cath else ('anodic' if in_anod else 'manual')
        if (p, cyc) not in self.peaks:
            self.peaks[(p, cyc)] = []
        self.peaks[(p, cyc)].append({'E': Ep, 'Y': Yp, 'type': typ})
        self.refresh_peak_panel()
        self.plot()

    # =================== Delete peak by click ===================
    def delete_peak_by_click(self, x, y, tol=0.02):
        nearest = None; best_d = None
        for key, lst in self.peaks.items():
            for pk in lst:
                d = ((pk['E']-x)**2 + (pk['Y']-y)**2)**0.5
                if best_d is None or d < best_d:
                    best_d = d; nearest = (key, pk)
        if nearest and best_d <= tol:
            key, pk = nearest
            self.peaks[key].remove(pk)
            if not self.peaks[key]:
                del self.peaks[key]
            self.refresh_peak_panel()
            self.plot()

    # =================== Peak list panel ===================
    def refresh_peak_panel(self):
        self.peak_list.delete(0, 'end')
        self.peak_index_map.clear()
        for (p, cyc), lst in self.peaks.items():
            for pk in lst:
                self.peak_list.insert(
                    'end',
                    f"{os.path.basename(p)} | Cycle {cyc} | {pk['type']:8s} | E={pk['E']:.3f} | I={pk['Y']:.3f}"
                )
                self.peak_index_map.append((p, cyc, pk['E'], pk['Y']))

    def delete_selected_peak(self):
        sel = list(self.peak_list.curselection())
        if not sel:
            messagebox.showwarning("No Selection", "Please select peaks to delete.")
            return
        for idx in sorted(sel, reverse=True):
            if idx >= len(self.peak_index_map):
                continue
            p, cyc, Ee, Yy = self.peak_index_map[idx]
            key = (p, cyc)
            if key in self.peaks:
                new_list = []
                removed = False
                for pk in self.peaks[key]:
                    if (not removed) and abs(pk['E'] - Ee) < 1e-9 and abs(pk['Y'] - Yy) < 1e-9:
                        removed = True
                        continue
                    new_list.append(pk)
                if removed:
                    if new_list:
                        self.peaks[key] = new_list
                    else:
                        del self.peaks[key]
            self.peak_list.delete(idx)
            del self.peak_index_map[idx]
        self.plot()

    # =================== Plot ===================
    def plot(self):
        self.ax.clear()
        try:
            pH = float(self.ph_var.get())
        except Exception:
            pH = 14.0
        try:
            fsz = int(float(self.font_var.get()))
        except Exception:
            fsz = 12
        try:
            lw = float(self.line_var.get())
        except Exception:
            lw = 2.0
    
        mode = self.plot_mode.get()
        self.last_debug = {}
        any_plot = False
    
        # draw onset markers only once per plot refresh
        onset_drawn = False
    
        # helper: compute onset E percentiles on forward scan (robust)
        def _compute_onset_forward(E_arr, I_arr, e_low, e_high, percentiles=(5,10,15,25)):
            E_arr = np.asarray(E_arr, dtype=float)
            I_arr = np.asarray(I_arr, dtype=float)
            out = {f'onset_E_p{p}': np.nan for p in percentiles}
            if len(E_arr) < 6:
                return out
    
            # forward scan = indices up to max E
            idx_sw = int(np.argmax(E_arr))
            fmask = np.arange(len(E_arr)) <= idx_sw
    
            Ef = E_arr[fmask]
            If = I_arr[fmask]
    
            m = (Ef >= e_low) & (Ef <= e_high) & np.isfinite(Ef) & np.isfinite(If)
            if np.sum(m) < 6:
                return out
    
            Ew = Ef[m]
            Iw = If[m]
    
            # sort by E so "first E crossing" is well-defined
            order = np.argsort(Ew)
            Ew = Ew[order]
            Iw = Iw[order]
    
            # robust thresholds: fraction of dynamic range (prevents identical values)
            Imin = float(np.nanmin(Iw))
            Imax = float(np.nanmax(Iw))
            if (not np.isfinite(Imin)) or (not np.isfinite(Imax)) or (Imax <= Imin):
                return out
    
            for p in percentiles:
                thr = Imin + (p / 100.0) * (Imax - Imin)
                idxs = np.where(Iw >= thr)[0]
                if len(idxs):
                    out[f'onset_E_p{p}'] = float(Ew[idxs[0]])
            return out
    
        for it, ch in self.checked.items():
            if not ch or self.tree.item(it, 'values')[1] != 'cycle':
                continue
    
            parent = self.tree.parent(it)
            p = self.tree.item(parent, 'values')[3] or self.file_nodes.get(parent)
            if not p:
                continue
            try:
                label = self.tree.item(it, 'text')
                cyc = int(label.split()[1])
            except Exception:
                continue
    
            E_rhe, t_sel, Y_raw, Y_corr, mask_corr, I_for_peaks, I_sel, dbg, is_ocv = self.get_display_series(p, cyc, pH)
            if E_rhe is None:
                continue
    
            if mode in ('IvsT', 'EvsT') and len(t_sel) > 0 and self.time_zero.get():
                t_plot = t_sel - t_sel[0]
            else:
                t_plot = t_sel
    
            removed_main = 0
            if mode in ('IvsE', 'QvsE'):
                Xr, Yr = (E_rhe, Y_raw)
                Xc, Yc = (E_rhe, Y_corr)
                if self.sanity_mode.get():
                    Xr, Yr, _, rem1 = clean_xy(Xr, Yr)
                    Xc, Yc, _, rem2 = clean_xy(Xc, Yc)
                    removed_main = rem1 + rem2
            elif mode == 'IvsT':
                Xr, Yr = (t_plot, Y_raw)
                Xc, Yc = (t_plot, Y_corr)
                if self.sanity_mode.get():
                    Xr, Yr, _, rem1 = clean_xy(Xr, Yr)
                    Xc, Yc, _, rem2 = clean_xy(Xc, Yc)
                    removed_main = rem1 + rem2
            elif mode == 'EvsT':
                Xr, Yr = (t_plot, E_rhe)
                Xc, Yc = (t_plot, E_rhe)
                if self.sanity_mode.get():
                    Xr, Yr, _, rem1 = clean_xy(Xr, Yr)
                    Xc, Yc = Xr, Yr
                    removed_main = rem1
            else:
                Xr, Yr = (E_rhe, Y_raw)
                Xc, Yc = (E_rhe, Y_corr)
    
            file_name = self.get_display_name(parent, os.path.basename(p))
            sr = self.scan_rates.get((p, cyc), np.nan)
            if np.isfinite(sr):
                cycle_name = self.get_display_name(it, f"Cycle {cyc} ({sr*1000:.2f} mV/s)")
                legend_suffix = f", {sr*1000:.2f} mV/s"
            else:
                cycle_name = self.get_display_name(it, f"Cycle {cyc}")
                legend_suffix = ""
    
            color = self.file_colors.get((p, cyc), self.file_colors.get(p, 'C0'))
    
            # RAW curve
            if len(Xr) and len(Yr):
                self.ax.plot(
                    Xr, Yr,
                    color=color, linewidth=lw,
                    label=f"{file_name} – {cycle_name} (raw{legend_suffix}{', OCV' if is_ocv else ''})"
                )
    
            # Corrected curve: only if checkbox ON and not OCV (same as your current behavior)
            if self.subtract_oer.get() and (not is_ocv):
                if len(Xc) and len(Yc):
                    self.ax.plot(
                        Xc, Yc,
                        color=color, linestyle='--', linewidth=lw+1,
                        label=f"{file_name} – {cycle_name} (corrected{legend_suffix})"
                    )
    
            # ---- Onset percentile markers on RAW plot (only once) ----
            if (not onset_drawn) and (mode == 'IvsE') and (not is_ocv) and getattr(self, "show_onset_lines", None) is not None and self.show_onset_lines.get():
                try:
                    e_an_low = float(self.an_low.get()); e_an_high = float(self.an_high.get())
                except Exception:
                    e_an_low, e_an_high = 1.35, 1.50
    
                # compute on forward scan only
                onset = _compute_onset_forward(Xr, Yr, e_an_low, e_an_high, percentiles=(5,10,15,25))
    
                # draw vertical lines + labels (no legend spam)
                y0, y1 = self.ax.get_ylim()
                y_text = y1 - 0.04*(y1 - y0)
    
                for pp in (5,10,15,25):
                    xline = onset.get(f'onset_E_p{pp}', np.nan)
                    if np.isfinite(xline):
                        self.ax.axvline(xline, linestyle=':', linewidth=1.2, alpha=0.7, label='_nolegend_')
                        self.ax.text(
                            xline, y_text, f"{pp}%",
                            rotation=90, va='top', ha='center',
                            fontsize=max(8, int(0.8*fsz)), alpha=0.85
                        )
    
                onset_drawn = True
                
            # ---- Onset kσ markers on RAW plot (only once) ----
            if (not onset_drawn) and (mode == 'IvsE') and (not is_ocv) and self.show_onset_k_lines.get():
                try:
                    e_an_low = float(self.an_low.get()); e_an_high = float(self.an_high.get())
                except Exception:
                    e_an_low, e_an_high = 1.35, 1.50
            
                ks = (2,3,4,5,8,12,15,20)
                onsetk = self._onset_k_sigma_integrator(
                    Xr, Yr, e_an_low, e_an_high,
                    ks=ks, forward_only=True, min_run=int(self.onset_min_run.get())
                )
            
                y0, y1 = self.ax.get_ylim()
                y_text = y1 - 0.04*(y1 - y0)
            
                for k in ks:
                    xline = onsetk.get(f'onset_E_k{k}', np.nan)
                    if np.isfinite(xline):
                        self.ax.axvline(xline, linestyle='--', linewidth=1.2, alpha=0.7, label='_nolegend_')
                        self.ax.text(
                            xline, y_text, f"k={k}",
                            rotation=90, va='top', ha='center',
                            fontsize=max(8, int(0.8*fsz)), alpha=0.85
                        )
            
                onset_drawn = True

    
            # Peaks overlay (IvsE/QvsE only)
            if (not is_ocv) and mode in ('IvsE', 'QvsE'):
                key = (p, cyc)
                for pk in self.peaks.get(key, []):
                    col = 'red' if pk['type'] == 'anodic' else ('blue' if pk['type'] == 'cathodic' else color)
                    self.ax.plot(pk['E'], pk['Y'], 'o', color=col, markersize=7)
                    self.ax.text(
                        pk['E'], pk['Y'],
                        f"E={pk['E']:.3f} V\nI={pk['Y']:.3f} mA",
                        color=col, fontsize=max(fsz, 12)
                    )
    
            any_plot = True
            self.last_debug[f"{os.path.basename(p)} | Cycle {cyc}"] = {
                **dbg,
                'mode': mode,
                'removed_points_plot': int(removed_main)
            }
    
        # axes labels, legend, etc.
        if any_plot:
            if mode == 'IvsE':
                self.ax.set_xlabel("E vs RHE (V)", fontsize=fsz)
                self.ax.set_ylabel("Current (mA)", fontsize=fsz)
            elif mode == 'QvsE':
                self.ax.set_xlabel("E vs RHE (V)", fontsize=fsz)
                self.ax.set_ylabel("Charge (mA·s)", fontsize=fsz)
            elif mode == 'IvsT':
                self.ax.set_xlabel("Time (s)", fontsize=fsz)
                self.ax.set_ylabel("Current (mA)", fontsize=fsz)
            elif mode == 'EvsT':
                self.ax.set_xlabel("Time (s)", fontsize=fsz)
                self.ax.set_ylabel("E vs RHE (V)", fontsize=fsz)
    
            self.ax.tick_params(labelsize=fsz)
            self.ax.grid(True)
    
            if not self.hide_legend.get():
                self.ax.legend(fontsize=fsz)
    
        self.canvas.draw()
        self.refresh_debug_panel()


    # =================== Debug panel ===================
    def refresh_debug_panel(self):
        self.debug_text.delete('1.0', 'end')
        if not self.last_debug:
            self.debug_text.insert('end', "No debug yet. Plot a selection to populate diagnostics.\n")
            return
        for name, info in self.last_debug.items():
            self.debug_text.insert('end', f"[{name}]\n")
            self.debug_text.insert('end', f"  Mode: {info.get('mode')}  |  OCV: {info.get('is_ocv')}\n")
            self.debug_text.insert('end', f"  E_len: {info.get('E_len')} | I_raw_len: {info.get('I_raw_len')} | I_corr_len: {info.get('I_corr_len')} | t_len: {info.get('t_len')}\n")
            self.debug_text.insert('end', f"  mask_corr_applied: {info.get('mask_corr_applied')}\n")
            oer = info.get('oer', {})
            self.debug_text.insert('end', f"  OER: enabled={oer.get('enabled')} | valid={oer.get('valid')} | fw_pos={oer.get('fw_pos')} | bw_pos={oer.get('bw_pos')} | applied={oer.get('applied')}\n")
            self.debug_text.insert('end', f"  removed_points_plot: {info.get('removed_points_plot')}\n\n")

    # =================== Tafel ===================
    def build_and_plot_tafel(self):
        try:
            E0 = float(self.tafel_E0.get())
            N = int(float(self.tafel_N.get()))
            logI_min = float(self.tafel_logI_min.get())
            logI_max = float(self.tafel_logI_max.get())
        except Exception:
            messagebox.showerror("Error", "Invalid Tafel parameters.")
            return

        pts = []
        for it, ch in self.checked.items():
            if not ch or self.tree.item(it, 'values')[1] != 'cycle':
                continue
            parent = self.tree.parent(it)
            p = self.tree.item(parent, 'values')[3] or self.file_nodes.get(parent)
            if not p:
                continue
            try:
                label_txt = self.tree.item(it, 'text')
                cyc = int(label_txt.split()[1])
            except Exception:
                continue
            try:
                pH = float(self.ph_var.get())
            except Exception:
                pH = 14.0

            E_rhe, t_sel, Y_raw, Y_corr, mask_corr, I_for_peaks, I_sel, dbg, is_ocv = self.get_display_series(p, cyc, pH)
            if E_rhe is None or len(E_rhe) == 0 or len(I_sel) == 0:
                continue
            if is_ocv:
                continue

            n = len(E_rhe)
            tail = max(1, int(0.2 * n))
            E_set = float(np.median(E_rhe[-tail:]))

            N_eff = min(max(1, N), len(I_sel))
            Iss_A = float(np.mean(I_sel[-N_eff:])) * 1e-3  # mA -> A
            if not np.isfinite(Iss_A) or Iss_A == 0.0:
                continue

            eta = E_set - E0
            logI = np.log10(abs(Iss_A))
            color = self.file_colors.get((p, cyc), self.file_colors.get(p, 'C0'))
            name = f"{os.path.basename(p)} – Cycle {cyc}"
            pts.append((logI, eta, name, color))

        if len(pts) < 3:
            messagebox.showwarning("Tafel", "Need at least 3 cycles (points).")
            return

        pts_win = [q for q in pts if (logI_min <= q[0] <= logI_max)]
        if len(pts_win) < 2:
            messagebox.showwarning("Tafel", "Not enough points inside the chosen log(I) window.")
            pts_win = pts

        X = np.array([p[0] for p in pts_win])
        Y = np.array([p[1] for p in pts_win])

        A = np.vstack([X, np.ones_like(X)]).T
        try:
            sol, *_ = np.linalg.lstsq(A, Y, rcond=None)
            a, b = sol[0], sol[1]
            Yhat = a * X + b
            ss_res = float(np.sum((Y - Yhat) ** 2))
            ss_tot = float(np.sum((Y - np.mean(Y)) ** 2))
            R2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        except Exception:
            a = b = R2 = np.nan

        self.ax.clear()
        self.ax.grid(True)
        for lg, et, nm, col in pts:
            self.ax.plot(lg, et, 'o', label=nm, color=col)

        if np.isfinite(a):
            xs = np.linspace(min(X), max(X), 100)
            ys = a * xs + b
            self.ax.plot(xs, ys, '--', linewidth=2.2, label=f"Fit: slope={a*1000:.1f} mV/dec, R²={R2:.3f}")

        self.ax.set_xlabel("log10(|I_ss|) [A]", fontsize=int(float(self.font_var.get())))
        self.ax.set_ylabel("Overpotential η = Eset - E0 (V)", fontsize=int(float(self.font_var.get())))
        if not self.hide_legend.get():
            self.ax.legend(fontsize=int(float(self.font_var.get())))
        self.canvas.draw()

        self.tafel_text.delete('1.0', 'end')
        self.tafel_text.insert('end', f"Tafel points (η vs log10|I|), E0={E0:.3f} V, N={N} tail-avg\n")
        for lg, et, nm, _ in sorted(pts, key=lambda q: q[0]):
            self.tafel_text.insert('end', f"  {nm:<40s}  logI={lg:>7.3f}   η={et:>7.4f} V\n")
        self.tafel_text.insert('end', "\nFit window:\n")
        self.tafel_text.insert('end', f"  logI_min={logI_min}, logI_max={logI_max}\n")
        if np.isfinite(a):
            self.tafel_text.insert('end', f"\nResult:\n  slope = {a*1000:.2f} mV/dec\n  intercept = {b:.4f} V\n  R^2 = {R2:.4f}\n")
        else:
            self.tafel_text.insert('end', "\nResult: fit failed.\n")

    # =================== Analysis tab ===================
    def _build_analysis_tab(self, tab):
        ctrl = tk.Frame(tab); ctrl.pack(fill='x', padx=8, pady=6)
        self.btn_run_analysis = tk.Button(ctrl, text="Run Analysis on Checked Cycles", command=self.run_analysis)
        self.btn_run_analysis.pack(side='left')
        self.btn_export = tk.Button(ctrl, text="Export Analysis CSV", command=self.export_analysis_csv)
        self.btn_export.pack(side='left', padx=(6,0))
        # "i": what every column of the table means
        self.btn_param_info = tk.Button(ctrl, text=" i ", font=('Georgia', 10, 'bold italic'),
                                        fg='#0a4d8c', command=self.show_param_help)
        self.btn_param_info.pack(side='left', padx=(10,0))

        # Table on top, statistics report below (drag the sash to resize).
        split = tk.PanedWindow(tab, orient='vertical', sashrelief='raised', sashwidth=6)
        split.pack(fill='both', expand=True, padx=8, pady=(0,8))

        table_frame = tk.Frame(split)
        split.add(table_frame, minsize=140)

        # v8.12: onset split into RAW and CORR
        cols = [
            'file', 'cycle',
            'anodic_E_raw', 'anodic_I_raw',
            'cathodic_E_raw', 'cathodic_I_raw',
            'anodic_E_corr', 'anodic_I_corr',
            'cathodic_E_corr', 'cathodic_I_corr',
            'I_ratio_raw', 'I_ratio_corr',
            'dEp_raw', 'dEp_corr',
            'Q_end_raw', 'Q_end_corr',
            'onset_raw_E_p5', 'onset_raw_E_p10', 'onset_raw_E_p15', 'onset_raw_E_p25',
            'onset_corr_E_p5', 'onset_corr_E_p10', 'onset_corr_E_p15', 'onset_corr_E_p25',
            'onset_raw_E_k2','onset_raw_E_k3','onset_raw_E_k4','onset_raw_E_k5','onset_raw_E_k8','onset_raw_E_k12','onset_raw_E_k15','onset_raw_E_k20',
            'onset_corr_E_k2','onset_corr_E_k3','onset_corr_E_k4','onset_corr_E_k5','onset_corr_E_k8','onset_corr_E_k12','onset_corr_E_k15','onset_corr_E_k20',

        ]
        
        fixed_cols = []
        for Et in TARGET_E:
            k = _ekey(Et)
            fixed_cols += [f"raw_I_fwd_{k}", f"raw_I_bwd_{k}", f"raw_I_for_back_{k}"]
        
        fixed_cols += ["raw_junction_E", "raw_junction_I_fwd", "raw_junction_I_bwd", "raw_junction_for_back"]
        
        cols = list(cols) + fixed_cols
        cols = tuple(cols)


        self.analysis_tree = ttk.Treeview(table_frame, columns=cols, show='headings', height=8)
        for c in cols:
            self.analysis_tree.heading(c, text=c)
            self.analysis_tree.column(c, width=120, anchor='center', stretch=False)

        self.analysis_tree.column('file', width=240, anchor='w', stretch=False)
        self.analysis_tree.column('cycle', width=60, anchor='center', stretch=False)

        asb_y = tk.Scrollbar(table_frame, orient='vertical', command=self.analysis_tree.yview)
        asb_y.pack(side='right', fill='y')
        asb_x = tk.Scrollbar(table_frame, orient='horizontal', command=self.analysis_tree.xview)
        asb_x.pack(side='bottom', fill='x')
        self.analysis_tree.config(yscrollcommand=asb_y.set, xscrollcommand=asb_x.set)

        self.analysis_tree.pack(side='left', fill='both', expand=True)

        # ---- Report pane ----
        report_frame = tk.Frame(split)
        split.add(report_frame, minsize=160)

        rhead = tk.Frame(report_frame)
        rhead.pack(fill='x', pady=(4,2))
        tk.Label(rhead, text="Report", font=('TkDefaultFont', 9, 'bold')).pack(side='left')
        tk.Button(rhead, text="Copy report", command=self._copy_report).pack(side='right')

        rsb = tk.Scrollbar(report_frame, orient='vertical')
        rsb.pack(side='right', fill='y')
        rsb_x = tk.Scrollbar(report_frame, orient='horizontal')
        rsb_x.pack(side='bottom', fill='x')
        # Monospace so the statistics columns line up.
        self.analysis_text = tk.Text(report_frame, height=12, wrap='none',
                                     font=('Consolas', 9),
                                     yscrollcommand=rsb.set, xscrollcommand=rsb_x.set)
        self.analysis_text.pack(side='left', fill='both', expand=True)
        rsb.config(command=self.analysis_text.yview)
        rsb_x.config(command=self.analysis_text.xview)


    def _first_checked_cycle(self):
        for it, ch in self.checked.items():
            if ch and self.tree.item(it, 'values')[1] == 'cycle':
                parent = self.tree.parent(it)
                p = self.tree.item(parent, 'values')[3] or self.file_nodes.get(parent)
                if not p:
                    continue
                try:
                    label = self.tree.item(it, 'text')
                    cyc = int(label.split()[1])
                except Exception:
                    continue
                return p, cyc
        return None, None

    def _all_checked_cycles(self):
        out = []
        for it, ch in self.checked.items():
            if not ch or self.tree.item(it, 'values')[1] != 'cycle':
                continue
            parent = self.tree.parent(it)
            p = self.tree.item(parent, 'values')[3] or self.file_nodes.get(parent)
            if not p:
                continue
            try:
                label = self.tree.item(it, 'text')
                cyc = int(label.split()[1])
            except Exception:
                continue
            out.append((p, cyc))
        return out

    def _find_peaks_in_window(self, E, I, e_low, e_high, kind='anodic'):
        E = np.asarray(E, dtype=float)
        I = np.asarray(I, dtype=float)
        m = (E >= e_low) & (E <= e_high) & np.isfinite(E) & np.isfinite(I)
        if np.sum(m) < 3:
            return np.nan, np.nan
        Ew = E[m]; Iw = I[m]
        if kind == 'anodic':
            idx = int(np.argmax(Iw))
        else:
            idx = int(np.argmin(Iw))
        return float(Ew[idx]), float(Iw[idx])
    
    
    def _onset_k_sigma_integrator(self, E, I, e_low, e_high, ks=(2,3,4,5,8,12,15,20), forward_only=True, min_run=None):
        """
        Returns dict: {'onset_E_k2': ..., 'onset_E_k3': ..., ...} using an integrator-on-d2 method.
        - Restricts to forward scan (up to max(E)) if forward_only=True.
        - Restricts to E window [e_low, e_high].
        - Uses smoothing for derivatives if enabled (same as derivatives tab).
        """
        E = np.asarray(E, dtype=float)
        I = np.asarray(I, dtype=float)
    
        out = {f'onset_E_k{int(k)}': np.nan for k in ks}
        if len(E) < 10:
            return out
    
        # forward-only branch
        if forward_only:
            idx_sw = int(np.argmax(E))
            E = E[:idx_sw+1]
            I = I[:idx_sw+1]
    
        # window + finite
        m = (E >= e_low) & (E <= e_high) & np.isfinite(E) & np.isfinite(I)
        if np.sum(m) < 10:
            return out
    
        Ew = E[m]
        Iw = I[m]
    
        # sort by E for stable derivative/integration
        o = np.argsort(Ew)
        Ew = Ew[o]
        Iw = Iw[o]
    
        # min_run from UI (default to your onset_min_run)
        if min_run is None:
            try:
                min_run = int(self.onset_min_run.get())
            except Exception:
                min_run = 4
        min_run = max(1, int(min_run))
    
        # Smooth current for derivative computation (uses your smoothing toggle)
        Is = self._maybe_smooth(Iw, for_what='derivatives')
    
        # 2nd derivative
        d1 = safe_gradient(Is, Ew)
        d2 = safe_gradient(d1, Ew)
    
        # integrate only positive curvature
        dE = np.gradient(Ew)
        dE[~np.isfinite(dE)] = 0.0
        integ = np.cumsum(np.maximum(d2, 0.0) * np.abs(dE))
    
        # baseline = first 20% of points (at least 10)
        n0 = max(10, int(0.2 * len(integ)))
        base = integ[:n0]
        mu = float(np.nanmean(base))
        sd = float(np.nanstd(base))
        if (not np.isfinite(mu)) or (not np.isfinite(sd)) or sd <= 1e-12:
            return out
    
        above_cache = {}  # reuse comparisons per k
        for k in ks:
            thr = mu + float(k) * sd
            above = integ >= thr
            above_cache[k] = above
    
            run = 0
            for i, a in enumerate(above):
                if a:
                    run += 1
                    if run >= min_run:
                        out[f'onset_E_k{int(k)}'] = float(Ew[i - min_run + 1])
                        break
                else:
                    run = 0
    
        return out


    def _onset_percentiles(self, E, I, e_low, e_high, percentiles=(5,10,15,25), forward_only=True):
        """
        Onset = first E within [e_low, e_high] where I crosses a threshold.
        Threshold is defined as: I_min + (p/100)*(I_max - I_min) computed within the selected region.
    
        forward_only=True restricts to forward scan (up to max(E)), preventing mixing with backward branch.
        """
        E = np.asarray(E, dtype=float)
        I = np.asarray(I, dtype=float)
        if len(E) < 5:
            return {f'onset_E_p{p}': np.nan for p in percentiles}
    
        # Optional: restrict to forward scan
        if forward_only:
            idx_sw = int(np.argmax(E))  # switching potential index
            E = E[:idx_sw+1]
            I = I[:idx_sw+1]
    
        m = (E >= e_low) & (E <= e_high) & np.isfinite(E) & np.isfinite(I)
        if np.sum(m) < 5:
            return {f'onset_E_p{p}': np.nan for p in percentiles}
    
        Ew = E[m]
        Iw = I[m]
    
        # Sort by E (monotonic in forward branch, so this is safe)
        order = np.argsort(Ew)
        Ew = Ew[order]
        Iw = Iw[order]
    
        Imin = float(np.nanmin(Iw))
        Imax = float(np.nanmax(Iw))
        if not np.isfinite(Imin) or not np.isfinite(Imax) or (Imax <= Imin):
            return {f'onset_E_p{p}': np.nan for p in percentiles}
    
        res = {}
        for p in percentiles:
            thr = Imin + (float(p)/100.0) * (Imax - Imin)
            idxs = np.where(Iw >= thr)[0]
            res[f'onset_E_p{p}'] = float(Ew[idxs[0]]) if len(idxs) else np.nan
    
        return res


    def run_analysis(self, pairs=None, interactive=True):
        """Analyze `pairs` [(path, cycle), ...].

        pairs=None       -> use whatever is checked in the tree (button use).
        interactive=False -> no prompts; used by the batch runner, which does
        its own confirming one file at a time.
        Returns the number of rows produced.
        """
        if pairs is None:
            pairs = self._all_checked_cycles()
        if not pairs:
            if interactive:
                messagebox.showwarning("Analysis", "Check at least one cycle in the tree.")
            return 0

        # Warn up front if these files were already exported, so a long run
        # is not spent on work that is about to collide on save.
        if interactive and not self._warn_if_already_exported(pairs):
            return 0
    
        try:
            pH = float(self.ph_var.get())
        except Exception:
            pH = 14.0
    
        try:
            e_cath_low = float(self.cath_low.get()); e_cath_high = float(self.cath_high.get())
            e_an_low = float(self.an_low.get()); e_an_high = float(self.an_high.get())
        except Exception:
            e_cath_low, e_cath_high, e_an_low, e_an_high = 1.25, 1.40, 1.35, 1.50
    
        # clear table
        for it in self.analysis_tree.get_children():
            self.analysis_tree.delete(it)
    
        self.analysis_rows.clear()
        self.analysis_map.clear()
        self.analysis_text.delete('1.0', 'end')
    
        def fmt(v):
            if isinstance(v, (int, np.integer)):
                return str(int(v))
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                return ''
            if isinstance(v, float):
                return f"{v:.5g}"
            return str(v)
    
        cols = self.analysis_tree['columns']
    
        for (p, cyc) in pairs:
            E_rhe, t_sel, Y_raw, Y_corr, mask_corr, I_for_peaks, I_sel, dbg, is_ocv = self.get_display_series(p, cyc, pH)
            if E_rhe is None or len(E_rhe) == 0:
                continue
    
            E_rhe = np.asarray(E_rhe, dtype=float)
            t_sel = np.asarray(t_sel, dtype=float)
            I_raw = np.asarray(I_sel, dtype=float)
    
            # enforce same length for E/I/t before any masking
            m0 = min(len(E_rhe), len(I_raw), len(t_sel))
            if m0 < 5:
                continue
            E0 = E_rhe[:m0]
            I0 = I_raw[:m0]
            t0 = t_sel[:m0]
    
            # ============================================================
            # CORR metrics: ALWAYS attempt OER correction (independent of checkbox)
            # -> become NaN if correction isn't valid (or OCV)
            # ============================================================
            corr_valid = False
            I_corr = np.full_like(I0, np.nan, dtype=float)
    
            if not is_ocv:
                I_corr_try, mask_try, dbg_try = oer_subtract_modeA(E0, I0)
                corr_valid = bool(dbg_try.get('valid', False)) and int(dbg_try.get('applied', 0)) >= 3
                if corr_valid:
                    I_corr = np.asarray(I_corr_try, dtype=float)
    
            # ---------------------------
            # Clean RAW and keep the mask
            # ---------------------------
            E_use, I_raw2, m_raw, _ = clean_xy(E0, I0)  # m_raw is boolean mask for the original E0/I0/t0
            if E_use is None or len(E_use) < 5:
                continue
            t_raw = t0[:m0][m_raw]
            
            # ============================================================
            # FIXED-E sampling (RAW): I at specific potentials for forward/backward/for-back
            # + junction around 1.30-1.37 V
            # ============================================================
            
            # Use time-ordered cleaned mask points (still time-order, because we take from E0/I0 with m_raw)
            E_seq_raw = E0[m_raw]
            I_seq_raw = I0[m_raw]
            
            fixed_raw = compute_fixedE_metrics(E_seq_raw, I_seq_raw, targets=TARGET_E, junction_range=JUNCTION_RANGE)

    
            # -----------------------------------------
            # Force CORR to use the SAME E-points as RAW
            # -----------------------------------------
            if corr_valid:
                E_corr_masked = E0[m_raw]
                I_corr_masked = I_corr[m_raw]
    
                E_corr_use, I_corr2, m_corr2, _ = clean_xy(E_corr_masked, I_corr_masked)
                if E_corr_use is None or len(E_corr_use) < 5:
                    corr_valid = False
                    E_corr_use = np.array([], dtype=float)
                    I_corr2 = np.array([], dtype=float)
                    t_corr = np.array([], dtype=float)
                else:
                    t_corr = t_raw[m_corr2]
            else:
                E_corr_use = np.array([], dtype=float)
                I_corr2 = np.array([], dtype=float)
                t_corr = np.array([], dtype=float)
    
            # ---------------------------
            # RAW peaks + RAW onset + RAW Q
            # ---------------------------
            anod_E_raw, anod_I_raw = self._find_peaks_in_window(E_use, I_raw2, e_an_low, e_an_high, kind='anodic')
            cath_E_raw, cath_I_raw = self._find_peaks_in_window(E_use, I_raw2, e_cath_low, e_cath_high, kind='cathodic')
            onset_raw = self._onset_percentiles(E_use, I_raw2, e_an_low, e_an_high, percentiles=(5,10,15,25))
            onset_raw_k = self._onset_k_sigma_integrator(E_use, I_raw2, e_an_low, e_an_high,
                                                         ks=(2,3,4,5,8,12,15,20),
                                                         forward_only=True,
                                                         min_run=int(self.onset_min_run.get()) if hasattr(self, "onset_min_run") else 4)
            Q_raw_curve = integrate_Q(I_raw2, t_raw)
            Q_end_raw = float(Q_raw_curve[-1]) if len(Q_raw_curve) else np.nan
    
            # ---------------------------
            # CORR peaks + CORR onset + CORR Q (or NaNs)
            # ---------------------------
            if corr_valid:
                anod_E_corr, anod_I_corr = self._find_peaks_in_window(E_corr_use, I_corr2, e_an_low, e_an_high, kind='anodic')
                cath_E_corr, cath_I_corr = self._find_peaks_in_window(E_corr_use, I_corr2, e_cath_low, e_cath_high, kind='cathodic')
                onset_corr = self._onset_percentiles(E_corr_use, I_corr2, e_an_low, e_an_high, percentiles=(5,10,15,25))
                onset_corr_k = self._onset_k_sigma_integrator(E_corr_use, I_corr2, e_an_low, e_an_high,
                                                              ks=(2,3,4,5,8,12,15,20),
                                                              forward_only=True,
                                                              min_run=int(self.onset_min_run.get()) if hasattr(self, "onset_min_run") else 4)
                Q_corr_curve = integrate_Q(I_corr2, t_corr)
                Q_end_corr = float(Q_corr_curve[-1]) if len(Q_corr_curve) else np.nan
            else:
                anod_E_corr = anod_I_corr = np.nan
                cath_E_corr = cath_I_corr = np.nan
                onset_corr = {f'onset_E_p{pp}': np.nan for pp in (5,10,15,25)}
                onset_corr_k = {f'onset_E_k{int(k)}': np.nan for k in (2,3,4,5,8,12,15,20)}
                Q_end_corr = np.nan
    
            # ---------------------------
            # ratios and deltaEp
            # ---------------------------
            I_ratio_raw = (anod_I_raw / abs(cath_I_raw)) if (np.isfinite(anod_I_raw) and np.isfinite(cath_I_raw) and cath_I_raw != 0) else np.nan
            I_ratio_corr = (anod_I_corr / abs(cath_I_corr)) if (np.isfinite(anod_I_corr) and np.isfinite(cath_I_corr) and cath_I_corr != 0) else np.nan
    
            dEp_raw = (anod_E_raw - cath_E_raw) if (np.isfinite(anod_E_raw) and np.isfinite(cath_E_raw)) else np.nan
            dEp_corr = (anod_E_corr - cath_E_corr) if (np.isfinite(anod_E_corr) and np.isfinite(cath_E_corr)) else np.nan
    
            # ---------------------------
            # Build row dict matching the columns list in _build_analysis_tab
            # ---------------------------
                        # ---------------------------
            # Build row dict matching the columns list in _build_analysis_tab
            # ---------------------------
            row = {
                'file': os.path.basename(p),
                'cycle': int(cyc),

                'anodic_E_raw': anod_E_raw, 'anodic_I_raw': anod_I_raw,
                'cathodic_E_raw': cath_E_raw, 'cathodic_I_raw': cath_I_raw,

                'anodic_E_corr': anod_E_corr, 'anodic_I_corr': anod_I_corr,
                'cathodic_E_corr': cath_E_corr, 'cathodic_I_corr': cath_I_corr,

                'I_ratio_raw': I_ratio_raw,
                'I_ratio_corr': I_ratio_corr,

                'dEp_raw': dEp_raw,
                'dEp_corr': dEp_corr,

                'Q_end_raw': Q_end_raw,
                'Q_end_corr': Q_end_corr,

                # Percentile onsets
                'onset_raw_E_p5':  onset_raw.get('onset_E_p5',  np.nan),
                'onset_raw_E_p10': onset_raw.get('onset_E_p10', np.nan),
                'onset_raw_E_p15': onset_raw.get('onset_E_p15', np.nan),
                'onset_raw_E_p25': onset_raw.get('onset_E_p25', np.nan),

                'onset_corr_E_p5':  onset_corr.get('onset_E_p5',  np.nan),
                'onset_corr_E_p10': onset_corr.get('onset_E_p10', np.nan),
                'onset_corr_E_p15': onset_corr.get('onset_E_p15', np.nan),
                'onset_corr_E_p25': onset_corr.get('onset_E_p25', np.nan),

                # k-sigma onsets (RAW)
                'onset_raw_E_k2':  onset_raw_k.get('onset_E_k2',  np.nan),
                'onset_raw_E_k3':  onset_raw_k.get('onset_E_k3',  np.nan),
                'onset_raw_E_k4':  onset_raw_k.get('onset_E_k4',  np.nan),
                'onset_raw_E_k5':  onset_raw_k.get('onset_E_k5',  np.nan),
                'onset_raw_E_k8':  onset_raw_k.get('onset_E_k8',  np.nan),
                'onset_raw_E_k12': onset_raw_k.get('onset_E_k12', np.nan),
                'onset_raw_E_k15': onset_raw_k.get('onset_E_k15', np.nan),
                'onset_raw_E_k20': onset_raw_k.get('onset_E_k20', np.nan),

                # k-sigma onsets (CORR)
                'onset_corr_E_k2':  onset_corr_k.get('onset_E_k2',  np.nan),
                'onset_corr_E_k3':  onset_corr_k.get('onset_E_k3',  np.nan),
                'onset_corr_E_k4':  onset_corr_k.get('onset_E_k4',  np.nan),
                'onset_corr_E_k5':  onset_corr_k.get('onset_E_k5',  np.nan),
                'onset_corr_E_k8':  onset_corr_k.get('onset_E_k8',  np.nan),
                'onset_corr_E_k12': onset_corr_k.get('onset_E_k12', np.nan),
                'onset_corr_E_k15': onset_corr_k.get('onset_E_k15', np.nan),
                'onset_corr_E_k20': onset_corr_k.get('onset_E_k20', np.nan),
                # Add fixed-E metrics into the row (RAW)
                # Prefix them to avoid collisions with other metrics
                

            }
            for k, v in fixed_raw.items():
                row["raw_" + k] = v

            # Keep row for later use
            self.analysis_rows.append(row)
            self.analysis_map[(p, cyc)] = row

            # Insert row into table (respect the visible columns order)
            vals = [fmt(row.get(c, np.nan)) for c in cols]
            self.analysis_tree.insert('', 'end', values=vals)

        # Statistics report
        self._build_analysis_report()

        # Second heads-up now the run is done and Export is the next step.
        if interactive:
            self._notify_after_analysis()

        return len(self.analysis_rows)


    
        
    

    # =================== Export / rename helpers ===================
    @staticmethod
    def _find_case_insensitive(folder, filename):
        """Return the real name of `filename` inside `folder`, matched without
        regard to case, or None. Windows file names are case-insensitive, so
        "Nitrite_A..." and "nitrite_a..." are the SAME file - comparing with
        == would miss the clash and silently overwrite."""
        try:
            if not os.path.isdir(folder):
                return None
            want = filename.lower()
            for existing in os.listdir(folder):
                if existing.lower() == want:
                    return existing
        except Exception:
            pass
        return None

    def _existing_export(self, base):
        """Real name of an already-exported CSV for this base name, or None."""
        return self._find_case_insensitive(EXPORT_FOLDER, base + '.csv')

    def _guess_base_for_path(self, path):
        """The export base name that `path` would produce, or '' if the file
        name does not yield the required parts."""
        g = guess_name_parts(os.path.basename(path))
        parts = [g['anion'], g['station'], g['electrode'], g['meas'], g['laborant']]
        if not all(str(x).strip() for x in parts):
            return ''
        return "_".join(ExportNameDialog._clean(x) for x in parts)

    def _already_exported_hits(self, paths):
        """[(source .mpt name, existing CSV name), ...] for files whose
        guessed export name already exists. Files whose name cannot be
        parsed are skipped rather than falsely flagged."""
        seen = []
        hits = []
        for p_ in paths:
            if p_ in seen:
                continue
            seen.append(p_)
            base = self._guess_base_for_path(p_)
            if not base:
                continue
            existing = self._existing_export(base)
            if existing:
                hits.append((os.path.basename(p_), existing))
        return hits

    @staticmethod
    def _fmt_export_hits(hits):
        return "\n".join(f"   {src}\n      -> already exported as  {dst}" for src, dst in hits)

    def _warn_if_already_exported(self, pairs):
        """CHECK 1 - before a run. If the guessed name for any source file
        already has a CSV in the export folder, say so and let the user
        decide. Returns True to continue, False to abort."""
        hits = self._already_exported_hits([p_ for (p_, _c) in pairs])
        if not hits:
            return True

        return messagebox.askyesno(
            "Already exported",
            f"{len(hits)} of the selected file(s) appear to be exported already:\n\n"
            f"{self._fmt_export_hits(hits)}\n\nRun the analysis anyway?")

    def _notify_after_analysis(self):
        """CHECK 2 - right after a run finishes, so the collision is fresh in
        mind before the user reaches Export."""
        hits = self._already_exported_hits(list(dict.fromkeys(
            p_ for (p_, _c) in self.analysis_map.keys())))
        if not hits:
            return
        messagebox.showwarning(
            "Already exported",
            f"Analysis finished.\n\n{len(hits)} of the analyzed file(s) already have "
            f"an export:\n\n{self._fmt_export_hits(hits)}\n\n"
            "You will be asked again before anything is overwritten.")

    def _rename_mpt_file(self, old_path, new_base):
        """Rename a .mpt on disk. Disk only - no in-memory state is patched
        here, because the caller reloads the library afterwards, which
        rebuilds every path, node and cache entry from scratch.
        Returns (ok, message)."""
        folder = os.path.dirname(old_path)
        # Keep the file's own extension (.mpt, or .csv / .txt from MPR Converter).
        new_name = new_base + (os.path.splitext(old_path)[1] or '.mpt')
        new_path = os.path.join(folder, new_name)

        if os.path.normcase(new_path) == os.path.normcase(old_path):
            return True, "already named correctly"

        # Case-insensitive: on Windows two names differing only in case are
        # the same file, so os.rename would silently clobber the other one.
        clash = self._find_case_insensitive(folder, new_name)
        if clash:
            return False, f"'{clash}' already exists in that folder"

        try:
            os.rename(old_path, new_path)
        except Exception as e:
            return False, str(e)

        return True, new_name

    def _reload_library(self):
        """Clear the tree and re-open every root folder. Used after a batch
        rename, where the paths held in memory no longer exist on disk."""
        roots = list(self.open_roots)
        if not roots:
            return False
        self.clear_all()
        for r in roots:
            if os.path.isdir(r):
                self.open_folder(path_override=r)
        return True

    def _collect_tree_mpt_files(self, top_level_only=False):
        """Every .mpt currently listed in the tree, in display order.

        top_level_only=True keeps only files sitting directly inside one of
        the opened root folders, ignoring anything nested in a subfolder.
        The filter is done on the path rather than on tree depth, because
        open_folder inserts every walked directory as a top-level node, so
        tree position does not reflect the folder hierarchy.
        """
        out = []

        def walk(parent):
            for ch in self.tree.get_children(parent):
                p = self.file_nodes.get(ch)
                if p:
                    if p not in [x[1] for x in out]:
                        out.append((ch, p))
                else:
                    walk(ch)

        walk('')

        if not top_level_only:
            return out

        roots = {os.path.normcase(os.path.normpath(r)) for r in self.open_roots}
        if not roots:
            return out

        kept = []
        for node, path in out:
            parent_dir = os.path.normcase(os.path.normpath(os.path.dirname(path)))
            if parent_dir in roots:
                kept.append((node, path))
        return kept

    def _overlay_candidates(self):
        """[(node, path, meta)] for every .mpt listed in the tree, with the
        name parts read from the file name. Nothing is parsed here - the
        overlay window only loads the files the user actually selects."""
        out = []
        # top_level_only: subfolders are not searched, matching the rename
        # and batch-analysis modes.
        for node, path in self._collect_tree_mpt_files(top_level_only=True):
            out.append((node, path, guess_name_parts(os.path.basename(path))))
        return out

    def open_overlay_window(self):
        if not self._collect_tree_mpt_files(top_level_only=True):
            messagebox.showwarning("Overlay", "No .mpt files listed.\nOpen a folder first.")
            return
        # Reuse the existing window instead of stacking duplicates.
        if getattr(self, '_overlay_win', None) is not None:
            try:
                if self._overlay_win.winfo_exists():
                    self._overlay_win.deiconify()
                    self._overlay_win.lift()
                    return
            except Exception:
                pass
        self._overlay_win = OverlayWindow(self)

    def _cycles_of(self, path, node):
        """Every cycle number in `path`, loading the file if needed."""
        if not self._ensure_file_loaded(path, node):
            return []
        try:
            return sorted(set(int(c) for c in self.cache[path]['cyc']))
        except Exception:
            return []

    def find_unanalyzed_mpt(self):
        """Find .mpt files in the opened folder with no matching CSV in the
        export folder, then offer to analyze them one at a time."""
        files = self._collect_tree_mpt_files(top_level_only=True)
        if not files:
            messagebox.showwarning("Unanalyzed files",
                                   "No .mpt files directly in the opened folder.")
            return

        queue = []       # (node, path) still to analyze
        unnamed = []     # file names we could not read the parts from
        done = 0
        for node, path in files:
            base = self._guess_base_for_path(path)
            if not base:
                # Name unreadable -> cannot tell if it was exported. Treat as
                # unanalyzed so it is not silently ignored.
                unnamed.append(os.path.basename(path))
                queue.append((node, path))
            elif self._existing_export(base):
                done += 1
            else:
                queue.append((node, path))

        if not queue:
            messagebox.showinfo(
                "Unanalyzed files",
                f"Nothing to do.\nAll {done} file(s) in the folder already have an export.")
            return

        listing = "\n".join("   " + os.path.basename(p) for _n, p in queue[:20])
        if len(queue) > 20:
            listing += f"\n   ... and {len(queue) - 20} more"

        extra = ""
        if unnamed:
            extra = (f"\n\n{len(unnamed)} of these have a file name the parts could not be "
                     f"read from - fill those in by hand.")

        msg = (f"{len(queue)} file(s) have no export yet"
               f"{f' ({done} already done)' if done else ''}:\n\n{listing}{extra}\n\n"
               "Analyze them one at a time?\nYou will review each report before it is saved.")
        if not messagebox.askyesno("Unanalyzed files", msg):
            return

        self._run_batch_analysis(queue)

    def _run_batch_analysis(self, queue):
        """Analyze each queued file, show its report for approval, and save
        only what the user approves."""
        saved = skipped = failed = 0
        problems = []

        for i, (node, path) in enumerate(queue, start=1):
            fname = os.path.basename(path)
            progress = f"{i} of {len(queue)}"

            cycles = self._cycles_of(path, node)
            if not cycles:
                failed += 1
                problems.append(f"{fname}  ->  no readable cycles")
                continue

            n = self.run_analysis(pairs=[(path, c) for c in cycles], interactive=False)
            if not n:
                failed += 1
                problems.append(f"{fname}  ->  analysis produced no rows")
                continue

            report = self.analysis_text.get('1.0', 'end-1c')

            dlg = ExportNameDialog(self.root, [fname], mode='batch',
                                   progress=progress, report=report)

            if dlg.action == 'cancel':
                break
            if dlg.action == 'skip' or not dlg.result:
                skipped += 1
                continue

            base = dlg.result
            fn = dlg.full_path()
            try:
                os.makedirs(os.path.dirname(fn), exist_ok=True)
            except Exception as e:
                failed += 1
                problems.append(f"{fname}  ->  {e}")
                continue

            clash = self._existing_export(base)
            if clash:
                if not messagebox.askyesno(
                        "File exists",
                        f"{clash}\n\nalready exists in the export folder.\n\nOverwrite it?"):
                    skipped += 1
                    continue
                fn = os.path.join(EXPORT_FOLDER, clash)

            ok, err = self._write_analysis_csv(fn)
            if ok:
                saved += 1
            else:
                failed += 1
                problems.append(f"{fname}  ->  {err}")

        summary = f"Saved: {saved}\nSkipped: {skipped}\nFailed: {failed}"
        if problems:
            summary += "\n\n" + "\n".join(problems[:15])
            if len(problems) > 15:
                summary += f"\n... and {len(problems) - 15} more"
        messagebox.showinfo("Batch analysis", summary)

    def batch_rename_mpt(self):
        """Walk the listed .mpt files one by one, propose a name from the
        current file name, and rename only what the user approves.
        No analysis is run in this mode."""
        # Only files sitting directly in the opened folder - subfolders are
        # left alone.
        all_files = self._collect_tree_mpt_files()
        files = self._collect_tree_mpt_files(top_level_only=True)
        in_subfolders = len(all_files) - len(files)

        if not files:
            if in_subfolders:
                messagebox.showwarning(
                    "Rename MPT files",
                    f"No .mpt files directly in the opened folder.\n\n"
                    f"{in_subfolders} file(s) are in subfolders and are not renamed.\n"
                    "Open the subfolder itself to rename those.")
            else:
                messagebox.showwarning("Rename MPT files",
                                       "No .mpt files listed.\nOpen a folder first.")
            return

        renamed = skipped = failed = 0
        problems = []

        for i, (node, path) in enumerate(files, start=1):
            if not os.path.exists(path):
                failed += 1
                problems.append(f"{os.path.basename(path)}  ->  file no longer on disk")
                continue

            dlg = ExportNameDialog(
                self.root, [os.path.basename(path)],
                mode='rename',
                folder=os.path.dirname(path),
                progress=f"{i} of {len(files)}",
            )

            if dlg.action == 'cancel':
                break
            if dlg.action == 'skip' or not dlg.result:
                skipped += 1
                continue

            ok, msg = self._rename_mpt_file(path, dlg.result)
            if ok:
                renamed += 1
            else:
                failed += 1
                problems.append(f"{os.path.basename(path)}  ->  {msg}")

        summary = f"Renamed: {renamed}\nSkipped: {skipped}\nFailed: {failed}"
        if in_subfolders:
            summary += f"\n\n{in_subfolders} file(s) in subfolders were not touched."
        if problems:
            summary += "\n\n" + "\n".join(problems[:15])
            if len(problems) > 15:
                summary += f"\n... and {len(problems) - 15} more"

        # Reload so the tree reflects what is actually on disk. Every cached
        # path is stale after a rename, so a reload is safer than patching.
        if renamed:
            if self._reload_library():
                summary += "\n\nLibrary reloaded."
            else:
                summary += "\n\nReload the folder to see the new names."

        messagebox.showinfo("Rename MPT files", summary)

    # =================== Analysis report ===================
    def show_param_help(self):
        win = getattr(self, '_param_help', None)
        try:
            if win is not None and win.winfo_exists():
                win.deiconify()
                win.lift()
                return
        except Exception:
            pass
        lang = overlay_load_settings().get('tutorial_lang', 'en')
        self._param_help = ParamHelpWindow(self, lang)

    def _copy_report(self):
        """Put the whole report on the clipboard (easy paste into the lab notebook)."""
        try:
            txt = self.analysis_text.get('1.0', 'end-1c')
            self.root.clipboard_clear()
            self.root.clipboard_append(txt)
            messagebox.showinfo("Report", "Report copied to clipboard.")
        except Exception as e:
            messagebox.showwarning("Report", f"Could not copy: {e}")

    def _build_analysis_report(self):
        """Report shown in the Analysis tab after a run.

        Order is deliberate: the few parameters that matter most come first
        in bold, then everything that went wrong (missing values, failed OER
        correction), and only then the full per-parameter statistics table.
        Problems should never be buried below ~90 rows of numbers.

        Per parameter: n (finite), missing count and %, mean, std (sample,
        ddof=1), variance, min, max, median, and CV% (std/|mean|, scale-free
        so currents at different potentials stay comparable).
        """
        T = self.analysis_text
        T.delete('1.0', 'end')

        rows = self.analysis_rows
        if not rows:
            T.insert('end', "No analysis results.\n")
            return

        # Text tags for emphasis
        T.tag_configure('bold', font=('Consolas', 9, 'bold'))
        T.tag_configure('big', font=('Consolas', 10, 'bold'))
        T.tag_configure('bad', foreground='#b00020')
        T.tag_configure('good', foreground='#0a7d33')

        cols = list(self.analysis_tree['columns'])
        n_rows = len(rows)

        # Header defined up front so every rule matches the table width.
        hdr = (f"{'parameter':<26}{'n':>5}{'miss':>6}{'miss%':>7}"
               f"{'mean':>14}{'std':>13}{'variance':>14}"
               f"{'min':>14}{'max':>14}{'median':>14}{'CV%':>9}")
        W = len(hdr)

        def stats_of(c):
            """(n, miss, mean, std, var, min, max, median, cv) or None."""
            vals = np.array([self._as_float(r.get(c, np.nan)) for r in rows], dtype=float)
            fin = vals[np.isfinite(vals)]
            if fin.size == 0:
                return None
            n = fin.size
            mean = float(np.mean(fin))
            std = float(np.std(fin, ddof=1)) if n > 1 else np.nan
            var = float(np.var(fin, ddof=1)) if n > 1 else np.nan
            cv = (abs(std / mean) * 100.0) if (n > 1 and np.isfinite(std) and mean != 0) else np.nan
            return (n, n_rows - n, mean, std, var,
                    float(np.min(fin)), float(np.max(fin)), float(np.median(fin)), cv)

        def stat_line(c, st):
            n, miss, mean, std, var, mn, mx, md, cv = st
            return (f"{c:<26}{n:>5}{miss:>6}{100.0 * miss / n_rows:>6.0f}%"
                    f"{self._g(mean):>14}{self._g(std):>13}{self._g(var):>14}"
                    f"{self._g(mn):>14}{self._g(mx):>14}{self._g(md):>14}{self._g(cv):>9}\n")

        files = []
        for r in rows:
            f = str(r.get('file', ''))
            if f and f not in files:
                files.append(f)

        # ======================= 1. header =======================
        T.insert('end', "=" * W + "\n")
        T.insert('end', f"ANALYSIS REPORT   -   {n_rows} cycle(s)   |   {len(files)} file(s)\n", 'bold')
        T.insert('end', "=" * W + "\n\n")

        for f in files:
            cy = sorted(int(r['cycle']) for r in rows if str(r.get('file', '')) == f
                        and r.get('cycle') is not None)
            if cy:
                T.insert('end', f"  {f}   ->  {len(cy)} cycle(s): {self._fmt_ranges(cy)}\n")
            else:
                T.insert('end', f"  {f}\n")
        T.insert('end', "\n")

        # ================== 2. key parameters ====================
        KEY_COLS = [
            ('cathodic_E_raw', 'cathodic peak potential (raw)'),
            ('raw_I_fwd_1p4', 'I forward @ 1.40 V'),
            ('raw_I_bwd_1p4', 'I backward @ 1.40 V'),
            ('raw_I_for_back_1p4', 'I fwd-bwd @ 1.40 V'),
            ('raw_I_fwd_1p45', 'I forward @ 1.45 V'),
            ('raw_I_bwd_1p45', 'I backward @ 1.45 V'),
            ('raw_I_for_back_1p45', 'I fwd-bwd @ 1.45 V'),
        ]
        T.insert('end', "KEY PARAMETERS\n", 'big')
        T.insert('end', "-" * W + "\n")
        T.insert('end', hdr + "\n", 'bold')
        T.insert('end', "-" * W + "\n")
        for c, label in KEY_COLS:
            if c not in cols:
                continue
            st = stats_of(c)
            if st is None:
                T.insert('end', f"{c:<26}{'- no finite value in any cycle -':>40}\n", 'bad')
            else:
                T.insert('end', stat_line(c, st), 'bold')
        T.insert('end', "-" * W + "\n")
        T.insert('end', "    cathodic_E_raw = cathodic peak potential (raw curve)\n")
        T.insert('end', "    raw_I_fwd_1p4 / _1p45 = current on the forward sweep at 1.40 / 1.45 V\n")
        T.insert('end', "    _bwd = backward sweep,  _for_back = forward minus backward\n\n")

        # ============ 3. errors and missing values ===============
        skip = {'file', 'cycle'}
        stats_rows = []
        all_missing = []
        for c in cols:
            if c in skip:
                continue
            st = stats_of(c)
            if st is None:
                all_missing.append(c)
            else:
                stats_rows.append((c, st))

        T.insert('end', "ERRORS & MISSING VALUES\n", 'big')
        T.insert('end', "-" * W + "\n")

        # Every *_corr metric is blank when the Mode-A correction did not
        # validate, so this single number explains most missing values.
        corr_ok = sum(1 for r in rows
                      if np.isfinite(self._as_float(r.get('anodic_E_corr', np.nan))))
        tag = 'good' if corr_ok == n_rows else 'bad'
        T.insert('end',
                 f"  OER (Mode-A) correction valid for {corr_ok}/{n_rows} cycle(s) "
                 f"({100.0 * corr_ok / n_rows:.0f}%)"
                 + ("" if corr_ok == n_rows else "   ->  all *_corr metrics blank for the rest")
                 + "\n\n", tag)

        partial = [(c, st[1]) for c, st in stats_rows if st[1] > 0]

        if all_missing:
            T.insert('end', f"  MISSING ENTIRELY - no finite value in any cycle "
                            f"({len(all_missing)} parameter(s)):\n", 'bad')
            for i in range(0, len(all_missing), 3):
                T.insert('end', "      " + "".join(f"{x:<38}" for x in all_missing[i:i + 3]).rstrip() + "\n")
            T.insert('end', "\n")

        if partial:
            T.insert('end', f"  PARTIALLY MISSING ({len(partial)} parameter(s)):\n", 'bad')
            for c, miss in sorted(partial, key=lambda x: -x[1]):
                T.insert('end', f"      {c:<40} {miss}/{n_rows} blank "
                                f"({100.0 * miss / n_rows:.0f}%)\n")
            T.insert('end', "\n")

        if not all_missing and not partial:
            T.insert('end', "  No missing values - every parameter is finite in every cycle.\n\n", 'good')

        # ================ 4. full statistics =====================
        T.insert('end', "ALL PARAMETERS\n", 'big')
        T.insert('end', "-" * W + "\n")
        T.insert('end', hdr + "\n", 'bold')
        T.insert('end', "-" * W + "\n")
        for c, st in stats_rows:
            T.insert('end', stat_line(c, st))

        T.insert('end', "\n" + "-" * W + "\n")
        T.insert('end', "Note: *_corr metrics come from the Mode-A OER correction attempt and are\n"
                        "independent of the 'Subtract OER' checkbox; an invalid correction leaves them blank.\n")
        T.see('1.0')

    @staticmethod
    def _as_float(v):
        try:
            return float(v)
        except Exception:
            return np.nan

    @staticmethod
    def _g(v):
        """Compact fixed-width number formatting for the report table."""
        if v is None or not np.isfinite(v):
            return "-"
        a = abs(v)
        if a != 0 and (a < 1e-3 or a >= 1e5):
            return f"{v:.3e}"
        return f"{v:.5g}"

    @staticmethod
    def _fmt_ranges(nums):
        """[1,2,3,7,8] -> '1-3, 7-8' so long cycle lists stay readable."""
        if not nums:
            return ""
        out = []
        start = prev = nums[0]
        for x in nums[1:]:
            if x == prev + 1:
                prev = x
                continue
            out.append(f"{start}" if start == prev else f"{start}-{prev}")
            start = prev = x
        out.append(f"{start}" if start == prev else f"{start}-{prev}")
        return ", ".join(out)

    def plot_derivatives_for_first_checked(self):
        p, cyc = self._first_checked_cycle()
        if p is None:
            messagebox.showwarning("Derivatives", "Check at least one cycle in the tree.")
            return
        self._plot_derivatives_pair(p, cyc)
        self._deriv_last_pair = (p, cyc)



    def export_analysis_csv(self):
        if not self.analysis_rows:
            messagebox.showwarning("Export", "Run Analysis first.")
            return

        # Names of the .mpt file(s) covered by this analysis, shown at the
        # top of the dialog so it is clear what is being exported.
        names = []
        for r in self.analysis_rows:
            nm = str(r.get('file', '')).strip()
            if nm and nm not in names:
                names.append(nm)

        # Collect the name parts: Anion_Station_Electrode_Measurement_Laborant
        dlg = ExportNameDialog(self.root, names)
        if not dlg.result:
            return  # cancelled

        base = dlg.result

        # Write straight into EXPORT_FOLDER - no file browser.
        fn = dlg.full_path()

        try:
            os.makedirs(os.path.dirname(fn), exist_ok=True)
        except Exception as e:
            messagebox.showerror("Export", f"Could not reach the export folder:\n{e}")
            return

        # Case-insensitive: on Windows "Nitrite_A..." and "nitrite_a..." are
        # the same file, so a plain os.path.exists would miss the clash.
        clash = self._existing_export(base)
        if clash:
            if not messagebox.askyesno(
                    "File exists",
                    f"{clash}\n\nalready exists in the export folder.\n\nOverwrite it?"):
                return
            # Write to the name that is actually on disk, so we overwrite it
            # instead of leaving two files differing only in case.
            fn = os.path.join(EXPORT_FOLDER, clash)

        ok, err = self._write_analysis_csv(fn)
        if not ok:
            messagebox.showerror("Export", f"Could not write the file:\n{err}")
            return

        messagebox.showinfo("Export", f"Saved:\n{fn}")

    def _write_analysis_csv(self, fn):
        """Write the current analysis_rows to `fn`. Returns (ok, error)."""
        cols = list(self.analysis_tree['columns'])
        try:
            with open(fn, 'w', encoding='utf-8') as f:
                f.write(','.join(cols) + '\n')
                for r in self.analysis_rows:
                    row = []
                    for c in cols:
                        v = r.get(c, '')
                        if v is None or (isinstance(v, float) and not np.isfinite(v)):
                            row.append('')
                        else:
                            row.append(str(v))
                    f.write(','.join(row) + '\n')
        except Exception as e:
            # e.g. the file is open in Excel, or the Drive folder is read-only
            return False, str(e)
        return True, ''

    # =================== Comparison tab ===================
    def _build_comparison_tab(self, tab):
        top = tk.Frame(tab); top.pack(fill='x', padx=8, pady=6)

        # Multi-select list
        left = tk.Frame(top); left.pack(side='left', fill='y')
        tk.Label(left, text="Parameters (multi-select):").pack(anchor='w')

        self.cmp_list = tk.Listbox(left, selectmode='extended', height=9, exportselection=False)
        self.cmp_list.pack(side='left', fill='y')

        sb = tk.Scrollbar(left, orient='vertical', command=self.cmp_list.yview)
        sb.pack(side='left', fill='y')
        self.cmp_list.config(yscrollcommand=sb.set)

        # Single param fallback (kept)
        right = tk.Frame(top); right.pack(side='left', fill='both', expand=True, padx=(12,0))
        tk.Label(right, text="Single param (optional):").pack(anchor='w')
        self.cmp_param_menu = ttk.Combobox(
            right,
            textvariable=self.cmp_param,
            width=22,
            values=[
                'anodic_E_raw','anodic_I_raw','cathodic_E_raw','cathodic_I_raw',
                'anodic_E_corr','anodic_I_corr','cathodic_E_corr','cathodic_I_corr',
                'I_ratio_raw','I_ratio_corr','dEp_raw','dEp_corr','Q_end_raw','Q_end_corr',
                'onset_raw_E_p5','onset_raw_E_p10','onset_raw_E_p15','onset_raw_E_p25',
                'onset_corr_E_p5','onset_corr_E_p10','onset_corr_E_p15','onset_corr_E_p25',
                'onset_raw_E_k2','onset_raw_E_k3','onset_raw_E_k4','onset_raw_E_k5','onset_raw_E_k8','onset_raw_E_k12','onset_raw_E_k15','onset_raw_E_k20',
                'onset_corr_E_k2','onset_corr_E_k3','onset_corr_E_k4','onset_corr_E_k5','onset_corr_E_k8','onset_corr_E_k12','onset_corr_E_k15','onset_corr_E_k20',

            ]
        )
        self.cmp_param_menu.pack(anchor='w', pady=(2,0))

        # populate multi list with same values
        for v in self.cmp_param_menu['values']:
            self.cmp_list.insert('end', v)

        btnrow = tk.Frame(tab); btnrow.pack(fill='x', padx=8, pady=(0,6))
        tk.Button(btnrow, text="Build from Analysis + Plot", command=self.plot_comparison).pack(side='left')
        tk.Button(btnrow, text="Run Analysis Now", command=self.run_analysis).pack(side='left', padx=(6,0))
        tk.Label(btnrow, text="Tip: select multiple params with Ctrl/Shift.").pack(side='left', padx=12)

        plot_frame = tk.Frame(tab)
        plot_frame.pack(fill='both', expand=True, padx=8, pady=(0,8))

        fig = Figure(figsize=(6, 4))
        self.cmp_ax = fig.add_subplot(111)
        self.cmp_canvas = FigureCanvasTkAgg(fig, plot_frame)
        self.cmp_canvas.get_tk_widget().pack(fill='both', expand=True)
        tb = NavigationToolbar2Tk(self.cmp_canvas, plot_frame)
        tb.update()

    def plot_comparison(self):
        """
        Comparison tab:
        - plots chosen analysis parameters vs cycle number
        - each (param, file) is its own series
        - robustly discovers available params from analysis_rows / analysis_map
        - supports listbox selection OR single variable selection
        """
    
        # ---------------------------
        # Guard: need analysis results
        # ---------------------------
        if not getattr(self, "analysis_rows", None):
            # Try to build analysis_map/rows if they exist but empty
            messagebox.showwarning("Comparison", "No analysis results found.\nRun Analysis first (Analysis tab).")
            return
    
        # ---------------------------
        # Group rows by file -> cycle
        # ---------------------------
        by_file = {}
        for r in self.analysis_rows:
            fn = r.get("file", "")
            cyc = r.get("cycle", None)
            if not fn or cyc is None:
                continue
            by_file.setdefault(fn, {})[int(cyc)] = r
    
        if not by_file:
            messagebox.showwarning("Comparison", "No usable (file, cycle) rows found in analysis results.")
            return
    
        # ---------------------------
        # Discover all available params (keys)
        # ---------------------------
        # Exclude identifiers
        exclude = {"file", "cycle"}
        all_params = set()
        for fn, cyc_map in by_file.items():
            for cy, row in cyc_map.items():
                all_params.update([k for k in row.keys() if k not in exclude])
    
        all_params = sorted(all_params)
    
        # ----------------------------------------------------------
        # If you have a listbox in UI, ensure it contains the params
        # ----------------------------------------------------------
        # We handle two common patterns:
        # 1) self.cmp_param_listbox (tk.Listbox)
        # 2) self.cmp_param_var (StringVar / single param)
        params = []
    
        if hasattr(self, "cmp_param_listbox"):
            lb = self.cmp_param_listbox
    
            # If listbox is empty, populate it with all params we found
            try:
                if lb.size() == 0 and len(all_params):
                    for p in all_params:
                        lb.insert("end", p)
            except Exception:
                pass
    
            # Read selected items
            try:
                sel = list(lb.curselection())
                if sel:
                    params = [lb.get(i) for i in sel]
                else:
                    # Nothing selected -> be helpful: plot a small default set if available
                    default_candidates = [
                        "anodic_I_raw", "anodic_E_raw", "cathodic_I_raw", "cathodic_E_raw",
                        "Q_end_raw", "dEp_raw",
                    ]
                    # Add also a couple of fixed-current params if they exist
                    fixed_candidates = [p for p in all_params if p.startswith("Iraw_for_E") or p.startswith("Iraw_back_E")]
                    default_candidates += fixed_candidates[:4]
    
                    params = [p for p in default_candidates if p in all_params]
                    if not params and all_params:
                        params = [all_params[0]]
            except Exception:
                # Fallback if listbox behaves unexpectedly
                params = [p for p in ("anodic_I_raw", "Q_end_raw") if p in all_params] or (all_params[:1] if all_params else [])
    
        elif hasattr(self, "cmp_param_var"):
            # Single selection (combobox / entry bound to StringVar)
            try:
                s = str(self.cmp_param_var.get()).strip()
            except Exception:
                s = ""
            if not s:
                # default
                s = "anodic_I_raw" if "anodic_I_raw" in all_params else (all_params[0] if all_params else "")
            # allow comma-separated
            params = [x.strip() for x in s.split(",") if x.strip()]
            params = [p for p in params if p in all_params]
    
        else:
            # No UI binding found -> default to a small set
            params = [p for p in ("anodic_I_raw", "Q_end_raw") if p in all_params] or (all_params[:1] if all_params else [])
    
        if not params:
            messagebox.showwarning("Comparison", "No parameters selected / available to plot.")
            return
    
        # ---------------------------
        # Plot
        # ---------------------------
        self.cmp_ax.clear()
        self.cmp_ax.grid(True)
    
        any_plotted = False
        all_cycles = set()
    
        for param in params:
            for fn, cyc_map in by_file.items():
                cycles = sorted(cyc_map.keys())
                if not cycles:
                    continue
    
                xs = []
                ys = []
                for cy in cycles:
                    v = cyc_map[cy].get(param, np.nan)
                    try:
                        v = float(v)
                    except Exception:
                        v = np.nan
                    xs.append(int(cy))
                    ys.append(v)
    
                xs = np.asarray(xs, dtype=int)
                ys = np.asarray(ys, dtype=float)
    
                # Skip series that has no finite values
                if not np.any(np.isfinite(ys)):
                    continue
    
                self.cmp_ax.plot(xs, ys, marker='o', linewidth=1.8, label=f"{fn} | {param}")
                any_plotted = True
                all_cycles.update(xs.tolist())
    
        if not any_plotted:
            self.cmp_ax.text(
                0.5, 0.5,
                "Nothing to plot.\n(Selected params exist but are all NaN for the chosen files/cycles.)",
                ha="center", va="center", transform=self.cmp_ax.transAxes
            )
            self.cmp_ax.set_xlabel("cycle")
            self.cmp_ax.set_ylabel("value")
            try:
                self.cmp_canvas.draw()
            except Exception:
                pass
            return
    
        # Nice labels
        self.cmp_ax.set_xlabel("cycle")
        if len(params) == 1:
            self.cmp_ax.set_ylabel(params[0])
            self.cmp_ax.set_title(f"Comparison: {params[0]} vs cycle")
        else:
            self.cmp_ax.set_ylabel("value")
            self.cmp_ax.set_title("Comparison: selected params vs cycle")
    
        # X ticks (optional: keep readable)
        try:
            if len(all_cycles) <= 25:
                self.cmp_ax.set_xticks(sorted(all_cycles))
        except Exception:
            pass
    
        # Legend
        try:
            self.cmp_ax.legend(loc="best", fontsize=8)
        except Exception:
            pass
    
        # Draw
        try:
            self.cmp_canvas.draw()
        except Exception:
            pass



    # =================== Save PNG ===================
    def save_png(self):
        fn = filedialog.asksaveasfilename(defaultextension='.png')
        if fn:
            self.ax.get_figure().savefig(fn, dpi=300, bbox_inches='tight')
            messagebox.showinfo("Saved", fn)


# =================== RUN ===================
if __name__ == '__main__':
    os.makedirs(EXPORT_FOLDER, exist_ok=True)
    root = tk.Tk()
    app = CVApp(root)
    root.mainloop()
