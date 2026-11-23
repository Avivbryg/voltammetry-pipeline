# Metric families: spec for the next CV Analyzer upgrade

Status: planned. Written as the starting point for the work in Claude Code.

## Goal

Today CV Analyzer computes one fixed set of about 70 parameters per cycle, with
the potentials, windows and onset definitions written into the code. The
upgrade turns every kind of parameter into a **metric family** that the user
can switch on or off and configure in the program, and adds new families.

In the program, the user picks what to analyze (families, variants, settings),
saves the choice as a preset, and runs the analysis. The export contains only
the chosen columns, plus a record of the settings that produced them.

## Ground rules

1. **Keep the onset types, not the fixed values.** The existing onset
   *types* (percentile, curvature k·σ) stay, but their values (which
   percentiles, which k) become settings the user chooses in the program,
   like every new onset type. The 16 fixed definitions are no longer a
   requirement. Other existing columns (peaks, fixed potentials, junction,
   charge) keep their names and meaning.
2. **The export columns are the contract with CV Explorer.** A column name
   says what was measured, in a fixed pattern (below), so the Explorer can
   group columns without a hard-coded list.
3. **One family = one self-contained unit.** Adding a family must not need
   changes anywhere else (analysis loop, export, report, Explorer).
4. **Each step is its own branch and pull request,** and `main` always runs.

## Architecture

### Cycle context

The analysis loop builds one object per cycle and hands it to every enabled
family. It holds what the families share, so nothing is computed twice:

- `E`, `I_raw`, `I_corr` (OER-corrected, may be missing), `t`, time-ordered
- forward / backward split (`split_forward_backward`)
- scan rate (detected), cycle number, file name
- the user's windows: anodic peak window, cathodic peak window, junction
  range (today `an_low/an_high`, `cath_low/cath_high`, `JUNCTION_RANGE`)
- smoothing settings (today `_maybe_smooth`)
- cached results that other families use (e.g. the peaks, for E½ or for
  integration limits)

### Family definition

```python
@dataclass
class MetricFamily:
    key: str                 # "onset", "fixedE", "peaks", ...
    label: str               # shown in the UI
    description: str         # one or two lines for the UI and the guide
    curves: tuple            # ("raw",), or ("raw", "corr")
    variants: list           # [Variant(key, label, params={...}, default_on=bool)]
    compute: Callable        # compute(ctx, variant, curve) -> {column: value}
    columns: Callable        # columns(variant, curve) -> [column names], for empty rows
```

A registry (`METRIC_FAMILIES = [...]`) lists the families in display order.
The analysis loop is then just: for each enabled family, variant and curve,
`row.update(family.compute(ctx, variant, curve))`. A family that fails on a
cycle returns NaN for its columns and adds a line to the report, as today.

### Column names

New columns follow `family_variant_detail_curve`:

- `onset_thr0p5mA_raw`, `onset_tangent_corr`
- `I_fwd_1p45_raw` (new-style name of an existing quantity: see legacy below)
- `Qa_lin_raw`, `Qc_lin_raw`, `QcQa_lin_raw`

Potentials and currents in names use the existing `_ekey` style
(1.375 → `1p375`). Units stay as today: V vs RHE, mA, mA·s.

**Legacy names.** Existing columns keep their current names exactly
(`anodic_E_raw`, `raw_I_fwd_1p6`, `onset_raw_E_k5`, ...). Each family can
declare a `legacy_names` map, used for the existing variants only.

### Presets and settings record

- Presets are saved in the existing settings file (`~/.cv_overlay_settings.json`,
  new key `metric_presets`). Two built-in presets:
  - **Default**: today's families with today's values as defaults (used for
    the step 1 comparison), which the user can then change freely.
  - **Quick**: peaks, charge and a few fixed potentials.
- Every export writes a sidecar file next to the CSV
  (`<name>.settings.json`) with the families, variants, parameters, windows,
  smoothing and Analyzer version that produced it.

### UI

A new "Metrics" panel (or dialog from the analysis tab):

- tree with checkboxes: family → variant, raw / corr toggles per family
- selecting a variant shows its parameters (list of potentials, threshold,
  window, baseline type) in an editable form, with the family description
- preset dropdown + Save as preset
- a live count of how many columns the selection will export

The existing "Parameter guide" text is generated from the family
descriptions, so it stays in sync.

## Families

### Existing (refactor only, same output)

| Family | Variants today | Columns today |
|---|---|---|
| Onset: percentile | p5, p10, p15, p25, threshold = I_min + p·(I_max − I_min) inside the anodic window, forward scan | `onset_raw_E_p*`, `onset_corr_E_p*` |
| Onset: curvature (k·σ) | k = 2, 3, 4, 5, 8, 12, 15, 20 on the integral of positive d²I/dE², baseline = first 20 % of points, `min_run` from the UI | `onset_raw_E_k*`, `onset_corr_E_k*` |
| Fixed potential | forward, backward, forward − backward at `TARGET_E` = 1.0 … 1.6 V | `raw_I_fwd_*`, `raw_I_bwd_*`, `raw_I_for_back_*` |
| Junction | E where I_fwd ≈ I_bwd inside `JUNCTION_RANGE` (1.30–1.37 V), with the currents there | `raw_junction_E`, `raw_junction_I_fwd`, `raw_junction_I_bwd`, `raw_junction_for_back` |
| Peaks | anodic max / cathodic min in their windows; ratio; ΔEp | `anodic_E_*`, `anodic_I_*`, `cathodic_E_*`, `cathodic_I_*`, `I_ratio_*`, `dEp_*` |
| Charge | net charge over the cycle | `Q_end_raw`, `Q_end_corr` |

The fixed-potential list becomes a user setting (default = today's list).
For the onset types, the percentiles and the k values become user settings
too (a list typed in the Metrics panel); today's values are only the
defaults.

### New

**Onset (more definitions)**

Every onset type is a variant with its own user settings (thresholds, n·σ,
windows), and the user can add several values of the same type. Column names
carry the setting, e.g. `onset_p10_raw`, `onset_k5_raw`, `onset_abs0p5mA_raw`.

- *Absolute threshold*: first E on the forward scan where I exceeds a set
  current (e.g. 0.1, 0.5, 1 mA). Simple and comparable between samples.
- *Threshold above baseline*: fit a line to a pre-peak window, onset = first
  E where I − baseline > n·σ of the baseline residuals.
- *Tangent intersection*: intersection of the baseline line with the tangent
  at the steepest point of the rising edge (max dI/dE). The textbook
  definition.
- *Derivative maximum*: E of max dI/dE and of max d²I/dE² in the window.
- Option for all onset variants: anodic window (the Ni²⁺/Ni³⁺ peak, as
  today) or a separate OER/UOR window, so the onset of catalysis can be
  defined separately from the onset of the redox peak.

**Peaks (more)**
- *E½* = (Epa + Epc) / 2: the formal-potential estimate of the couple; often
  steadier than either peak.
- *FWHM* of each peak (after a linear baseline).
- *Edge flag*: 1 if the found peak sits at the edge of its window, which
  means the real peak is probably outside the window.

**Peak integration**
- Anodic charge Qa and cathodic charge Qc over each peak's window, after
  subtracting the **capacitive baseline** (the double-layer current measured
  in a non-faradaic window by the capacitive family; forward branch for Qa,
  backward branch for Qc). This is the default and the main method. "None"
  and "linear between the window ends" can remain as options for checking.
- `QcQa` = |Qc| / Qa: reversibility of the redox couple / how much of the
  oxidized charge comes back.
- Not planned: electrons per Ni. All electrodes carry the same amount of
  active material and the comparison is relative, so absolute values per Ni
  are not needed.

**Capacitive / double layer**
- Current at a set potential in a non-faradaic window (forward and backward),
  and half their difference.
- If a file contains several scan rates: Cdl from the slope of that
  half-difference against scan rate, as an ECSA proxy.

**OER / UOR activity**
- *Potential at fixed current*: E on the forward scan where I reaches a set
  value (the inverse of "current at fixed potential"; the usual way activity
  is compared in the literature).
- *Tafel slope from the scan*: slope of E vs log I in a user window on the
  forward scan, with R².
- *Current at the top of the scan* (max E), so a scan that stops at 1.6 V
  still gives a value.

**Normalization (low priority)**
- All electrodes have the same amount of active material, so per-area or
  per-mass columns are not needed now. Left out of the plan; the family
  structure allows adding it later.
- Possibly useful later: normalization by the redox-peak charge (activity
  per electrochemically active Ni).

**Quality**
- Detected scan rate, E range actually reached, number of points.
- Flags: requested potential outside the scanned range (why the 1.6 V columns
  were empty in the old demo), peak at window edge, too few points in a
  window.

### Better in CV Explorer than in the Analyzer

Quantities that compare cycles or files are easier in the Explorer, from the
per-cycle columns:

- value relative to cycle N or to the maximum (retention, %)
- drift per 100 cycles (slope of a parameter against cycle number)
- difference between two experiments (e.g. with and without urea)

## Work plan

Each step is a separate branch and pull request.

1. **Refactor, no new output.** Cycle context + registry + the existing
   families. Acceptance: on the demo set (`sample_data/make_demo_dataset.py`)
   and on two real files, the exported CSV is identical to the one from
   `main` before the change (same columns, same order, values equal within
   1e-12). Keep that comparison as a script, `tests/compare_exports.py`.
   This check is only to prove the refactor changed nothing; after it,
   the onset values are free to change.
2. **Metrics panel + presets + settings sidecar.** Acceptance: the "Default"
   preset reproduces step 1's output; turning a family off removes exactly
   its columns.
3. **New families, one per pull request,** in this order: onset (user
   settings for the existing types, then the new definitions), capacitive
   (needed first, as the baseline for integration), peak integration,
   peaks (E½, FWHM, edge flag), OER/UOR activity, quality. Each comes with a
   short entry in the parameter guide and a check on the demo data.
4. **Explorer reads the families.** Build the parameter tree from the column
   names (family prefix → variant) and show the sidecar settings for a
   loaded file. `PARAM_TREE` stays only as a fallback for old exports.

## Decisions (Aviv, 2026-11-23)

- Onset: keep the types; the user sets the values in the program. The 16
  fixed definitions do not need to be kept as such.
- Normalization: not needed now. All electrodes carry the same amount of
  active material and the comparison is relative.
- Peak-integration baseline: the capacitive current.

## Open questions

- Non-faradaic window for the capacitive current in the α-Ni(OH)₂ data
  (default proposal: 1.00–1.10 V vs RHE, below the redox peaks), and
  whether one window serves both the forward and the backward branch.
