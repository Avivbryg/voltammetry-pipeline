# Voltammetry Pipeline

Tools that take cyclic voltammetry (CV) data from raw instrument files to comparison between samples, built during my M.Sc. research on electrocatalysts.

An electrochemistry experiment can produce thousands of CV cycles. This pipeline turns them into a consistent dataset, and then into plots that compare parameters across samples, anions and cycle number.

| Stage | Tool | What it does |
|---|---|---|
| 1. Extract | **CV Analyzer** (`analyzer/`) | Reads EC-Lab `.mpt` files, or `.mpr` files converted with [MPR Converter](https://mpr-converter.vercel.app), and turns every cycle into one row of metrics: peaks, charge, OER-corrected values, onsets, fixed-potential currents. Exports CSV. |
| 2. Compare | **CV Explorer** (`explorer/`) | Reads the exported CSV files and plots one parameter against cycle number for many experiments at once, on a Y axis shared by every experiment of the same measurement type, so anions can be compared directly. Graphs can be arranged on a printable appendix page. |

![CV Analyzer main window: 900 synthetic cycles from three files, analyzed in about 10 s](docs/screenshot_main.png)

## Stage 1: CV Analyzer

### The problem

An electrocatalysis experiment produces many `.mpt` files, each with dozens to hundreds of CV cycles. Extracting the same numbers from every cycle (peak potentials and currents, charge, currents at fixed potentials, onset potentials) by hand in Excel or EC-Lab is slow, hard to keep consistent between people, and easy to get wrong.

Some features are also hard to reduce to a single number. An onset potential, for example, depends on how it is defined. For these, a visual side-by-side comparison of the curves is more honest than a column in a table.

### What it does

CV Analyzer turns that into a repeatable workflow, the first stage of a pipeline from raw instrument files to comparison:

**Open a folder → tick cycles → Run Analysis → Export CSV**

- **File and cycle browser**: every `.mpt` in a folder tree, with each cycle and its scan rate. Files are parsed lazily, so large folders open instantly.
- **Binary `.mpr` files via MPR Converter**: CSV or TXT files exported by [MPR Converter](https://mpr-converter.vercel.app) are listed and analyzed like `.mpt` files, so there is no need to export from EC-Lab first. Other CSV files in the folder, including this program's own exports, are left out of the tree.
- **Reference conversion**: potentials converted to RHE from common reference electrodes and the electrolyte pH.
- **Peak analysis**: anodic and cathodic peak potentials and currents within configurable windows.
- **Charge integration**: a Q vs E plot alongside I vs E, I vs t and E vs t.
- **OER background subtraction**: a corrected curve and corrected metrics, with a validity flag per cycle.
- **Onset potentials and derivatives**: first and second derivatives with live smoothing (moving average or Savitzky–Golay).
- **Currents at fixed potentials**: forward, backward and forward-minus-backward currents at set potentials.
- **Analysis report**: statistics for every parameter (n, missing values, mean, std, variance), with missing values flagged.
- **Batch processing**: *Find unanalyzed…* steps through every file without an export yet.
- **Overlay window**: plots the same cycle from many files on a common Y axis, matched by name (for example by anion and measurement type). It is built for features that resist a fixed definition, such as onsets, where comparing the curves themselves tells more than a number.
- **Appendix page builder**: arranges graphs on a printable page and saves it as PNG or PDF.
- **Guided tutorial (English / Hebrew)**: runs on first launch on synthetic sample files, and every step waits until the user has actually done it.

![Overlay window comparing cycles 10 and 300 across three anions (synthetic data)](docs/screenshot_overlay.png)

### Quick start

Requires Python 3.9+ with Tkinter.

```bash
git clone https://github.com/Avivbryg/voltammetry-pipeline.git
cd voltammetry-pipeline
pip install -r requirements.txt
python analyzer/cv_analyzer.py
```

On first launch the tutorial starts by itself. To try the program on the included data instead, click **Open Folder** and choose `sample_data/`.

> **Linux:** Tkinter may need to be installed separately, for example `sudo apt install python3-tk`.
> SciPy is optional; without it, smoothing falls back to a moving average.

### Input files

- **EC-Lab ASCII** (`.mpt`), as exported by EC-Lab.
- **MPR Converter output** (`.csv` or `.txt`): open the binary `.mpr` in [MPR Converter](https://mpr-converter.vercel.app), export CSV or TXT, and put the file in the data folder. The converter keeps EC-Lab's column names (`Ewe/V`, `<I>/mA`, `cycle number`, `time/s`), which is what CV Analyzer looks for. Its XLSX output is not read; use CSV or TXT.

### Sample data

`sample_data/` contains three **synthetic** CV files in EC-Lab ASCII format: a Ni(OH)₂-like redox couple vs RHE at 10 mV/s, 3 cycles each, with an OER tail above 1.5 V. They are generated by the program itself (`make_tutorial_samples`), so they contain no real experimental data.

For a run at a realistic scale, `sample_data/make_demo_dataset.py` writes three synthetic activation experiments of 300 cycles each (about 7 MB per file) in MPR Converter's CSV format. Two slow changes are built in, with a different size for each anion, so that the CV Explorer has something to show: the cathodic peak drifts to higher potential, and the current at 1.6 V falls after activation. The sizes are invented and say nothing about real anions.

```bash
python sample_data/make_demo_dataset.py            # writes sample_data/demo_large/
```

Then open `sample_data/demo_large/` in CV Analyzer. These files are also synthetic, and CSV files are kept out of git.

File names follow the scheme `<measurement>_<station>_<electrode>_<anion>_<…>.mpt`. The export dialog and the overlay window read these parts from the name; any other name still works for analysis.

### Naming settings

The names in that scheme differ from lab to lab, so they are set in the program, under **File → Naming settings…**:

- **Measurement types**, one per line as `code | name | wording in graph titles`, for example `act | ACT-night | activation, no additive`. The code is what appears in file names.
- **Stations**, **anions** and **laborants**, one per line.

The lists are saved in `~/.cv_overlay_settings.json`. The defaults (`ACT-night`, `STB-day`, stations A–D, and the laborants Maya R., Tom K. and Noa S.) are invented examples that match the synthetic sample files.

### Output

**Export Analysis CSV** writes one row per cycle with all computed metrics: raw and OER-corrected peak values, charges, onsets, and fixed-potential currents. The default export folder is `~/CV_Analyzer_exports`. Both the data folder and the export folder can be changed under **File → Set default folders…**, and the choice is remembered in `~/.cv_overlay_settings.json`.

### Background

I built this tool during my M.Sc. in Chemistry at the Technion, while studying how intercalated anions affect α-Ni(OH)₂ activity toward the urea oxidation reaction. It grew over many versions from a plotting script into the tool I used to analyze all of my CV data, hundreds to thousands of cycles per experiment. The guided tutorial and folder setup were added so that another researcher can pick it up and run it without help.

Built with Python, Tkinter, NumPy, Matplotlib and SciPy.

## Stage 2: CV Explorer

### The problem

After stage 1 there is one CSV per experiment, with one row per cycle and about 70 parameters per row. The question is then which of these parameters actually separates the samples, for example which ones follow the stability of the phase over hundreds of cycles. Answering that means looking at the same parameter for many experiments at once, again and again, on axes that can be compared.

### What it does

**Load the export folder → pick a parameter → tick experiments (or a whole anion) → compare**

- **Experiments vs Cycles** (the main view): one parameter against cycle number for many experiments, one colour per anion. Clicking an anion selects all of its experiments, and the plot updates as the selection changes.
- **Shared Y range**: by default the Y axis spans the parameter's range over every loaded experiment of the same measurement type. Graphs of different anions then sit on the same scale and can be placed side by side.
- **Metadata from the file name**: anion, station, electrode, measurement type and laborant are read from the name that CV Analyzer's export dialog writes (`Anion_Station_Electrode_Measurement_Laborant.csv`). The measurement wording in graph titles comes from CV Analyzer's naming settings.
- **Parameter tree with favorites**: parameters grouped by kind (onsets, peaks, fixed-potential currents, charge; raw and OER-corrected). Columns the tree does not know are listed under *Other columns*.
- **Appendix page**: *→ Appendix page* sends the current graph to a page layout (rows × columns, A4 or Letter), where graphs can be reordered and retitled, and saved as PNG, PDF or TIFF at up to 1200 DPI.
- **Experimental tabs**: *Params vs Cycles* (several parameters of one experiment), *Scatter*, *Correlation Matrix* and *Target Correlations*. They work, but they were not part of the analysis I relied on.

### Quick start

```bash
pip install -r requirements.txt
python explorer/cv_metrics_explorer.py
```

Click **Load a folder…** and choose the folder of CSV files exported by CV Analyzer (by default `~/CV_Analyzer_exports`). The last folder is loaded again on the next start.

To try it on synthetic data: run `sample_data/make_demo_dataset.py`, open `sample_data/demo_large/` in CV Analyzer, analyze and export each file, then load the export folder in CV Explorer. In the demo, `cathodic_E_raw` and `raw_I_fwd_1p6` under *Experiments vs Cycles* show the two built-in trends.

Favorites and the last folder are saved next to the script (`cv_explorer_favorites.json`, `cv_explorer_config.json`).

### Units

Potentials are in V (vs RHE, as exported), currents in mA and charges in mA·s (= mC), as in the CV Analyzer export.

### Background

I used CV Explorer to go through the parameters of many samples in my M.Sc. work, looking for the ones that track phase stability during CV cycling. Two stood out: the cathodic peak potential, as a proxy for the α → β transition of Ni(OH)₂, and the current at a high anodic potential, which shows the loss of activity over time.

Built with Python, Tkinter, pandas, NumPy and Matplotlib.

## License

MIT. See [LICENSE](LICENSE).
