"""
Generate a large synthetic CV dataset in MPR Converter CSV format.

The three files in sample_data/ have 3 cycles each, enough for the tutorial.
This script writes a long activation experiment instead, so the pipeline can be
tried at the scale it was built for: hundreds of cycles per file.

Every number here is invented. The curves imitate a Ni(OH)2-like redox couple
vs RHE at 10 mV/s with an OER tail, and nothing is taken from real measurements.

Two slow changes are built in, so the CV Explorer has something to show:
- the cathodic peak drifts to higher potential over the run (the kind of shift
  that can serve as a proxy for an alpha -> beta phase change), and
- the current near the top of the scan (e.g. at 1.6 V) falls after activation.
How large each change is differs between the three anions. The sizes are
chosen for the demo and say nothing about real anions.

Usage:
    python sample_data/make_demo_dataset.py                  # 3 files x 300 cycles
    python sample_data/make_demo_dataset.py --cycles 100 --out my_folder

The output columns and number format match what MPR Converter
(https://mpr-converter.vercel.app) writes for a .mpr file, so CV Analyzer reads
these files the same way it reads converter output. CSV files are ignored by
git (.gitignore), so the generated data never ends up in the repository.
"""

import argparse
import os

import numpy as np

# name, anodic peak (V), cathodic peak (V), final peak height (mA),
# activation time constant (cycles), anodic peak drift over the run (V),
# cathodic peak shift (V), loss of high-potential activity (fraction),
# time constant of both slow changes (cycles)
ANIONS = [
    ("Nitrate",   1.420, 1.330, 1.30, 60.0,  0.012, 0.022, 0.25, 160.0),
    ("Carbonate", 1.438, 1.335, 1.05, 90.0,  0.018, 0.040, 0.40, 120.0),
    ("Chloride",  1.410, 1.322, 1.50, 40.0,  0.008, 0.006, 0.08, 220.0),
]

STEP = 0.002    # V per point
RATE = 0.010    # V/s
E_LOW, E_HIGH = 1.00, 1.62   # top just above 1.6 V, so the current at 1.6 V is defined in every cycle
C_DL = 1.5      # double-layer capacitance (mF): capacitive current = C * scan rate

HEADER = "mode,time/s,control/V,Ewe/V,<I>/mA,cycle number"


def peak(E, center, width, skew):
    """Asymmetric peak: a Gaussian with a slower tail on one side (skew > 0:
    tail toward higher potential, as on a real anodic wave)."""
    x = (E - center) / width
    tail = np.where(x * skew > 0, 1.0 / (1.0 + abs(skew)), 1.0)
    return np.exp(-(x * tail) ** 2)


def one_file(path, anion, n_cycles, rng):
    name, Ea0, Ec0, height, tau, drift, c_shift, loss, tau_slow = anion
    up = np.arange(E_LOW, E_HIGH + 1e-9, STEP)
    down = up[-2::-1]
    E = np.concatenate([up, down])
    forward = np.arange(len(E)) < len(up)
    n = len(E)
    dt = STEP / RATE

    i_cap = C_DL * RATE                       # mA, sign follows the sweep
    rows = []
    t0 = 0.0
    for cyc in range(1, n_cycles + 1):
        act = 1.0 - 0.75 * np.exp(-(cyc - 1) / tau)       # activation: grows, then levels off
        Ea = Ea0 + drift * (1.0 - np.exp(-(cyc - 1) / (2 * tau)))
        slow = 1.0 - np.exp(-(cyc - 1) / tau_slow)          # 0 -> 1 over the run
        Ec = Ec0 + c_shift * slow                           # cathodic peak moves up
        decay = 1.0 - loss * slow                           # high-potential activity falls
        h = height * act * (1.0 + rng.normal(0, 0.01))     # small cycle-to-cycle scatter

        I = np.where(forward,
                     i_cap + h * peak(E, Ea, 0.024, +0.6),
                     -i_cap - 0.78 * h * peak(E, Ec, 0.022, -0.6))
        I += 0.008 * (1.0 + 0.3 * act) * decay * np.exp((E - 1.50) / 0.028)   # OER tail
        I += 0.004 * np.sin(np.pi * (E - E_LOW) / (E_HIGH - E_LOW))   # sloping baseline
        I += rng.normal(0, 0.004, n)                                   # noise

        t = t0 + dt * np.arange(n)
        t0 = t[-1] + dt
        Ewe = E + rng.normal(0, 0.0002, n)                             # potentiostat jitter
        block = np.column_stack([np.full(n, 2), t, E, Ewe, I, np.full(n, cyc)])
        rows.append(block)

    table = np.vstack(rows)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(HEADER + "\n")
        np.savetxt(fh, table, delimiter=",",
                   fmt=["%d", "%.10g", "%.10g", "%.7g", "%.7g", "%d"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cycles", type=int, default=300, help="cycles per file (default 300)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo_large"),
                    help="output folder (default sample_data/demo_large)")
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    for k, anion in enumerate(ANIONS, start=1):
        fname = f"ACT_{'ABC'[k - 1]}_{20000 + k}_{anion[0]}_demo.csv"
        path = os.path.join(args.out, fname)
        one_file(path, anion, args.cycles, rng)
        size = os.path.getsize(path) / 1e6
        print(f"{fname}: {args.cycles} cycles, {size:.1f} MB")
    print(f"\nWritten to {args.out}\nOpen this folder in CV Analyzer.")


if __name__ == "__main__":
    main()
