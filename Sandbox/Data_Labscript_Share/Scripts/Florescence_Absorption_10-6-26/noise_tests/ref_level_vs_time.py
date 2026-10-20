"""Probe light at the 0 mm reference visits versus time, for the three iris scans on one plot.

Usage: python ref_level_vs_time.py <figure.png>
(Runs: 2026-10-15 0017-0033, 2026-10-16 0023-0035, 2026-10-17 0003-0013, as in
beam_profile_softedge_fit.py.)

These visits are what the drift correction in beam_profile_softedge_fit.py interpolates between.
Per shot: median of the PDA record (TOF_florescence, the whole 0.2 s; the probe is on throughout and
the median ignores the few-ms atom dip) minus the dark level, 52.7 mV.
Per 0 mm folder: the mean of its shot medians, error bar = their standard deviation (ddof 1; 8 shots
per folder on 10/16-17, 4 on 10/15), at the mean shot time.
Time: minutes since the run's first shot (any position), from each shot's labscript 'run time'
attribute (the shot start), not the file modification time, which later analysis can change.
No drift correction here: this is the raw reference level.
"""
import sys
from datetime import datetime
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = Path(sys.argv[1])
ROOT = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10")
DARK = 0.0527                                   # V at 70 dB (10/15 0010, probe blocked)
RUNS = [("10/15", ROOT / "15", 17, 33), ("10/16", ROOT / "16", 23, 35), ("10/17", ROOT / "17", 3, 13)]
# Categorical slots 1-3 of the reference palette, in order; ink, muted grey, surface
COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
INK, MUTED, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"


def run_time(h5):
    """Shot start time from labscript's 'run time' attribute, e.g. '20261017T002702.781890'."""
    return datetime.strptime(h5.attrs["run time"], "%Y%m%dT%H%M%S.%f").timestamp()


fig, ax = plt.subplots(figsize=(8, 4.8), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
ax.tick_params(colors=MUTED)
for spine in ax.spines.values():
    spine.set_color(MUTED)
ax.grid(True, color=MUTED, alpha=0.15)
print("run    folder  t (min)  shots  mean (V)   sd (V)   sd %")
for (name, day, first, last), color in zip(RUNS, COLORS):
    t0, visits = None, []
    for n in range(first, last + 1):
        times, levels, x = [], [], None
        for p in sorted((day / f"{n:04d}").glob("*.h5")):
            with h5py.File(p, "r") as h5:
                x = float(h5["globals"].attrs["iris_x_mm"])
                times.append(run_time(h5))
                levels.append(np.median(h5["data/traces/TOF_florescence"][()]["values"]) - DARK)
        if not times:
            continue
        t0 = min(times) if t0 is None else min(t0, min(times))   # first shot of the run, any position
        if x == 0:
            visits.append((n, np.mean(times), np.mean(levels), np.std(levels, ddof=1), len(levels)))
    t = np.array([(v[1] - t0) / 60 for v in visits])
    m = np.array([v[2] for v in visits])
    s = np.array([v[3] for v in visits])
    for (n, _, mean, sd, k), tt in zip(visits, t):
        print(f"{name}  {n:04d}   {tt:6.1f}   {k:4d}   {mean:7.3f}   {sd:7.3f}  {100*sd/mean:5.1f}")
    ax.plot(t, m, color=color, lw=1.2, alpha=0.6)          # thin guide between visits
    ax.errorbar(t, m, yerr=s, fmt="o", ms=8, color=color, mec=SURFACE, mew=1.5, capsize=4, lw=1.5,
                label=f"{name} ({visits[0][4]} shots/visit)")
ax.set_xlabel("Time since the run's first shot (min)", color=INK)
ax.set_ylabel("Probe light at 0 mm (V)", color=INK)
ax.set_title("0 mm reference level vs time (mean ± sd of shot medians)", color=INK)
ax.set_ylim(bottom=0)
ax.legend(frameon=False, fontsize=9, loc="center right", labelcolor=INK)
fig.tight_layout()
fig.savefig(FIG, dpi=130, facecolor=SURFACE)
print(f"figure: {FIG}")
