"""Which 0 mm visits does each position's drift correction use, and with what weights?
Same correction as beam_profile_softedge_fit.py: corrected = raw * mean(visits) / interp(t; visits),
visit time = mean file mtime of that folder's shots. Also compares mtime with the labscript
'run time' (shot start), to see whether later file writes have moved the mtimes.
"""
import os
from datetime import datetime
from pathlib import Path
import h5py
import numpy as np

ROOT = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10")
RUNS = [("10/15", ROOT / "15", 17, 33), ("10/16", ROOT / "16", 23, 35), ("10/17", ROOT / "17", 3, 13)]

for name, day, first, last in RUNS:
    folders = []
    for n in range(first, last + 1):
        mt, rt, x = [], [], None
        for p in sorted((day / f"{n:04d}").glob("*.h5")):
            with h5py.File(p, "r") as h5:
                x = float(h5["globals"].attrs["iris_x_mm"])
                rt.append(datetime.strptime(h5.attrs["run time"], "%Y%m%dT%H%M%S.%f").timestamp())
            mt.append(os.path.getmtime(p))
        folders.append((n, x, np.array(mt), np.array(rt)))
    t0 = min(f[3].min() for f in folders)
    lag = np.concatenate([f[2] - f[3] for f in folders])
    print(f"\n{name}: mtime - run time per shot: min {lag.min():.1f} s, median {np.median(lag):.1f} s, max {lag.max():.1f} s")
    refs = [(n, mt.mean()) for n, x, mt, rt in folders if x == 0]
    ref_t = np.array([t for _, t in refs])
    print("  folder  x mm   t (min, mtime)  bracketing visits (weights)")
    for n, x, mt, rt in folders:
        tm = mt.mean()
        if x == 0:
            print(f"  {n:04d}  {x:+5.1f}   {(tm - t0)/60:6.1f}         (0 mm visit)")
            continue
        k = np.searchsorted(ref_t, tm)
        if k == 0:
            desc = f"before first visit -> {refs[0][0]:04d} only (held flat)"
        elif k == len(ref_t):
            desc = f"after last visit -> {refs[-1][0]:04d} only (held flat)"
        else:
            w = (tm - ref_t[k - 1]) / (ref_t[k] - ref_t[k - 1])
            desc = f"{refs[k-1][0]:04d} ({1-w:.2f}) + {refs[k][0]:04d} ({w:.2f})"
        print(f"  {n:04d}  {x:+5.1f}   {(tm - t0)/60:6.1f}         {desc}")
