"""Pass 1 of iris_atom_scan.py (free centre and width per shot), 10/17 0003-0013, to see which shots the
8-sigma cut removes and whether weighting the widths by their errors changes the common width.
Same preprocessing and model as iris_atom_scan.py.
Usage: python shape_width_weighting.py <helper dir>
"""
import sys
from pathlib import Path
import h5py
import numpy as np
from scipy import optimize

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse

DAY = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\17")
DARK = 0.0527


def dip_model(t, tc, c0, c1, c2, depth, centre, width):
    u = t - tc
    return (c0 + c1 * u + c2 * u ** 2) * (1 - depth * np.exp(-0.5 * ((t - centre) / width) ** 2))


rows = []
for n in range(3, 14):
    for p in sorted((DAY / f"{n:04d}").glob("*.h5")):
        shot = load_shot(p)
        with h5py.File(p, "r") as h5:
            x = float(h5["globals"].attrs["iris_x_mm"])
        pmt = cse.pmt_cloud_fit(shot, shot["mod_freq"] or 6125.0)
        if shot["pda_v"].max() >= 9.95 or not pmt["ok"]:
            continue
        k = int(round(shot["fs_pda"] / 10e3))
        m = len(shot["pda_v"]) // k * k
        tb = shot["t_pda"][:m].reshape(-1, k).mean(axis=1)
        vb = shot["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
        win = (tb > shot["enable_start_s"] + 0.002) & (tb < shot["enable_stop_s"] - 0.002)
        tw, vw, tc = tb[win], vb[win], pmt["center_s"]
        p0 = [np.median(vw), 0, 0, 0.02, tc, 1.3e-3]
        lo = [-np.inf, -np.inf, -np.inf, -0.5, tc - 0.005, 0.3e-3]
        hi = [np.inf, np.inf, np.inf, 0.5, tc + 0.005, 5e-3]
        try:
            pf, cov = optimize.curve_fit(lambda t, *a: dip_model(t, tc, *a), tw, vw, p0=p0, bounds=(lo, hi))
            e = np.sqrt(np.diag(cov))
            rows.append((x, pf[3], e[3], pf[4] - tc, pf[5], e[5]))
        except (RuntimeError, ValueError):
            rows.append((x, np.nan, np.nan, np.nan, np.nan, np.nan))

a = np.array(rows)
x, d, ed, off, w, ew = a.T
clean = (d > 8 * ed) & (w > 0.35e-3) & (w < 4.9e-3)
print(f"{clean.sum()} clean of {len(a)}")
print("  x mm  shots  clean  median D %   median width (clean, ms)  median width error (ms)")
for xx in np.unique(x):
    s = x == xx
    c = s & clean
    print(f"  {xx:+5.1f}  {s.sum():4d}  {c.sum():4d}     {100*np.nanmedian(d[s]):5.2f}        "
          f"{1e3*np.median(w[c]) if c.any() else np.nan:5.2f}                    "
          f"{1e3*np.median(ew[c]) if c.any() else np.nan:5.2f}")
wc, ewc = w[clean], ew[clean]
wt = 1 / ewc ** 2
print(f"median width (as used): {1e3*np.median(wc):.2f} ms; inverse-variance weighted mean: "
      f"{1e3*np.sum(wt*wc)/np.sum(wt):.2f} ms; |x| <= 2 mm only, median: "
      f"{1e3*np.median(w[clean & (np.abs(x) <= 2)]):.2f} ms")
print(f"median centre offset (as used): {1e3*np.median(off[clean]):+.2f} ms; weighted: "
      f"{1e3*np.sum(wt*off[clean])/np.sum(wt):+.2f} ms")
