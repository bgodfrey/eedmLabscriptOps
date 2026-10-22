"""Is the leading-edge tail in the vertical dip profile (10/17 0 mm shots) real or a baseline artifact?

Usage: python vertical_tail_check.py <helper dir> <dark mV>
Compares, for the 0 mm shots of 0003/0008/0013:
 - dip profile with baselines from different far windows (quadratic |u| > 6 ms whole window;
   quadratic 6-15 ms; linear 6-12 ms);
 - the mean PMT rate envelope (background-subtracted, peak-normalised), which has no light baseline.
Prints the profile at +-4.6, +-6.9, +-9.2 mm (u = -+2, -+3, -+4 ms) for each.
"""
import sys
import warnings
from pathlib import Path
import numpy as np
from scipy import signal

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DARK = float(sys.argv[2]) * 1e-3
DAY = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\17")
V = 2.30
grid = np.arange(-0.012, 0.012 + 1e-9, 2e-4)
prof = {"quad |u|>6": [], "quad 6-15": [], "lin 6-12": []}
pmt_tr = []
for f in ("0003", "0008", "0013"):
    for p in sorted((DAY / f).glob("*.h5")):
        s = load_shot(p)
        pm = cse.pmt_cloud_fit(s, 6125.0)
        tc = pm["center_s"]
        k = int(round(s["fs_pda"] / 10e3))
        m = len(s["pda_v"]) // k * k
        tb = s["t_pda"][:m].reshape(-1, k).mean(axis=1)
        vb = s["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
        w = (tb > s["enable_start_s"] + 0.002) & (tb < s["enable_stop_s"] - 0.002)
        u, v = tb[w] - tc, vb[w]
        for name, lo, hi, deg in (("quad |u|>6", 0.006, 1.0, 2), ("quad 6-15", 0.006, 0.015, 2),
                                  ("lin 6-12", 0.006, 0.012, 1)):
            sel = (np.abs(u) > lo) & (np.abs(u) < hi)
            prof[name].append(np.interp(grid, u, 1 - v / np.polyval(np.polyfit(u[sel], v[sel], deg), u)))
        # PMT envelope: smoothed rate, minus the fitted offset, over the fitted amplitude
        kern, _ = cse.triangular_kernel(s["fs_pmt"], 6125.0, 2)
        r = signal.fftconvolve(s["pmt_rate"], kern, mode="same")
        pmt_tr.append(np.interp(grid, s["t_pmt"] - tc, (r - pm["offset_cps"]) / pm["amp_cps"]))
z = -1e3 * V * grid
probe = [2, 3, 4]
print("u ms:      " + "  ".join(f"{-t:+d}ms(z{V*t:+.1f})  {t:+d}ms(z{-V*t:+.1f})" for t in probe))
for name, tr in prof.items():
    mtr = np.mean(tr, axis=0)
    pk = mtr[np.abs(grid) < 1e-3].max()
    vals = [(np.interp(-t * 1e-3, grid, mtr) / pk, np.interp(t * 1e-3, grid, mtr) / pk) for t in probe]
    print(f"{name:11s} " + "  ".join(f"   {a:5.3f}           {b:5.3f}     " for a, b in vals)
          + f"  floor at -10/+10 ms: {100*np.interp(-0.010, grid, mtr):+.3f}/{100*np.interp(0.010, grid, mtr):+.3f}%")
mp = np.mean(pmt_tr, axis=0)
vals = [(np.interp(-t * 1e-3, grid, mp), np.interp(t * 1e-3, grid, mp)) for t in probe]
print(f"{'PMT':11s} " + "  ".join(f"   {a:5.3f}           {b:5.3f}     " for a, b in vals))
g = [np.exp(-0.5 * (t / 1.40) ** 2) for t in probe]
print(f"{'Gauss 1.40':11s} " + "  ".join(f"   {a:5.3f}           {a:5.3f}     " for a in g))
