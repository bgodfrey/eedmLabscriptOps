"""Cross-check of the 2026-10-17 line-centre scan (0003-0013) against slow lock drift.

Usage: python scan_drift_check.py <helper dir> <dark mV>

The 0 mm visits fell 4.5 -> 4.0 -> 3.3% over the run. Three corrections, each followed by the same
iris-averaged Gaussian fit (|x| <= 8 mm):
  raw        - no correction
  pmt        - D * mean(F/L~)/(F/L~), as in iris_atom_scan.py (under-corrects: the PMT is
               compressed by dead time at ~11 Mcps)
  reference  - D / (0 mm dip interpolated linearly in time between the three visits / their mean)
Also prints the PMT-normalised 0 mm visits, to see how much drift the PMT leaves.
"""
import sys
import warnings
from pathlib import Path
import h5py
import numpy as np
from scipy import optimize
from scipy.special import erf

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DARK = float(sys.argv[2]) * 1e-3
DAY = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\17")
OFFSET, WIDTH = 0.05e-3, 1.51e-3            # this run's common shape (iris_atom_scan.py output)

rows = []
for n in range(3, 14):
    for p in sorted((DAY / f"{n:04d}").glob("*.h5")):
        s = load_shot(p)
        with h5py.File(p, "r") as h5:
            x = float(h5["globals"].attrs["iris_x_mm"])
        pmt = cse.pmt_cloud_fit(s, 6125.0)
        tc = pmt["center_s"]
        k = int(round(s["fs_pda"] / 10e3))
        m = len(s["pda_v"]) // k * k
        tb = s["t_pda"][:m].reshape(-1, k).mean(axis=1)
        vb = s["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
        w = (tb > s["enable_start_s"] + 0.002) & (tb < s["enable_stop_s"] - 0.002)

        def model(t, c0, c1, c2, d):
            u = t - tc
            return (c0 + c1 * u + c2 * u ** 2) * (1 - d * np.exp(-0.5 * ((t - tc - OFFSET) / WIDTH) ** 2))

        pf, _ = optimize.curve_fit(model, tb[w], vb[w], p0=[np.median(vb[w]), 0, 0, 0.03])
        rows.append({"f": n, "x": x, "t": p.stat().st_mtime, "D": pf[3], "L": pf[0],
                     "F": pmt["amp_cps"] * pmt["sigma_s"]})

xs = np.array([r["x"] for r in rows])
t = np.array([r["t"] for r in rows])
t -= t.min()
D = np.array([r["D"] for r in rows])
L = np.array([r["L"] for r in rows])
F = np.array([r["F"] for r in rows])
pos = np.unique(xs)
lmean = {x: L[xs == x].mean() for x in pos}
q = F / (L / np.array([lmean[x] for x in xs]))
d_pmt = D * q.mean() / q
# Reference: 0 mm visits, one point per folder at its mean time
ref_f = sorted({r["f"] for r in rows if r["x"] == 0})
fold = np.array([r["f"] for r in rows])
rt = np.array([t[fold == f].mean() for f in ref_f])
rv = np.array([D[fold == f].mean() for f in ref_f])
d_ref = D / (np.interp(t, rt, rv) / rv.mean())
print("0 mm visits (raw / PMT-normalised): " + ", ".join(
    f"{f:04d} at {tt/60:4.1f} min: {100*D[fold == f].mean():.2f} / {100*d_pmt[fold == f].mean():.2f}%"
    for f, tt in zip(ref_f, rt)))


def profile(x, a, xc, sx):
    h = 1.0
    return a * sx * np.sqrt(np.pi / 2) / 2.0 * (erf((x + h - xc) / (np.sqrt(2) * sx)) - erf((x - h - xc) / (np.sqrt(2) * sx)))


for name, d in (("raw", D), ("pmt", d_pmt), ("reference", d_ref)):
    mean = np.array([d[xs == x].mean() for x in pos])
    sem = np.array([d[xs == x].std(ddof=1) / np.sqrt((xs == x).sum()) for x in pos])
    p, c = optimize.curve_fit(profile, pos, mean, p0=(mean.max(), 0, 3), sigma=sem, absolute_sigma=True)
    chi2 = np.sum(((mean - profile(pos, *p)) / sem) ** 2)
    e = np.sqrt(np.diag(c) * max(1, chi2 / (len(pos) - 3)))
    print(f"{name:9s}: chi2/dof {chi2:5.1f}/{len(pos)-3}, peak {100*p[0]:.2f} +- {100*e[0]:.2f}%, "
          f"centre {p[1]:+.2f} +- {e[1]:.2f} mm, sigma_x {abs(p[2]):.2f} +- {e[2]:.2f} mm")
