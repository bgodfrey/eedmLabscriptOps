"""Which functional form fits the vertical (transit-time) dip profile of the 10/17 0 mm shots?

Usage: python vertical_model_compare.py <helper dir> <dark mV>
Profile: fractional dip vs z = v (t_c - t), 0.2 ms bins, mean +- sem over the 0 mm shots of
0003/0008/0013 (no PMT scaling: shape only). Two baselines: quadratic fitted at |u| > 6 ms (as in the
figure) and at |u| > 10 ms (to see whether the tail runs into the 6 ms baseline region).
Models (tail toward +z = early arrival):
  gauss   A G(z; mu, s)
  split   two-sided Gaussian (s below / above the peak)       [the current figure]
  emg     exponentially modified Gaussian: Gaussian (*) one-sided exponential toward +z
  2gauss  core Gaussian + second Gaussian (free centre and width)
chi2 uses the per-bin sem, which ignores the bin-to-bin correlation, so compare chi2 between models
rather than reading it absolutely.
"""
import sys
import warnings
from pathlib import Path
import numpy as np
from scipy import optimize
from scipy.stats import exponnorm

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DARK = float(sys.argv[2]) * 1e-3
DAY = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\17")
V = 2.30
edges = np.arange(-0.012, 0.012 + 1e-9, 2e-4)
cen = 0.5 * (edges[1:] + edges[:-1])
z = -1e3 * V * cen
traces = {6: [], 10: []}
for f in ("0003", "0008", "0013"):
    for p in sorted((DAY / f).glob("*.h5")):
        s = load_shot(p)
        pm = cse.pmt_cloud_fit(s, 6125.0)
        k = int(round(s["fs_pda"] / 10e3))
        m = len(s["pda_v"]) // k * k
        tb = s["t_pda"][:m].reshape(-1, k).mean(axis=1)
        vb = s["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
        w = (tb > s["enable_start_s"] + 0.002) & (tb < s["enable_stop_s"] - 0.002)
        u, v = tb[w] - pm["center_s"], vb[w]
        idx = np.digitize(u, edges) - 1
        for lo in traces:
            far = np.abs(u) > lo * 1e-3
            fr = 1 - v / np.polyval(np.polyfit(u[far], v[far], 2), u)
            traces[lo].append([fr[idx == i].mean() for i in range(len(cen))])


def gauss(z, a, mu, s):
    return a * np.exp(-0.5 * ((z - mu) / s) ** 2)


def split(z, a, mu, s1, s2):
    return a * np.exp(-0.5 * ((z - mu) / np.where(z < mu, s1, s2)) ** 2)


def emg(z, area, mu, s, tau):
    return area * exponnorm.pdf(z, abs(tau) / abs(s), loc=mu, scale=abs(s))


def two(z, a1, m1, s1, a2, m2, s2):
    return gauss(z, a1, m1, s1) + gauss(z, a2, m2, s2)


MODELS = {"gauss": (gauss, lambda y: (y.max(), 0, 3)),
          "split": (split, lambda y: (y.max(), 0, 3, 3)),
          "emg": (emg, lambda y: (y.sum() * 0.46, -1, 2.5, 2)),
          "2gauss": (two, lambda y: (y.max(), -0.5, 3, 0.1 * y.max(), 5, 5))}
for lo, zmax in ((6, 13.0), (10, 22.0)):
    tr = np.array(traces[lo])
    y, e = tr.mean(axis=0), tr.std(axis=0, ddof=1) / np.sqrt(len(tr))
    sel = np.abs(z) <= zmax
    print(f"\nbaseline |u| > {lo} ms, fit |z| <= {zmax:.0f} mm, {sel.sum()} bins, {len(tr)} shots")
    print("   profile at z = +9/+12/+15/+18/+21 mm: " + ", ".join(
        f"{100*np.interp(zz, z[::-1], y[::-1]):.2f}%" for zz in (9, 12, 15, 18, 21))
          + f";  at -9/-12: {100*np.interp(-9, z[::-1], y[::-1]):.2f}%, {100*np.interp(-12, z[::-1], y[::-1]):.2f}%")
    for name, (fn, p0) in MODELS.items():
        try:
            p, c = optimize.curve_fit(fn, z[sel], y[sel], p0=p0(y[sel]), sigma=e[sel], maxfev=20000)
        except RuntimeError as err:
            print(f"   {name:7s} failed: {err}")
            continue
        r = (y[sel] - fn(z[sel], *p)) / e[sel]
        chi2, dof = np.sum(r ** 2), sel.sum() - len(p)
        peak_i = np.argmax(fn(z, *p))
        extra = ""
        if name == "emg":
            extra = f"  gauss s {abs(p[2]):.2f} mm, exp tail length {abs(p[3]):.2f} mm, mu {p[1]:+.2f}"
        if name == "2gauss":
            ar1, ar2 = abs(p[0] * p[2]), abs(p[3] * p[5])
            extra = (f"  core s {abs(p[2]):.2f} mm at {p[1]:+.2f}; second s {abs(p[5]):.2f} mm at {p[4]:+.2f}, "
                     f"{ar2/(ar1+ar2):.0%} of the area")
        if name == "split":
            extra = f"  s below {abs(p[2]):.2f}, above {abs(p[3]):.2f} mm"
        if name == "gauss":
            extra = f"  s {abs(p[2]):.2f} mm"
        print(f"   {name:7s} chi2/dof {chi2:7.1f}/{dof}  peak {100*fn(z, *p).max():.2f}% at {z[peak_i]:+.1f} mm "
              f"(data {100*y.max():.2f}% at {z[np.argmax(y)]:+.1f}){extra}")
