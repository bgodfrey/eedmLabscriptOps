"""How far off-centre must the beam be? Chi2 of the plain-Gaussian fit (|x| <= 6 mm) with the centre
fixed at a series of values, width and amplitudes free. Same data/errors as beam_profile_gauss_fit.py.
Also prints the model-free numbers: light ratios across the window and the log-slope/curvature.
"""
import os
from pathlib import Path
import h5py
import numpy as np
from scipy import optimize

ROOT = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10")
DARK = 0.0527


def load_run(day, first, last):
    x, t, v, f = [], [], [], []
    for n in range(first, last + 1):
        for p in sorted((day / f"{n:04d}").glob("*.h5")):
            with h5py.File(p, "r") as h5:
                x.append(float(h5["globals"].attrs["iris_x_mm"]))
                v.append(np.median(h5["data/traces/TOF_florescence"][()]["values"]) - DARK)
            t.append(os.path.getmtime(p)); f.append(n)
    x, t, v, f = map(np.array, (x, t, v, f))
    ref_f = sorted(set(f[x == 0]))
    rt = np.array([t[f == k].mean() for k in ref_f]); rv = np.array([v[f == k].mean() for k in ref_f])
    vc = v / (np.interp(t, rt, rv) / rv.mean())
    floor = rv.std() / rv.mean()
    pos = np.unique(x)
    mean = np.array([vc[x == p].mean() for p in pos])
    sem = np.array([vc[x == p].std(ddof=1) / np.sqrt((x == p).sum()) for p in pos])
    return pos, mean, np.sqrt(sem ** 2 + (floor * mean) ** 2 + 0.005 ** 2)


def gauss(x, a, x0, w):
    return a * np.exp(-2 * (x - x0) ** 2 / w ** 2)


sets = [load_run(ROOT / "16", 23, 35), load_run(ROOT / "17", 3, 13)]
sel = [np.abs(x) <= 6 for x, _, _ in sets]
for (x, y, _), s in zip(sets, sel):
    ly = np.log(y[s])
    c = np.polyfit(x[s], ly, 2)
    print(f"ln(light) quadratic over |x|<=6: slope at 0 {c[1]:+.3f}/mm, curvature {2*c[0]:+.4f}/mm^2 "
          f"-> Gaussian w {np.sqrt(-4/(2*c[0])) if c[0] < 0 else np.inf:.1f} mm, x0 {c[1]/(-2*c[0]) if c[0] < 0 else np.nan:+.1f} mm; "
          f"light(-6)/light(+6) = {y[x == -6][0]/y[x == 6][0]:.1f}")
print("\nfixed centre x0 -> best w, chi2 (10 dof at the free optimum, 11 here)")
for x0 in (-4, -6, -8, -10, -12, -14, -16, -20, -25, -30):
    def resid(p):
        return np.concatenate([(y[s] - gauss(x[s], p[i], x0, p[2])) / e[s] for i, ((x, y, e), s) in enumerate(zip(sets, sel))])
    r = optimize.least_squares(resid, [10, 10, 15], bounds=([0, 0, 2], [np.inf, np.inf, 80]))
    print(f"  x0 {x0:+4d} mm: w {r.x[2]:5.1f} mm, chi2 {np.sum(r.fun**2):6.1f}")
