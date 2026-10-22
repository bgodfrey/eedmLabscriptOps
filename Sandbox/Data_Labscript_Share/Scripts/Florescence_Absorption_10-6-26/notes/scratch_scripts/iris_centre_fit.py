"""Find the probe-beam centre from no-launch iris shots (2026-10-15 0012-0016).

Usage: python iris_centre_fit.py <day dir> <dark mV> <folder>:<offset mm> [...]

For each folder: globals, PDA level per file (median over the record; also the last 100 ms, after
any shutter transient) minus the dark level. Then a Gaussian fit I(x) = A exp(-2 (x - x0)^2 / w^2)
to the per-folder means, with w the 1/e^2 radius. With five points and three parameters the fit
has two degrees of freedom, so the quoted errors are rough.
"""
import sys
from pathlib import Path
import h5py
import numpy as np
from scipy.optimize import curve_fit

DAY, DARK = Path(sys.argv[1]), float(sys.argv[2]) * 1e-3
G = ("MRR_TRIG_BOOL", "pda_gain_db", "iris_x_mm", "mod_amp", "t_PMT_delay")
xs, levels, spreads = [], [], []
for arg in sys.argv[3:]:
    folder, offset = arg.split(":")
    vals = []
    for p in sorted((DAY / folder).glob("*.h5")):
        with h5py.File(p, "r") as h5:
            g = {k: h5["globals"].attrs.get(k) for k in G}
            v = h5["data/traces/TOF_florescence"][()]["values"]
        vals.append((np.median(v) - DARK, np.median(v[len(v) // 2:]) - DARK, v.std()))
        print(f"{folder} {p.name[-8:]}: level {1e3*vals[-1][0]:8.1f} mV (last 100 ms {1e3*vals[-1][1]:8.1f}), "
              f"std {1e3*vals[-1][2]:5.2f} mV  " + ", ".join(f"{k}={g[k]}" for k in G))
    a = np.array(vals)
    xs.append(float(offset))
    levels.append(a[:, 0].mean())
    spreads.append(a[:, 0].std() if len(a) > 1 else np.nan)

xs, levels = np.array(xs), np.array(levels)
print("\noffset mm  level above dark (mV)")
for x, l, s in zip(xs, levels, spreads):
    print(f"{x:+8.1f}  {1e3*l:9.1f}" + (f"  (+-{1e3*s:.1f} across shots)" if np.isfinite(s) else ""))


def gauss(x, a, x0, w):
    return a * np.exp(-2 * (x - x0) ** 2 / w ** 2)


p, cov = curve_fit(gauss, xs, levels, p0=(levels.max(), xs[np.argmax(levels)], 9.0))
err = np.sqrt(np.diag(cov))
print(f"\nGaussian fit: peak {1e3*p[0]:.0f} mV, centre {p[1]:+.2f} +- {err[1]:.2f} mm (stage offset), "
      f"1/e^2 radius w = {abs(p[2]):.2f} +- {err[2]:.2f} mm (1/e^2 diameter {2*abs(p[2]):.1f} mm)")
resid = levels - gauss(xs, *p)
print("residuals (% of peak): " + ", ".join(f"{x:+.0f}: {100*r/p[0]:+.1f}" for x, r in zip(xs, resid)))
