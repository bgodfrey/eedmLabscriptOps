"""Horizontal probe profile through a ~2 mm iris, drift-corrected, fitted as a clipped Gaussian.

Usage: python iris_profile_fit.py <day dir> <dark mV> <first folder> <last folder> <figure.png>
Example (2026-10-15): ... 2026/10/15 52.7 0017 0033 iris_profile_fit_0017-0033.png

Data: no-launch shots, one folder per stage position, position read from the global iris_x_mm.
Positions were visited in shuffled order, with returns to 0 mm (the reference) at the start, end
and twice in between.

Steps:
1. Level per shot = median PDA voltage over the 0.2 s record minus the dark level (run 0010).
2. Drift: the reference shots give the light level at 0 mm versus time. Linear interpolation
   between them gives a drift factor for every shot (reference level at that time / mean
   reference level); each shot is divided by it. Shots outside the reference span use the
   nearest reference.
3. Fit: the beam is a Gaussian I(u) = A exp(-2 (u - x0)^2 / w^2) (w = 1/e^2 radius), cut by an
   opening that passes lo < u < hi, seen through the iris (a top-hat of width D_IRIS_MM). The
   model level at stage offset x is the clipped beam averaged over x - D/2 .. x + D/2, written
   with erf so it varies smoothly with the edge positions (a sampled version gives the fitter
   flat steps and it cannot move the edges).
   Five free parameters (A, x0, w, lo, hi). Weights: each position's standard error combined with
   a 6% relative term (the rms spread of the reference visits: references are about 9 min apart
   while the power swings with a period of about 5 min, so the interpolation cannot follow it) and
   a 5 mV absolute term (dark-level drift). Parameter errors come from the fit covariance, scaled
   by the reduced chi-square when it exceeds 1.
"""
from scipy.special import erf
import os
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

DAY, DARK = Path(sys.argv[1]), float(sys.argv[2]) * 1e-3
FIRST, LAST, FIG = int(sys.argv[3]), int(sys.argv[4]), Path(sys.argv[5])
D_IRIS_MM = 2.0
# Colors as in the other noise_tests figures: data blue, model orange, ink, muted grey, surface
BLUE, ORANGE, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"

# --- 1. per-shot levels, positions and times ------------------------------------------------------
rows = []                                   # (position mm, time s, level V, folder)
for n in range(FIRST, LAST + 1):
    folder = DAY / f"{n:04d}"
    for p in sorted(folder.glob("*.h5")):
        with h5py.File(p, "r") as h5:
            x = float(h5["globals"].attrs["iris_x_mm"])
            level = np.median(h5["data/traces/TOF_florescence"][()]["values"]) - DARK
        rows.append((x, os.path.getmtime(p), level, folder.name))
x_all = np.array([r[0] for r in rows])
t_all = np.array([r[1] for r in rows])
t_all -= t_all.min()
v_all = np.array([r[2] for r in rows])

# --- 2. drift correction from the 0 mm references --------------------------------------------------
ref = x_all == 0
ref_folders = sorted({r[3] for r, is_ref in zip(rows, ref) if is_ref})
ref_t = np.array([t_all[ref & (np.array([r[3] for r in rows]) == f)].mean() for f in ref_folders])
ref_v = np.array([v_all[ref & (np.array([r[3] for r in rows]) == f)].mean() for f in ref_folders])
drift = np.interp(t_all, ref_t, ref_v) / ref_v.mean()
v_corr = v_all / drift
print("reference (0 mm) visits: " + ", ".join(
    f"{f} at {t/60:4.1f} min: {1e3*v:6.0f} mV" for f, t, v in zip(ref_folders, ref_t, ref_v)))
print(f"reference spread: {100*ref_v.std()/ref_v.mean():.1f}% rms, "
      f"range {100*(ref_v.max()-ref_v.min())/ref_v.mean():.1f}% of the mean")

positions = np.unique(x_all)
mean_raw = np.array([v_all[x_all == x].mean() for x in positions])
mean_c = np.array([v_corr[x_all == x].mean() for x in positions])
sem_c = np.array([v_corr[x_all == x].std(ddof=1) / np.sqrt((x_all == x).sum()) for x in positions])
print("\n  x mm   raw mV  corrected mV  (sem)  shots")
for x, r, c, e in zip(positions, mean_raw, mean_c, sem_c):
    print(f"{x:+6.1f}  {1e3*r:7.1f}  {1e3*c:9.1f}  ({1e3*e:5.1f})  {(x_all == x).sum():3d}")


# --- 3. clipped-Gaussian fit -----------------------------------------------------------------------
def model(x, a, x0, w, lo, hi):
    """Clipped Gaussian averaged over the iris width, evaluated at each stage offset x.

    The light inside the iris spans u from max(x - D/2, lo) to min(x + D/2, hi); the integral of
    exp(-2 (u - x0)^2 / w^2) over it is (w/2) sqrt(pi/2) [erf(sqrt2 (b - x0)/w) - erf(sqrt2 (a - x0)/w)].
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    lo_lim = np.maximum(x - D_IRIS_MM / 2, lo)
    hi_lim = np.minimum(x + D_IRIS_MM / 2, hi)
    hi_lim = np.maximum(hi_lim, lo_lim)                       # empty interval -> zero light
    s = np.sqrt(2) / w
    integral = (w / 2) * np.sqrt(np.pi / 2) * (erf(s * (hi_lim - x0)) - erf(s * (lo_lim - x0)))
    return a * integral / D_IRIS_MM


sigma = np.sqrt(sem_c ** 2 + (0.06 * mean_c) ** 2 + 0.005 ** 2)
p0 = (mean_c.max(), positions[np.argmax(mean_c)], 12.0, positions.min() + 1.5, positions.max() - 1.5)
p, cov = curve_fit(model, positions, mean_c, p0=p0, sigma=sigma, absolute_sigma=True, maxfev=20000)
chi2 = np.sum(((mean_c - model(positions, *p)) / sigma) ** 2)
dof = len(positions) - len(p)
err = np.sqrt(np.diag(cov) * max(1.0, chi2 / dof))
names = ("peak A (V)", "beam centre x0 (mm)", "1/e^2 radius w (mm)", "lower edge (mm)", "upper edge (mm)")
print(f"\nclipped-Gaussian fit (iris {D_IRIS_MM:.0f} mm top-hat), chi2/dof = {chi2:.1f}/{dof}:")
for name, val, e in zip(names, p, err):
    print(f"  {name:22s} {val:8.2f} +- {e:.2f}")
width = p[4] - p[3]
centre = (p[4] + p[3]) / 2
print(f"  opening width {width:.1f} mm, opening centre {centre:+.1f} mm; "
      f"beam centre is {p[1]-centre:+.1f} mm from the opening centre")
resid = (mean_c - model(positions, *p)) / mean_c.max()
print("  residuals (% of max): " + ", ".join(f"{x:+.0f}: {100*r:+.1f}" for x, r in zip(positions, resid)))

# --- figure: drift (top), profile and fit (bottom) -------------------------------------------------
fig, (ax0, ax) = plt.subplots(2, 1, figsize=(8, 7.5), height_ratios=(1, 2.6), facecolor=SURFACE)
for a_ in (ax0, ax):
    a_.set_facecolor(SURFACE)
    a_.tick_params(colors=MUTED)
    for spine in a_.spines.values():
        spine.set_color(MUTED)
    a_.grid(True, color=MUTED, alpha=0.15)
ax0.plot(t_all[ref] / 60, 1e3 * v_all[ref], "o", ms=5, color=BLUE, alpha=0.6, mec="none")
tt = np.linspace(t_all.min(), t_all.max(), 200)
ax0.plot(tt / 60, 1e3 * np.interp(tt, ref_t, ref_v), color=ORANGE, lw=2)
ax0.set_xlabel("Time (min)", color=INK)
ax0.set_ylabel("Level at 0 mm (mV)", color=INK)
ax0.set_title("Drift reference", color=INK)

xx = np.linspace(positions.min() - 2, positions.max() + 2, 600)
unclipped = p[0] * np.exp(-2 * (xx - p[1]) ** 2 / p[2] ** 2)
ax.axvspan(p[3], p[4], color=MUTED, alpha=0.07, lw=0)
ax.plot(xx, unclipped, ls="--", color=ORANGE, lw=1.5, alpha=0.7, label="Beam without clipping")
ax.plot(xx, model(xx, *p), color=ORANGE, lw=2.2,
        label=f"Fit: w {p[2]:.1f} mm, edges {p[3]:+.1f} / {p[4]:+.1f} mm")
ax.plot(x_all, v_corr, "o", ms=4, color=BLUE, alpha=0.35, mec="none")
ax.errorbar(positions, mean_c, yerr=sem_c, fmt="o", ms=8, color=BLUE, mec=SURFACE, mew=1.5,
            capsize=3, lw=1.5, label="Mean, drift-corrected")
ax.axhline(0, color=MUTED, lw=1)
ax.set_xlabel("Stage offset (mm; + = left, facing N)", color=INK)
ax.set_ylabel("PDA above dark (V, 70 dB)", color=INK)
ax.set_title(f"Iris profile, no launch ({FIRST:04d}–{LAST:04d})", color=INK)
ax.legend(frameon=False, fontsize=9, loc="upper left", labelcolor=INK)
fig.tight_layout()
fig.savefig(FIG, dpi=130, facecolor=SURFACE)
print(f"\nfigure: {FIG}")
