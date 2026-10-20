"""Quick-and-dirty check of the probe light profile: a plain Gaussian fitted inside the opening,
extrapolated outward, to see where the light falls off.

Usage: python beam_profile_gauss_fit.py <figure.png>
(Runs as in beam_profile_softedge_fit.py: 2026-10-15 0017-0033, 10/16 0023-0035, 10/17 0003-0013.)

Data: identical to beam_profile_softedge_fit.py (per-shot PDA median minus the 52.7 mV dark level,
drift-corrected by the 0 mm visits; error = sem + drift floor + 5 mV, combined in quadrature).
Model: I(x) = A exp(-2 (x - x0)^2 / w^2), w = 1/e^2 radius. No edges and no iris averaging (the
2 mm iris adds only ~0.6 mm rms to a w ~ 15-20 mm beam).
Fit: weighted least squares on the points with |x| <= FIT_MM only (well inside the cut-offs seen
at about -7.5 and +9 mm), so the edges do not drag the Gaussian. 10/16 and 10/17 share x0 and w
(same alignment) with their own amplitudes. Several starting centres; the best is kept. Errors from
the covariance, scaled by sqrt(chi2/dof) when that exceeds 1. Also refitted with |x| <= 4 mm as a
check on how much the window matters.
Bottom panels: data / Gaussian. Inside the window it should sit near 1; where it drops toward 0 is
where something cuts the beam, so this ratio is roughly the opening's transmission. Far outside
the window it leans on the extrapolation, which is uncertain because x0 and w trade off (only one
flank of the beam is seen on each date).

What to trust (2026-10-19): where the data fall below the curve (the cut-offs near -7.5 and +9 mm),
and the lopsidedness inside the window. The fitted centre and width are NOT physical (a 30-40 mm
beam centred 8-16 mm from a 16 mm opening, moving ~25 mm after an input nudge): the light inside
the opening is not a clean Gaussian beam, so the Gaussian is only a smooth reference curve here.
They are printed for the record but kept off the figure, and the curve is drawn only to +-10 mm.
"""
import os
import sys
import warnings
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import optimize

warnings.filterwarnings("ignore")
FIG = Path(sys.argv[1])
ROOT = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10")
DARK = 0.0527                                   # V at 70 dB (10/15 0010, probe blocked)
RUNS = [("10/15", ROOT / "15", 17, 33), ("10/16", ROOT / "16", 23, 35), ("10/17", ROOT / "17", 3, 13)]
FIT_MM = 6.0                                    # fit window |x| <= FIT_MM (the figure)
CHECK_MM = 4.0                                  # narrower window, printed as a check
BLUE, ORANGE, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"


def load_run(day, first, last):
    """Drift-corrected mean light per position and its error (same as beam_profile_softedge_fit.py)."""
    x, t, v, f = [], [], [], []
    for n in range(first, last + 1):
        for p in sorted((day / f"{n:04d}").glob("*.h5")):
            with h5py.File(p, "r") as h5:
                x.append(float(h5["globals"].attrs["iris_x_mm"]))
                v.append(np.median(h5["data/traces/TOF_florescence"][()]["values"]) - DARK)
            t.append(os.path.getmtime(p))
            f.append(n)
    x, t, v, f = map(np.array, (x, t, v, f))
    ref_f = sorted(set(f[x == 0]))
    rt = np.array([t[f == k].mean() for k in ref_f])
    rv = np.array([v[f == k].mean() for k in ref_f])
    vc = v / (np.interp(t, rt, rv) / rv.mean())
    floor = rv.std() / rv.mean()
    pos = np.unique(x)
    mean = np.array([vc[x == p].mean() for p in pos])
    sem = np.array([vc[x == p].std(ddof=1) / np.sqrt((x == p).sum()) for p in pos])
    return pos, mean, np.sqrt(sem ** 2 + (floor * mean) ** 2 + 0.005 ** 2)


def gauss(x, a, x0, w):
    return a * np.exp(-2 * (x - x0) ** 2 / w ** 2)


def fit_gauss(sets, window, x0_starts):
    """Joint weighted fit; sets = [(x, y, e), ...] share x0 and w, each has its own amplitude.
    Parameters: [a_1 .. a_n, x0, w]."""
    n = len(sets)
    sel = [np.abs(x) <= window for x, _, _ in sets]

    def resid(p):
        return np.concatenate([(y[s] - gauss(x[s], p[i], p[n], p[n + 1])) / e[s]
                               for i, ((x, y, e), s) in enumerate(zip(sets, sel))])
    best = None
    for x0 in x0_starts:
        p0 = [y.max() for _, y, _ in sets] + [x0, 15.0]
        r = optimize.least_squares(resid, p0, bounds=([0] * n + [-60, 2], [np.inf] * n + [60, 80]),
                                   x_scale="jac", max_nfev=20000)
        if best is None or r.cost < best.cost:
            best = r
    chi2 = np.sum(best.fun ** 2)
    dof = sum(s.sum() for s in sel) - (n + 2)
    cov = np.linalg.pinv(best.jac.T @ best.jac)
    err = np.sqrt(np.diag(cov) * max(1.0, chi2 / max(dof, 1)))
    return best.x, err, chi2, dof


data = {name: load_run(day, a, b) for name, day, a, b in RUNS}
groups = [("Before realignment (10/15)", ["10/15"], [BLUE], [0, 5, 10, 20]),
          ("After realignment (10/16, 10/17)", ["10/16", "10/17"], [BLUE, ORANGE], [0, -5, -10, -20])]
fits = []
for title, names, _, starts in groups:
    sets = [data[nm] for nm in names]
    for window in (FIT_MM, CHECK_MM):
        p, e, chi2, dof = fit_gauss(sets, window, starts)
        k = len(names)
        print(f"{title}, |x| <= {window:.0f} mm: chi2/dof {chi2:.1f}/{dof}; (not physical:) centre x0 {p[k]:+.1f} +- {e[k]:.1f} mm, "
              f"w (1/e^2 radius) {p[k+1]:.1f} +- {e[k+1]:.1f} mm, i.e. sigma {p[k+1]/2:.1f} mm, "
              f"FWHM {p[k+1]*1.177:.1f} mm; amplitudes " + ", ".join(f"{a:.2f} V" for a in p[:k]))
        if window == FIT_MM:
            fits.append((p, e, chi2, dof))
    # Where the data fall below the window's Gaussian: ratio at every position
    p = fits[-1][0]
    for i, nm in enumerate(names):
        x, y, er = data[nm]
        g = gauss(x, p[i], p[len(names)], p[len(names) + 1])
        print(f"   {nm} data / Gaussian: " + ", ".join(f"{xx:+.0f}: {yy/gg:.2f}" for xx, yy, gg in zip(x, y, g)))

# --- figure: top = volts with the Gaussian, bottom = data / Gaussian ------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 7.2), height_ratios=(1.6, 1), sharex=True, facecolor=SURFACE)
xx = np.linspace(-13, 13, 600)
for col, ((title, names, colors, _), (p, e, chi2, dof)) in enumerate(zip(groups, fits)):
    top, bot = axes[0, col], axes[1, col]
    k = len(names)
    x0, w = p[k], p[k + 1]
    for a_ in (top, bot):
        a_.set_facecolor(SURFACE)
        a_.tick_params(colors=MUTED)
        for spine in a_.spines.values():
            spine.set_color(MUTED)
        a_.grid(True, color=MUTED, alpha=0.15)
        a_.axvspan(-FIT_MM, FIT_MM, color=MUTED, alpha=0.07, lw=0)   # fit window
        a_.axhline(0, color=MUTED, lw=1)
    # Legend entry for the shared fit; no centre or width, which are not physical (see docstring)
    top.plot([], [], " ", label=f"Gaussian fitted to |x| ≤ {FIT_MM:.0f} mm (shaded), χ²/dof {chi2:.0f}/{dof}")
    y_top = 1.25 * max((y + er).max() for y, er in ((data[nm][1], data[nm][2]) for nm in names))
    near = np.abs(xx) <= 10                     # extend the curve only a little past the window
    for i, (nm, c) in enumerate(zip(names, colors)):
        x, y, er = data[nm]
        inw = np.abs(x) <= FIT_MM
        top.plot(xx[near], gauss(xx[near], p[i], x0, w), color=c, lw=1.5, ls="--")   # extrapolated
        top.plot(xx[np.abs(xx) <= FIT_MM], gauss(xx[np.abs(xx) <= FIT_MM], p[i], x0, w), color=c, lw=2.2)
        top.errorbar(x[inw], y[inw], yerr=er[inw], fmt="o", ms=7, color=c, mec=SURFACE, mew=1.5,
                     capsize=3, lw=1.5, label=f"{nm} (fitted)")
        top.errorbar(x[~inw], y[~inw], yerr=er[~inw], fmt="o", ms=7, mfc=SURFACE, color=c, mew=1.5,
                     capsize=3, lw=1.5, label=f"{nm} (not fitted)")
        g = gauss(x, p[i], x0, w)
        bot.errorbar(x, y / g, yerr=er / g, fmt="o-", ms=6, color=c, mec=SURFACE, mew=1.2, capsize=3, lw=1.2,
                     label=nm)
    bot.axhline(1, color=MUTED, lw=1, ls=":")
    top.set_title(title, color=INK)
    top.set_ylabel("Probe light (V)", color=INK)
    top.set_ylim(-0.05 * y_top, y_top)          # set by the data; the extrapolation may run off the top
    # Legend on the dim side of each panel, clear of the steep curve
    top.legend(frameon=False, fontsize=8, loc="upper left" if col == 0 else "upper right", labelcolor=INK)
    bot.set_ylim(-0.1, 1.4)
    bot.set_ylabel("Data / Gaussian", color=INK)
    bot.set_xlabel("Stage offset (mm; + = left, facing N)", color=INK)
    bot.set_xlim(-13, 13)
fig.tight_layout()
fig.savefig(FIG, dpi=130, facecolor=SURFACE)
print(f"figure: {FIG}")
