"""Horizontal absorption profile from the time-integrated dip (no dip template), and the atom number.

Usage: python integrated_dip_profile.py <helper dir> <day dir> <dark mV> <first folder> <last folder>
                                        <figure.png> [<dropout fraction>] [--show]
Example (2026-10-17, line centre): ... 2026/10/17 52.7 3 13 integrated_dip_1017.png 0.6
Example (2026-10-16, mod_offset 0): ... 2026/10/16 52.7 23 35 integrated_dip_1016.png

Why (2026-10-20): iris_atom_scan.py measures each shot with a fixed-shape dip D (one Gaussian in
time, width fixed to the run's median, 1.51 ms on 10/17). Free fits show the dip lasts longer away
from the cloud centre (1.41 ms at 0 mm, ~1.6 ms at +-4 mm, 2.2-3 ms at +-6-8 mm), so one fixed width
under-reads the absorption in the wings, and with it the profile width and the atom number. Here
each shot's fractional dip is summed over time instead, with no shape assumed.

Per shot:
1. PMT Gaussian fit -> cloud time t_c (as iris_atom_scan.py).
2. PDA averaged to 10 kHz inside the jump-enabled window (2 ms guards), dark subtracted.
3. Baseline: quadratic in u = t - t_c, fitted only where |u| > HALF_S (both sides of the cloud).
4. Fractional dip alpha(t) = 1 - V / baseline (positive = absorption). Integrated dip
   I = sum of alpha * dt over |u| <= HALF_S, in % x ms. Model-free apart from the baseline.
5. The fixed-shape D (iris_atom_scan.py pass 2, the run's shape from SHAPES) for comparison:
   the old analysis used D * sqrt(2 pi) * width as the integrated dip. The effective width
   I / (sqrt(2 pi) D) equals the template width only if the dip has the template's shape.
6. PMT normalisation as in iris_atom_scan.py: x mean(F / L~) / (F / L~), with F the PMT integral
   and L~ the shot's light over the mean light at its position. Keeps each position's mean.

Limits: the leading tail runs to about -16 ms (notebook 2026-10-19), past HALF_S, so the far tail
falls in the baseline region and I is a lower bound; the printout repeats the totals at several
HALF_S to show how much that matters. Positions outside the lit range (|x| > FIT_MAX_MM) are
plotted hollow and left out of the fit and the sum.

Atom number (as iris_atom_scan.py): optically thin, double pass, cloud crossing the sheet at v,
    N = v / (2 sigma_eff) * integral over x of I(x),   sigma_eff = sigma_0 / (1 + s).
The x integral is the sum over the measured positions times the 2 mm step (model-free, nothing
beyond the lit range) and, for comparison, the area of the iris-averaged Gaussian fitted to I(x).
Off resonance (10/16, mod_offset 0) sigma_eff is lower; multiply that N by the 0036 ratio, 3.6 +- 0.4.
"""
import sys
import warnings
from pathlib import Path
import h5py
import numpy as np
import matplotlib
SHOW = "--show" in sys.argv
if not SHOW:
    matplotlib.use("Agg")                      # file only, unless asked to open a window
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.special import erf

ARGS = [a for a in sys.argv[1:] if a != "--show"]
sys.path.insert(0, ARGS[0])
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DAY, DARK = Path(ARGS[1]), float(ARGS[2]) * 1e-3
FIRST, LAST, FIG = int(ARGS[3]), int(ARGS[4]), Path(ARGS[5])
DROP_FRAC = float(ARGS[6]) if len(ARGS) > 6 else 0.0
D_IRIS_MM = 2.0
FIT_MAX_MM = 8.0                               # lit range; beyond it the light sits on a hard edge
V_PROBE = 2.30                                 # m/s through the probe sheet
HALF_S = 0.010                                 # integrate |t - t_c| <= this; baseline outside it
HALF_S_CHECK = (0.008, 0.010, 0.012)           # repeated in the printout as a check
MIN_BASE_S = 0.003                             # need at least this much baseline on each side
S_PROBE = 0.2                                  # saturation parameter, as iris_atom_scan.py
SIGMA0_CM2 = {"isotropic": 1.405e-9, "cycling": 3.468e-9}   # Cs D2 F=4 -> F'=5 (Steck)
# Fixed dip shape per day (centre from the PMT t_c, width sigma), from iris_atom_scan.py pass 1
SHAPES = {"16": (-0.13e-3, 1.48e-3), "17": (0.05e-3, 1.51e-3)}
OFFSET, WIDTH = SHAPES[DAY.name]
BLUE, ORANGE, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"


def integrated_dip(u, v, half):
    """Quadratic baseline fitted at |u| > half, then the sum of the fractional dip over |u| <= half.
    Returns (I in s, baseline at the cloud time in V, fractional-dip trace, baseline trace)."""
    far = np.abs(u) > half
    c = np.polyfit(u[far], v[far], 2)
    base = np.polyval(c, u)
    alpha = 1 - v / base
    dt = np.median(np.diff(u))
    return np.sum(alpha[~far]) * dt, c[-1], alpha, base


# --- 1. per shot: integrated dip at each half-window, fixed-shape D, light, PMT integral ------------
shots, short = [], 0
for n in range(FIRST, LAST + 1):
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
        u = tw - tc
        # The widest half-window must still leave baseline on both sides of the cloud
        if (u.min() > -max(HALF_S_CHECK) - MIN_BASE_S) or (u.max() < max(HALF_S_CHECK) + MIN_BASE_S):
            short += 1
            continue
        ints = {h: integrated_dip(u, vw, h)[0] for h in HALF_S_CHECK}
        _, light, _, _ = integrated_dip(u, vw, HALF_S)

        # Fixed-shape D, exactly as iris_atom_scan.py pass 2 (quadratic baseline x (1 - D Gaussian))
        def fixed(t, c0, c1, c2, d):
            uu = t - tc
            return (c0 + c1 * uu + c2 * uu ** 2) * (1 - d * np.exp(-0.5 * ((t - tc - OFFSET) / WIDTH) ** 2))
        pf, _ = optimize.curve_fit(fixed, tw, vw, p0=[np.median(vw), 0, 0, 0.01])
        shots.append({"x": x, "I": ints, "D": pf[3], "light": light, "F": pmt["amp_cps"] * pmt["sigma_s"]})
print(f"{len(shots)} shots from {DAY.name}/{FIRST:04d}-{LAST:04d}"
      + (f" ({short} skipped: cloud too close to the window edge)" if short else ""))
if DROP_FRAC > 0:                              # lock-slip dropouts, as iris_atom_scan.py
    med = np.median([s["F"] for s in shots])
    n0 = len(shots)
    shots = [s for s in shots if s["F"] >= DROP_FRAC * med]
    print(f"dropout cut at {DROP_FRAC:.2f} x median PMT: {n0 - len(shots)} dropped")

# --- 2. PMT normalisation (keeps each position's mean) ---------------------------------------------
xs = np.array([s["x"] for s in shots])
positions = np.unique(xs)
light_mean = {x: np.mean([s["light"] for s in shots if s["x"] == x]) for x in positions}
q = np.array([s["F"] / (s["light"] / light_mean[s["x"]]) for s in shots])
scale = q.mean() / q
I_all = {h: np.array([s["I"][h] for s in shots]) * scale for h in HALF_S_CHECK}
D_all = np.array([s["D"] for s in shots]) * scale
I_shot = I_all[HALF_S]


def per_pos(vals):
    """Mean and standard error at each iris position."""
    mean = np.array([vals[xs == x].mean() for x in positions])
    sem = np.array([vals[xs == x].std(ddof=1) / np.sqrt((xs == x).sum()) for x in positions])
    return mean, sem


I_mean, I_sem = per_pos(I_shot)
D_mean, D_sem = per_pos(D_all)
# The old integrated-dip estimate: D x sqrt(2 pi) x template width
T_mean, T_sem = D_mean * np.sqrt(2 * np.pi) * WIDTH, D_sem * np.sqrt(2 * np.pi) * WIDTH
# Effective width: the Gaussian sigma that would give the measured I at the measured D
w_eff = I_mean / (np.sqrt(2 * np.pi) * D_mean)
w_err = np.abs(w_eff) * np.sqrt((I_sem / I_mean) ** 2 + (D_sem / D_mean) ** 2)   # ignores I-D correlation
inside = np.abs(positions) <= FIT_MAX_MM

print(f"\nIntegrated dip over |t - t_c| <= {1e3*HALF_S:.0f} ms (PMT-normalised), vs the fixed-shape estimate "
      f"D x sqrt(2pi) x {1e3*WIDTH:.2f} ms")
# Medians as a robustness check: the integrated dip picks up slow baseline errors, and one bad shot
# out of 8 can move a mean (10/17 -6 mm has a shot at -9.7 % ms)
I_med = np.array([np.median(I_shot[xs == x]) for x in positions])
print("  x mm  shots   I (% ms)  (sem)   I median   D x sqrt(2pi) w (% ms)   I / that   effective width (ms)")
for i, x in enumerate(positions):
    flag = "" if inside[i] else "   (outside lit range)"
    print(f"  {x:+5.1f}  {(xs == x).sum():4d}   {1e5*I_mean[i]:6.3f} ({1e5*I_sem[i]:.3f})   {1e5*I_med[i]:6.3f}     "
          f"{1e5*T_mean[i]:6.3f} ({1e5*T_sem[i]:.3f})          {I_mean[i]/T_mean[i]:5.2f}     "
          f"{1e3*w_eff[i]:5.2f} +- {1e3*w_err[i]:.2f}{flag}")
# Shot-to-shot noise of I against that of the old estimate, at the well-measured 0 mm visits
z0 = xs == 0
print(f"  shot sd at 0 mm: I {1e5*I_shot[z0].std(ddof=1):.2f} % ms, old {1e5*np.sqrt(2*np.pi)*WIDTH*D_all[z0].std(ddof=1):.2f} % ms")


# --- 3. profile fits: iris-averaged Gaussian, to I(x) and to the old estimate ---------------------
def profile(x, a, xc, sx):
    """Gaussian of peak a averaged over the iris width (analytic)."""
    x = np.asarray(x, dtype=float)
    h = D_IRIS_MM / 2
    return a * sx * np.sqrt(np.pi / 2) / D_IRIS_MM * (erf((x + h - xc) / (np.sqrt(2) * sx))
                                                     - erf((x - h - xc) / (np.sqrt(2) * sx)))


fits = {}
for label, mean, sem in (("integrated I", I_mean, I_sem), ("fixed-shape D", T_mean, T_sem)):
    xf, mf, ef = positions[inside], mean[inside], sem[inside]
    pp, cov = optimize.curve_fit(profile, xf, mf, p0=(mf.max(), xf[np.argmax(mf)], 3.0), sigma=ef,
                                 absolute_sigma=True)
    chi2 = np.sum(((mf - profile(xf, *pp)) / ef) ** 2)
    dof = len(xf) - 3
    err = np.sqrt(np.diag(cov) * max(1.0, chi2 / dof))
    fits[label] = (pp, err, chi2, dof)
    print(f"\nprofile ({label}): chi2/dof {chi2:.1f}/{dof}; centre {pp[1]:+.2f} +- {err[1]:.2f} mm, "
          f"sigma_x {abs(pp[2]):.2f} +- {err[2]:.2f} mm")
# Model-free rms width of I(x) over the lit points, minus the iris's own (d^2/12) variance
xi, Ii = positions[inside], I_mean[inside]
xbar = np.sum(xi * Ii) / np.sum(Ii)
rms = np.sqrt(np.sum((xi - xbar) ** 2 * Ii) / np.sum(Ii) - D_IRIS_MM ** 2 / 12)
print(f"model-free (lit points only, so truncated at +-{FIT_MAX_MM:.0f} mm): centroid {xbar:+.2f} mm, "
      f"rms width {rms:.2f} mm")

# --- 4. atom number ---------------------------------------------------------------------------------
step = np.diff(positions).min()


def atoms(area_s_mm):
    """N for an x-integrated dip in s*mm, both cross-sections."""
    return {name: V_PROBE * area_s_mm * 1e-3 / (2 * s0 * 1e-4 / (1 + S_PROBE)) for name, s0 in SIGMA0_CM2.items()}


print(f"\nAtom number (s = {S_PROBE}; isotropic / cycling):")
pI = fits["integrated I"][0]
area_rows = [("integrated I, sum of lit points", np.sum(I_mean[inside]) * step),
             ("integrated I, sum of per-position medians", np.sum(I_med[inside]) * step),
             ("integrated I, Gaussian area", pI[0] * abs(pI[2]) * np.sqrt(2 * np.pi)),
             ("old: fixed-shape D, sum of lit points", np.sum(T_mean[inside]) * step)]
for label, area in area_rows:
    nn = atoms(area)
    print(f"  {label:40s} {nn['isotropic']:.2e} / {nn['cycling']:.2e}")
print(f"  integrated / old (sums): {area_rows[0][1]/area_rows[-1][1]:.2f}")
print("\nCheck, integration half-window (I at 0 mm; sum over lit points; N isotropic):")
for h in HALF_S_CHECK:
    m_h, _ = per_pos(I_all[h])
    area = np.sum(m_h[inside]) * step
    print(f"  +-{1e3*h:4.0f} ms: I(0) {1e5*m_h[positions == 0][0]:.3f} % ms; N {atoms(area)['isotropic']:.2e}")

# --- figure: profile (top), effective width (bottom) --------------------------------------------------
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8, 7.5), height_ratios=(2, 1), sharex=True, facecolor=SURFACE)
for a_ in (ax, ax2):
    a_.set_facecolor(SURFACE)
    a_.tick_params(colors=MUTED)
    for spine in a_.spines.values():
        spine.set_color(MUTED)
    a_.grid(True, color=MUTED, alpha=0.15)
xx = np.linspace(positions.min() - 1, positions.max() + 1, 400)
in_shot = np.abs(xs) <= FIT_MAX_MM
ax.axhline(0, color=MUTED, lw=1)
ax.plot(xs[in_shot], 1e5 * I_shot[in_shot], "o", ms=4, color=BLUE, alpha=0.25, mec="none")
ax.errorbar(positions[inside], 1e5 * I_mean[inside], yerr=1e5 * I_sem[inside], fmt="o", ms=8, color=BLUE,
            mec=SURFACE, mew=1.5, capsize=3, lw=1.5, label="Integrated dip (mean ± sem)")
if (~inside).any():
    ax.errorbar(positions[~inside], 1e5 * I_mean[~inside], yerr=1e5 * I_sem[~inside], fmt="o", ms=8,
                mfc=SURFACE, color=MUTED, mew=1.5, capsize=3, lw=1.2, label="Outside lit range (not used)")
ax.errorbar(positions[inside], 1e5 * T_mean[inside], yerr=1e5 * T_sem[inside], fmt="s", ms=6, mfc=SURFACE,
            color=MUTED, mew=1.5, capsize=3, lw=1.2, label=f"Old: D × √(2π) × {1e3*WIDTH:.2f} ms")
pp, err, chi2, dof = fits["integrated I"]
ax.plot(xx, 1e5 * profile(xx, *pp), color=ORANGE, lw=2.2,
        label=f"Gaussian fit to integrated dip: σx {abs(pp[2]):.1f} mm (χ²/dof {chi2:.0f}/{dof})")
ax.set_ylabel("Integrated dip (% · ms)", color=INK)
ax.set_title(f"Upward cloud, time-integrated dip vs iris ({DAY.parent.name}/{DAY.name} "
             f"{FIRST:04d}–{LAST:04d})", color=INK)
ax.legend(frameon=False, fontsize=9, loc="upper left", labelcolor=INK)
ax2.axhline(1e3 * WIDTH, color=MUTED, lw=1, ls="--")
ax2.text(positions.max() + 0.6, 1e3 * WIDTH, "template", color=MUTED, fontsize=8, va="bottom", ha="right")
ax2.errorbar(positions[inside], 1e3 * w_eff[inside], yerr=1e3 * w_err[inside], fmt="o-", ms=7, color=BLUE,
             mec=SURFACE, mew=1.2, capsize=3, lw=1.2)
ax2.set_ylabel("Effective width (ms)", color=INK)
ax2.set_xlabel("Stage offset (mm; + = left, facing N)", color=INK)
ax2.set_ylim(bottom=0)
fig.tight_layout()
fig.savefig(FIG, dpi=130, facecolor=SURFACE)
print(f"\nfigure: {FIG}")
if SHOW:
    plt.show()
