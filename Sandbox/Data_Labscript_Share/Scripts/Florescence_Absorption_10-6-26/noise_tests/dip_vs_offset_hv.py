"""Absorption dip (%) of the upward cloud versus horizontal and vertical offset, on the same mm scale.

Usage: python dip_vs_offset_hv.py <helper dir> <day dir> <dark mV> <first folder> <last folder>
                                  <figure.png> [<dropout fraction>]
Example (2026-10-17, line centre): ... 2026/10/17 52.7 3 13 dip_vs_offset_hv_1017.png 0.6

There is no vertical stage scan: only the horizontal iris stage (global iris_x_mm) was moved, and
the vertical stage stayed put. The vertical profile comes from the transit through the sheet
instead: the cloud rises through the 1 mm probe sheet at v = 2.30 m/s, so the dip versus time at a
fixed iris position is the dip versus height in the cloud, z = v (t_c - t) (the top arrives first).

Left, horizontal: per shot the fixed-shape peak dip D (quadratic baseline x (1 - D Gaussian), the
run's common shape), PMT-normalised as in iris_atom_scan.py (D * mean(F/L~) / (F/L~)), averaged
per iris position, and fitted with a Gaussian averaged over the 2 mm iris (|x| <= 8 mm).

Right, vertical: the shots at iris 0 mm from the same run. Each trace is divided by a quadratic
baseline fitted only more than BASE_FAR_S from the PMT cloud time, scaled by the same PMT factor as
its D, and averaged in 1 mm bins (the sheet height; finer bins only oversample).
The profile is not Gaussian: it has a core plus a long leading tail (+z, atoms above the core that
arrive early), still ~3% of the peak at +20 mm. Checks (2026-10-19, scratch scripts):
- the tail is the same for quadratic or linear baselines over several far windows;
- the PMT, dead-time corrected (tau 60 ns), shows the same early tail out to -16 ms (z ~ +37 mm),
  while a no-atom shot (10/17 0000 #7, lock lost) is flat to +-0.2% of the peak. The PMT background
  is scattered probe light, so a probe-power hump would move it opposite to the PDA dip; atoms move
  both the same way. So the tail is atoms.
- Model comparison over z <= 23 mm (chi2 with per-bin sem, which ignores the bin-to-bin
  correlation, so compare between models only): Gaussian 1224/93, two-sided Gaussian 688/92,
  exponentially modified Gaussian 286/92, two Gaussians 202/90. The figure uses two Gaussians: a
  core plus a broad component above it.
The tail runs into the baseline window (the record starts 25 ms before the cloud), so the far tail
is partly absorbed into the baseline and the broad component's share is a lower bound.
"""
import sys
import warnings
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.special import erf

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DAY, DARK = Path(sys.argv[2]), float(sys.argv[3]) * 1e-3
FIRST, LAST, FIG = int(sys.argv[4]), int(sys.argv[5]), Path(sys.argv[6])
DROP_FRAC = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
OFFSET, WIDTH = 0.05e-3, 1.51e-3          # common dip shape of the 10/17 run (iris_atom_scan.py output)
D_IRIS_MM, FIT_MAX_MM = 2.0, 8.0
V_PROBE = 2.30                            # m/s through the probe (0130 crossing times)
BASE_FAR_S = 0.010                        # vertical: baseline fitted only at |t - t_c| > this
BIN_S = 1e-3 / V_PROBE                    # vertical bin = 1 mm of height (the sheet height)
Z_RANGE_MM = (-12.0, 30.0)                # vertical panel range; horizontal panel uses H_RANGE_MM
H_RANGE_MM = (-12.0, 12.0)
# Colors as in the other noise_tests figures: data blue, model orange, ink, muted grey, surface
BLUE, ORANGE, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"

# --- per shot: fixed-shape D, light, PMT integral, and the fractional-dip trace ----------------------
shots = []
for n in range(FIRST, LAST + 1):
    for p in sorted((DAY / f"{n:04d}").glob("*.h5")):
        shot = load_shot(p)
        with h5py.File(p, "r") as h5:
            x = float(h5["globals"].attrs["iris_x_mm"])
        pmt = cse.pmt_cloud_fit(shot, shot["mod_freq"] or 6125.0)
        if shot["pda_v"].max() >= 9.95 or not pmt["ok"]:
            continue
        # PDA at 10 kHz inside the jump-enabled window (2 ms guards), dark subtracted
        k = int(round(shot["fs_pda"] / 10e3))
        m = len(shot["pda_v"]) // k * k
        tb = shot["t_pda"][:m].reshape(-1, k).mean(axis=1)
        vb = shot["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
        win = (tb > shot["enable_start_s"] + 0.002) & (tb < shot["enable_stop_s"] - 0.002)
        tw, vw, tc = tb[win], vb[win], pmt["center_s"]

        def fixed(t, c0, c1, c2, d):
            u = t - tc
            return (c0 + c1 * u + c2 * u ** 2) * (1 - d * np.exp(-0.5 * ((t - tc - OFFSET) / WIDTH) ** 2))

        pf, _ = optimize.curve_fit(fixed, tw, vw, p0=[np.median(vw), 0, 0, 0.01])
        # Fractional dip versus time from the PMT centre, baseline from well outside the cloud
        u = tw - tc
        far = np.abs(u) > BASE_FAR_S
        frac = 1 - vw / np.polyval(np.polyfit(u[far], vw[far], 2), u)   # positive = absorption
        shots.append({"x": x, "D": pf[3], "light": pf[0], "F": pmt["amp_cps"] * pmt["sigma_s"],
                      "u": u, "frac": frac})
print(f"{len(shots)} shots from {FIRST:04d}-{LAST:04d}; window {1e3*shots[0]['u'].min():+.1f} to "
      f"{1e3*shots[0]['u'].max():+.1f} ms from the cloud")
if DROP_FRAC > 0:                          # lock-slip dropouts (see iris_atom_scan.py)
    med = np.median([s["F"] for s in shots])
    n0 = len(shots)
    shots = [s for s in shots if s["F"] >= DROP_FRAC * med]
    print(f"dropout cut at {DROP_FRAC:.2f} x median PMT: {n0 - len(shots)} dropped")

# PMT normalisation, keeping each position's mean (see iris_atom_scan.py)
xs = np.array([s["x"] for s in shots])
positions = np.unique(xs)
light_mean = {x: np.mean([s["light"] for s in shots if s["x"] == x]) for x in positions}
q = np.array([s["F"] / (s["light"] / light_mean[s["x"]]) for s in shots])
scale = q.mean() / q
dn = np.array([s["D"] for s in shots]) * scale

# --- horizontal profile ----------------------------------------------------------------------------
h_mean = np.array([dn[xs == x].mean() for x in positions])
h_sem = np.array([dn[xs == x].std(ddof=1) / np.sqrt((xs == x).sum()) for x in positions])


def profile(x, a, xc, sx):
    """Gaussian averaged over the iris width (analytic)."""
    h = D_IRIS_MM / 2
    return a * sx * np.sqrt(np.pi / 2) / D_IRIS_MM * (erf((x + h - xc) / (np.sqrt(2) * sx))
                                                     - erf((x - h - xc) / (np.sqrt(2) * sx)))


inside = np.abs(positions) <= FIT_MAX_MM
ph, ch = optimize.curve_fit(profile, positions[inside], h_mean[inside], p0=(h_mean.max(), 0, 3),
                            sigma=h_sem[inside], absolute_sigma=True)
chi2 = np.sum(((h_mean[inside] - profile(positions[inside], *ph)) / h_sem[inside]) ** 2)
eh = np.sqrt(np.diag(ch) * max(1, chi2 / (inside.sum() - 3)))
print(f"horizontal: peak {100*ph[0]:.2f} +- {100*eh[0]:.2f}%, centre {ph[1]:+.2f} +- {eh[1]:.2f} mm, "
      f"sigma_x {abs(ph[2]):.2f} +- {eh[2]:.2f} mm (chi2/dof {chi2:.1f}/{inside.sum()-3})")
print("  0 mm mean: " + f"{100*h_mean[positions == 0][0]:.2f} +- {100*h_sem[positions == 0][0]:.2f}%")

# --- vertical profile: the 0 mm shots, 1 mm bins ----------------------------------------------------
# Bin edges in time, chosen so the bins run over Z_RANGE_MM (z = -v u, so u runs the other way)
edges = np.arange(-Z_RANGE_MM[1] / 1e3 / V_PROBE, -Z_RANGE_MM[0] / 1e3 / V_PROBE + 1e-9, BIN_S)
centres = 0.5 * (edges[1:] + edges[:-1])
binned = []
for s, sc in zip(shots, scale):
    if s["x"] != 0:
        continue
    idx = np.digitize(s["u"], edges) - 1
    binned.append([sc * s["frac"][idx == i].mean() if np.any(idx == i) else np.nan
                   for i in range(len(centres))])
binned = np.array(binned)
v_mean = np.nanmean(binned, axis=0)
v_sem = np.nanstd(binned, axis=0, ddof=1) / np.sqrt(np.sum(np.isfinite(binned), axis=0))
z = 1e3 * V_PROBE * (-centres)                       # mm, + = upper part of the cloud


def gauss(zz, a, z0, sz):
    return a * np.exp(-0.5 * ((zz - z0) / sz) ** 2)


def core_plus_broad(zz, a1, z1, s1, a2, z2, s2):
    """Core Gaussian plus a broad Gaussian (the leading tail); see the docstring for the comparison."""
    return gauss(zz, a1, z1, s1) + gauss(zz, a2, z2, s2)


# Fit only where the baseline was not fitted (|u| < BASE_FAR_S), weighted by the per-bin sem
fit = np.abs(centres) < BASE_FAR_S
pv, _ = optimize.curve_fit(core_plus_broad, z[fit], v_mean[fit], sigma=v_sem[fit], maxfev=20000,
                           p0=(v_mean.max(), -0.5, 3, 0.1 * v_mean.max(), 5, 5))
area_core, area_broad = abs(pv[0] * pv[2]), abs(pv[3] * pv[5])
broad_frac = area_broad / (area_core + area_broad)
res = (v_mean[fit] - core_plus_broad(z[fit], *pv)) / v_sem[fit]
print(f"vertical: {len(binned)} shots at 0 mm, peak {100*np.nanmax(v_mean):.2f}% at "
      f"{z[np.nanargmax(v_mean)]:+.1f} mm; core sigma_z {abs(pv[2]):.2f} mm at {pv[1]:+.2f} mm; broad part "
      f"sigma {abs(pv[5]):.2f} mm at {pv[4]:+.2f} mm, {broad_frac:.0%} of the area "
      f"(chi2/dof {np.sum(res**2):.0f}/{fit.sum()-6}, sem ignores bin correlation)")
for zz in (-9, 9, 15, 20):
    i = np.argmin(np.abs(z - zz))
    print(f"  z {z[i]:+5.1f} mm: {100*v_mean[i]:.2f}% ({v_mean[i]/np.nanmax(v_mean):.1%} of peak)")

# --- figure: panels sized so 1 mm is the same width in both -------------------------------------------
span_h, span_v = np.diff(H_RANGE_MM)[0], np.diff(Z_RANGE_MM)[0]
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True, width_ratios=(span_h, span_v),
                              facecolor=SURFACE)
for a_ in (ax, ax2):
    a_.set_facecolor(SURFACE)
    a_.tick_params(colors=MUTED)
    for spine in a_.spines.values():
        spine.set_color(MUTED)
    a_.grid(True, color=MUTED, alpha=0.15)
    a_.axhline(0, color=MUTED, lw=1)
ax.set_xlim(*H_RANGE_MM)
ax2.set_xlim(*Z_RANGE_MM)
xh = np.linspace(*H_RANGE_MM, 400)
ax.plot(xs, 100 * dn, "o", ms=4, color=BLUE, alpha=0.3, mec="none")
ax.errorbar(positions, 100 * h_mean, yerr=100 * h_sem, fmt="o", ms=8, color=BLUE, mec=SURFACE, mew=1.5,
            capsize=3, lw=1.5, label="Mean ± sem")
ax.plot(xh, 100 * profile(xh, *ph), color=ORANGE, lw=2.2,
        label=f"Gaussian: σx {abs(ph[2]):.1f} mm")
ax.set_xlabel("Horizontal offset (mm, iris stage)", color=INK)
ax.set_ylabel("Peak dip (%)", color=INK)
ax.set_title("Horizontal: iris scan", color=INK)
ax.legend(frameon=False, fontsize=9, loc="upper left", labelcolor=INK)
xv = np.linspace(*Z_RANGE_MM, 600)
ax2.axvspan(1e3 * V_PROBE * BASE_FAR_S, Z_RANGE_MM[1], color=MUTED, alpha=0.08, lw=0)
ax2.text(1e3 * V_PROBE * BASE_FAR_S + 0.5, 0.55 * 100 * np.nanmax(v_mean), "baseline region\n(not fitted)",
         color=MUTED, fontsize=8, va="top")                 # below the legend
ax2.errorbar(z, 100 * v_mean, yerr=100 * v_sem, fmt="o", ms=5, color=BLUE, mec=SURFACE, mew=1,
             capsize=2, lw=1.2, label=f"Mean of {len(binned)} shots ± sem (1 mm bins)")
ax2.plot(xv, 100 * core_plus_broad(xv, *pv), color=ORANGE, lw=2.2, label="Fit: core + broad")
ax2.plot(xv, 100 * gauss(xv, *pv[:3]), color=ORANGE, lw=1.2, ls="--",
         label=f"Core: σz {abs(pv[2]):.1f} mm")
ax2.plot(xv, 100 * gauss(xv, *pv[3:]), color=ORANGE, lw=1.2, ls=":",
         label=f"Broad: σ {abs(pv[5]):.0f} mm, {broad_frac:.0%} of atoms")
ax2.set_xlabel("Vertical offset (mm, transit time × 2.30 m/s; + = up)", color=INK)
ax2.set_title("Vertical: dip vs time at iris 0 mm", color=INK)
ax2.legend(frameon=False, fontsize=9, loc="upper right", labelcolor=INK)
fig.tight_layout()
fig.savefig(FIG, dpi=130, facecolor=SURFACE)
print(f"figure: {FIG}")
