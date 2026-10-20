"""Horizontal absorption profile of the upward cloud from an iris scan (launched shots, DC dip).

Usage: python iris_atom_scan.py <helper dir> <day dir> <dark mV> <first folder> <last folder> <figure.png>
                                [<dropout fraction>]
Example (2026-10-16): ... 2026/10/16 52.7 23 35 iris_atom_scan_0023-0035.png
Example (2026-10-17, line centre): ... 2026/10/17 52.7 3 13 iris_atom_scan_1017_0003-0013.png 0.6

Dropouts (optional last argument, default 0 = keep every shot): shots whose PMT integral is below
that fraction of the run's median are dropped before any fitting. The PMT sees the whole cloud
whatever the iris position, so a low PMT marks a shot where the laser slipped off resonance or
fewer atoms arrived (2026-10-16/17: lock slips). Use it at line centre, where normal shots vary
about 10%; at an off-resonance offset the normal PMT scatter is 20-30% and a cut would bias.

Data: exp12 launched shots, modulation off (mod_amp = 0), PDA100A2 at 70 dB, a ~2 mm iris on the
horizontal stage in front of the PDA. One folder per stage position, read from the global iris_x_mm.
Positions were visited in shuffled order with returns to 0 mm.

Per shot (same dip model as absorption_dip_atoms.py, the 0130 analysis):
1. PMT Gaussian fit -> the time the cloud is centred in the sheet (t_c), independent of the PDA.
2. PDA block-averaged to 10 kHz inside the jump-enabled window (2 ms guard at each edge), dark level
   (blocked-probe run) subtracted, so the fit sees only probe light.
3. Fit  V(t) = quadratic baseline * (1 - D exp(-(t - t0)^2 / (2 s_t^2))).  D is the peak fractional
   dip of the double-passed light through the iris. Because D is a ratio within one shot, slow probe
   power drift (the PM-fiber polarization swing) cancels.
   - Pass 1: centre and width free, to learn the dip's shape.
   - Pass 2: the shape fixed to the median of the clean pass-1 fits (centre relative to the PMT t_c,
     and width), so only the baseline and D are free. A common shape keeps D comparable between
     bright and dim positions, where a free width would chase noise.
4. Time-integrated fractional dip I = D sqrt(2 pi) s_t.
5. Normalisation. The MOT fluorescence at coils-off (ai2) is printed as an atom-number monitor
   (2026-10-16: steady to 2.5% and uncorrelated with D, so it is not used). The PMT fluorescence
   integral F (amplitude x transit width) sees the whole cloud in the sheet and shares the
   shot-to-shot changes of D (launched atoms reaching the probe, probe detuning; correlation
   +0.7 to +0.9). F also scales with probe power, so it is divided by L~, the shot's light level
   relative to the mean at its position. D_pmt = D * mean(F / L~) / (F / L~) removes those common
   changes and keeps each position's mean.

Profile: the per-position means are fitted with a Gaussian in x averaged over the iris width
(analytic erf form): D(x) = A * mean over [x - d/2, x + d/2] of exp(-(u - xc)^2 / (2 sx^2)).
Only positions inside the lit window (|x| <= FIT_MAX_MM) enter the fit and the sum; beyond it the
light is a few percent of the peak and sits on a hard edge, so a tiny pointing change mimics a dip.

Atom number (optically thin, cloud crossing the sheet vertically at speed v, double pass): at a
horizontal position x the light crosses the column density along the beam; integrating the
fractional dip over the transit integrates over the vertical too, so
    integral(alpha dt) = 2 sigma_eff lambda(x) / v,   lambda(x) = atoms per unit horizontal length,
    N = v / (2 sigma_eff) * integral over x of I(x).
The sheet height and the iris size drop out. sigma_eff = sigma_0 / (1 + s), and any detuning from
line centre lowers sigma_eff further, so N here is a lower bound in that respect.
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
from fluorescence_absorption_lockin import load_shot, load_mot_trace, fit_mot_loading
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DAY, DARK = Path(sys.argv[2]), float(sys.argv[3]) * 1e-3
FIRST, LAST, FIG = int(sys.argv[4]), int(sys.argv[5]), Path(sys.argv[6])
DROP_FRAC = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
D_IRIS_MM = 2.0
FIT_MAX_MM = 8.0                     # lit window on 2026-10-15/16: edges near -8 to -9.5 and +9 mm
V_PROBE = 2.30                       # m/s through the probe, from the up/down crossing times (0130)
S_PROBE = 0.2                        # saturation parameter assumed for the sigma_eff correction
SIGMA0_CM2 = {"isotropic": 1.405e-9, "cycling": 3.468e-9}   # Cs D2 F=4 -> F'=5 (Steck)
# Colors as in the other noise_tests figures: data blue, model orange, ink, muted grey, surface
BLUE, ORANGE, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"


def dip_model(t, tc, c0, c1, c2, depth, centre, width):
    """Quadratic baseline (in t - tc, for good conditioning) times a Gaussian dip."""
    u = t - tc
    return (c0 + c1 * u + c2 * u ** 2) * (1 - depth * np.exp(-0.5 * ((t - centre) / width) ** 2))


# --- 1. load every shot: PMT timing, the PDA window, the MOT signal --------------------------------
shots = []
for n in range(FIRST, LAST + 1):
    folder = DAY / f"{n:04d}"
    for p in sorted(folder.glob("*.h5"), key=lambda q: int(q.stem.split("_")[-1])):
        shot = load_shot(p)
        with h5py.File(p, "r") as h5:
            x = float(h5["globals"].attrs["iris_x_mm"])
        if shot["pda_v"].max() >= 9.95:
            print(f"{folder.name}/{p.name}: PDA clipped, skipped")
            continue
        pmt = cse.pmt_cloud_fit(shot, shot["mod_freq"] or 6125.0)   # f only sets the smoothing
        if not pmt["ok"]:
            print(f"{folder.name}/{p.name}: no PMT fit, skipped")
            continue
        # 10 kHz blocks inside the enabled window, dark subtracted
        k = int(round(shot["fs_pda"] / 10e3))
        m = len(shot["pda_v"]) // k * k
        tb = shot["t_pda"][:m].reshape(-1, k).mean(axis=1)
        vb = shot["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
        win = (tb > shot["enable_start_s"] + 0.002) & (tb < shot["enable_stop_s"] - 0.002)
        try:
            mot_v = fit_mot_loading(load_mot_trace(p))["atoms_v_at_coils_off"]
        except Exception:                                  # missing trace or failed loading fit
            mot_v = np.nan
        shots.append({"folder": folder.name, "x": x, "time": p.stat().st_mtime, "tw": tb[win],
                      "vw": vb[win], "tc": pmt["center_s"], "pmt_sigma": pmt["sigma_s"],
                      "pmt_int": pmt["amp_cps"] * pmt["sigma_s"], "mot_v": mot_v})
print(f"{len(shots)} shots loaded from {FIRST:04d}-{LAST:04d}")
# Optional dropout cut on the PMT integral (see the docstring)
if DROP_FRAC > 0:
    pmt_med = np.median([s["pmt_int"] for s in shots])
    dropped = [s for s in shots if s["pmt_int"] < DROP_FRAC * pmt_med]
    shots = [s for s in shots if s["pmt_int"] >= DROP_FRAC * pmt_med]
    print(f"dropout cut at {DROP_FRAC:.2f} x median PMT ({pmt_med:.0f} counts): {len(dropped)} dropped"
          + "".join(f"\n   {s['folder']} x {s['x']:+.0f} mm: PMT {s['pmt_int']:.0f}, MOT {1e3*s['mot_v']:.1f} mV"
                    for s in dropped))

# --- 2. pass 1: free centre and width, to learn the dip shape -------------------------------------
for s in shots:
    tw, vw, tc = s["tw"], s["vw"], s["tc"]
    p0 = [np.median(vw), 0, 0, 0.02, tc, 1.3e-3]
    lo = [-np.inf, -np.inf, -np.inf, -0.5, tc - 0.005, 0.3e-3]
    hi = [np.inf, np.inf, np.inf, 0.5, tc + 0.005, 5e-3]
    try:
        pf, cov = optimize.curve_fit(lambda t, *a: dip_model(t, tc, *a), tw, vw, p0=p0, bounds=(lo, hi))
        s["free"] = (pf[3], np.sqrt(cov[3, 3]), pf[4] - tc, pf[5])
    except (RuntimeError, ValueError):
        s["free"] = (np.nan, np.nan, np.nan, np.nan)
# Clean shots: dip resolved at 8 sigma and not pinned at a bound
clean = [s for s in shots if s["free"][0] > 8 * s["free"][1] and 0.35e-3 < s["free"][3] < 4.9e-3]
OFFSET = np.median([s["free"][2] for s in clean])
WIDTH = np.median([s["free"][3] for s in clean])
print(f"dip shape from {len(clean)} clean shots: centre {1e3*OFFSET:+.2f} ms from the PMT centre, "
      f"width sigma {1e3*WIDTH:.2f} ms (spread {1e3*np.std([s['free'][3] for s in clean]):.2f} ms); "
      f"PMT width {1e3*np.median([s['pmt_sigma'] for s in shots]):.2f} ms")

# --- 3. pass 2: shape fixed, fit baseline and depth -------------------------------------------------
for s in shots:
    tw, vw, tc = s["tw"], s["vw"], s["tc"]
    centre = tc + OFFSET
    pf, cov = optimize.curve_fit(lambda t, c0, c1, c2, d: dip_model(t, tc, c0, c1, c2, d, centre, WIDTH),
                                 tw, vw, p0=[np.median(vw), 0, 0, 0.01])
    s["depth"], s["depth_err"], s["light"] = pf[3], np.sqrt(cov[3, 3]), pf[0]
xs = np.array([s["x"] for s in shots])
positions = np.unique(xs)
# PMT normalisation: F / L~ per shot, L~ = light relative to the mean light at that position
light_mean = {x: np.mean([s["light"] for s in shots if s["x"] == x]) for x in positions}
for s in shots:
    s["q"] = s["pmt_int"] / (s["light"] / light_mean[s["x"]])
q_mean = np.mean([s["q"] for s in shots])
for s in shots:
    s["depth_norm"] = s["depth"] * q_mean / s["q"]

# --- 4. per-position summary --------------------------------------------------------------------------


def per_pos(key):
    vals = [np.array([s[key] for s in shots if s["x"] == x]) for x in positions]
    vals = [v[np.isfinite(v)] for v in vals]
    return (np.array([v.mean() for v in vals]), np.array([v.std(ddof=1) / np.sqrt(len(v)) for v in vals]),
            np.array([v.std(ddof=1) for v in vals]), np.array([len(v) for v in vals]))


d_mean, d_sem, d_sd, n_pos = per_pos("depth")
dn_mean, dn_sem, dn_sd, _ = per_pos("depth_norm")
l_mean, l_sem, _, _ = per_pos("light")
fit_err = np.array([np.median([s["depth_err"] for s in shots if s["x"] == x]) for x in positions])
print(f"\n  x mm  shots  light V   dip %  (sem)  shot sd  fit err   PMT-normed %  (sem)  shot sd")
for i, x in enumerate(positions):
    print(f"{x:+6.0f}  {n_pos[i]:4d}  {l_mean[i]:7.3f}  {100*d_mean[i]:6.3f} ({100*d_sem[i]:.3f})  "
          f"{100*d_sd[i]:6.3f}  {100*fit_err[i]:6.3f}    {100*dn_mean[i]:6.3f} ({100*dn_sem[i]:.3f})  {100*dn_sd[i]:6.3f}")

# Reference visits (0 mm) in time order: atom-number drift and MOT monitor
print("\n0 mm visits:")
for f in sorted({s["folder"] for s in shots if s["x"] == 0}):
    grp = [s for s in shots if s["folder"] == f]
    print(f"  {f}: dip {100*np.mean([s['depth'] for s in grp]):.3f}%, MOT {1e3*np.nanmean([s['mot_v'] for s in grp]):.1f} mV, "
          f"light {np.mean([s['light'] for s in grp]):.3f} V")
mots = np.array([s["mot_v"] for s in shots])
pmts = np.array([s["pmt_int"] for s in shots])
print(f"MOT signal: mean {1e3*np.nanmean(mots):.1f} mV, shot-to-shot {100*np.nanstd(mots)/np.nanmean(mots):.1f}% rms; "
      f"PMT integral: {100*pmts.std()/pmts.mean():.0f}% rms, corr with MOT {np.corrcoef(pmts, mots)[0,1]:+.2f}")
bright = [i for i, x in enumerate(positions) if d_mean[i] > 0.5 * d_mean.max()]
print(f"relative shot scatter at positions above half the peak: raw {100*np.median(d_sd[bright]/d_mean[bright]):.1f}%, "
      f"PMT-normalised {100*np.median(dn_sd[bright]/dn_mean[bright]):.1f}% (median over positions)")
tails = [i for i, x in enumerate(positions) if 6 <= abs(x) <= FIT_MAX_MM]
print(f"absolute shot scatter where the dip is small (|x| 6-{FIT_MAX_MM:.0f} mm): "
      f"{100*np.median(d_sd[tails]):.2f}% of the light; median fit error {100*np.median(fit_err[tails]):.2f}%")


# --- 5. Gaussian profile fit, averaged over the iris width -------------------------------------------
def profile(x, a, xc, sx):
    x = np.asarray(x, dtype=float)
    h = D_IRIS_MM / 2
    return a * sx * np.sqrt(np.pi / 2) / D_IRIS_MM * (erf((x + h - xc) / (np.sqrt(2) * sx))
                                                     - erf((x - h - xc) / (np.sqrt(2) * sx)))


inside = np.abs(positions) <= FIT_MAX_MM
results = {}
for label, mean, sem in (("raw", d_mean, d_sem), ("PMT-normalised", dn_mean, dn_sem)):
    xf, mf, ef = positions[inside], mean[inside], sem[inside]
    p0 = (mf.max(), xf[np.argmax(mf)], 3.0)
    p, cov = optimize.curve_fit(profile, xf, mf, p0=p0, sigma=ef, absolute_sigma=True)
    chi2 = np.sum(((mf - profile(xf, *p)) / ef) ** 2)
    dof = len(xf) - 3
    err = np.sqrt(np.diag(cov) * max(1.0, chi2 / dof))
    results[label] = (p, err, chi2, dof)
    print(f"\nprofile fit ({label}), chi2/dof {chi2:.1f}/{dof}: peak {100*p[0]:.2f} +- {100*err[0]:.2f}%, "
          f"centre {p[1]:+.2f} +- {err[1]:.2f} mm, sigma_x {abs(p[2]):.2f} +- {err[2]:.2f} mm")

# --- 6. atom number -------------------------------------------------------------------------------------
p, err, _, _ = results["PMT-normalised"]
area_fit = p[0] * np.sqrt(2 * np.pi) * WIDTH * abs(p[2]) * np.sqrt(2 * np.pi) * 1e-3   # m*s, whole Gaussian
step = np.diff(positions).min()
area_sum = np.sum(dn_mean[inside] * np.sqrt(2 * np.pi) * WIDTH) * step * 1e-3           # m*s, lit points
print(f"\nintegral of I(x) dx: {1e6*1e3*area_fit:.1f} um*mm (Gaussian), {1e6*1e3*area_sum:.1f} um*mm (sum of points)")
for name, s0 in SIGMA0_CM2.items():
    sig = s0 * 1e-4 / (1 + S_PROBE)
    print(f"  N ({name}, s = {S_PROBE}): {V_PROBE*area_fit/(2*sig):.2e} (Gaussian), {V_PROBE*area_sum/(2*sig):.2e} (sum)")

# --- figure: absorption profile (top), probe light at the cloud time (bottom) ---------------------
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8, 7.5), height_ratios=(2.2, 1), sharex=True, facecolor=SURFACE)
for a_ in (ax, ax2):
    a_.set_facecolor(SURFACE)
    a_.tick_params(colors=MUTED)
    for spine in a_.spines.values():
        spine.set_color(MUTED)
    a_.grid(True, color=MUTED, alpha=0.15)
    a_.axvline(p[1], color=MUTED, lw=1, ls=":")
xx = np.linspace(positions.min() - 1, positions.max() + 1, 400)
in_shot = np.abs(xs) <= FIT_MAX_MM
dn_all = 100 * np.array([s["depth_norm"] for s in shots])
ax.plot(xs[in_shot], dn_all[in_shot], "o", ms=4, color=BLUE, alpha=0.3, mec="none")
ax.errorbar(positions[inside], 100 * dn_mean[inside], yerr=100 * dn_sem[inside], fmt="o", ms=8, color=BLUE,
            mec=SURFACE, mew=1.5, capsize=3, lw=1.5, label="Mean, PMT-normalised")
# Outside the lit window: hollow, not fitted
ax.errorbar(positions[~inside], 100 * dn_mean[~inside], yerr=100 * dn_sem[~inside], fmt="o", ms=8,
            mfc=SURFACE, color=MUTED, mew=1.5, capsize=3, lw=1.2, label="Outside lit window (not fitted)")
ax.plot(xx, 100 * profile(xx, *p), color=ORANGE, lw=2.2,
        label=f"Gaussian: centre {p[1]:+.1f} mm, σ {abs(p[2]):.1f} mm")
ax.axhline(0, color=MUTED, lw=1)
ax.set_ylabel("Peak dip (%)", color=INK)
ax.set_title(f"Upward cloud, iris scan ({FIRST:04d}–{LAST:04d})", color=INK)
ax.legend(frameon=False, fontsize=9, loc="upper left", labelcolor=INK)
ax2.plot(xs, [s["light"] for s in shots], "o", ms=4, color=BLUE, alpha=0.3, mec="none")
ax2.plot(positions, l_mean, "o-", ms=6, color=BLUE, lw=1.5)
ax2.axhline(0, color=MUTED, lw=1)
ax2.set_ylabel("Probe light (V)", color=INK)
ax2.set_xlabel("Stage offset (mm; + = left, facing N)", color=INK)
fig.tight_layout()
fig.savefig(FIG, dpi=130, facecolor=SURFACE)
print(f"\nfigure: {FIG}")
