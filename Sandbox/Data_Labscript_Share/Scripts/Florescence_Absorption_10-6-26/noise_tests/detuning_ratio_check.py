"""Dip and PMT fluorescence at two probe offsets, interleaved, iris fixed (2026-10-16 0036).

Usage: python detuning_ratio_check.py <helper dir> <shot folder> <dark mV>

Purpose: the iris atom scan (0023-0035) ran at mod_offset 0 instead of line centre (9.5 mV). This
folder alternates the two offsets in shuffled order at iris 0 mm, everything else as in the scan,
so their ratio converts the scan to line centre.

Per shot (same dip model as iris_atom_scan.py):
1. PMT Gaussian fit -> cloud time t_c and the fluorescence integral F = amplitude x transit width.
2. PDA at 10 kHz in the jump-enabled window, dark subtracted, fitted as
   quadratic baseline * (1 - D exp(-(t - t0)^2 / (2 s_t^2))).
   Pass A: shape fixed to the scan's common shape (t0 = t_c - 0.13 ms, s_t = 1.48 ms), the same
   estimator as the scan, so the ratio applies to it directly.
   Pass B: centre and width free, to check the shape does not change with detuning.
3. MOT fluorescence at coils-off as the atom-number monitor.

Ratios (line centre / offset 0):
- dip ratio from raw D (not PMT-normalised: the PMT carries the same detuning factor and would
  cancel it). Error from the two standard errors.
- PMT ratio: fluorescence also scales with the scattering rate per atom, so it is an independent
  estimate of the same factor (both at low saturation).
- Shot-to-shot scatter of D and F at each setting: probe-frequency jitter is first order off
  resonance and second order at line centre, launch variation is the same at both.
- PMT dead time: the H7421 loses counts at Mcps rates (notebook 2026-10-15: about 67 ns,
  non-paralysable, from "linear to 1.5 Mcps"). The brighter line-centre shots lose more, which
  lowers the PMT ratio and widens the fitted transit. F_dt corrects the smoothed rate,
  n = m / (1 - m tau), then integrates it above the background over +-8 ms of the cloud.
  The correction grows steeply near 1/tau = 15 Mcps, so it only indicates the size of the effect.
"""
import sys
import warnings
from pathlib import Path
import h5py
import numpy as np
from scipy import optimize

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot, load_mot_trace, fit_mot_loading
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
FOLDER, DARK = Path(sys.argv[2]), float(sys.argv[3]) * 1e-3
OFFSET, WIDTH = -0.13e-3, 1.48e-3          # common dip shape from the scan (iris_atom_scan.py)
TAU_DT = 67e-9                              # PMT dead time (s), notebook estimate, not measured


def pmt_integral_deadtime(shot, tc, tau=TAU_DT, smooth_s=0.2e-3, half_s=0.008):
    """Counts from the cloud, with the measured rate corrected for a non-paralysable dead time.

    The rate is smoothed first (0.2 ms boxcar, short against the 2 ms transit) because the
    correction is nonlinear and would amplify the counting noise of single 10 us bins.
    Background = median corrected rate 8-30 ms from the cloud centre.
    """
    t, m = shot["t_pmt"], shot["pmt_rate"]
    k = max(1, int(round(smooth_s * shot["fs_pmt"])))
    ms = np.convolve(m, np.ones(k) / k, mode="same")
    n = ms / np.clip(1 - ms * tau, 0.05, None)          # clip: guard against rates near 1/tau
    u = t - tc
    bg = np.median(n[(np.abs(u) > half_s) & (np.abs(u) < 0.030)])
    sel = np.abs(u) <= half_s
    return np.sum(n[sel] - bg) / shot["fs_pmt"], ms[sel].max()


def dip_model(t, tc, c0, c1, c2, depth, centre, width):
    """Quadratic baseline (in t - tc) times a Gaussian dip."""
    u = t - tc
    return (c0 + c1 * u + c2 * u ** 2) * (1 - depth * np.exp(-0.5 * ((t - centre) / width) ** 2))


rows = []
for p in sorted(FOLDER.glob("*.h5"), key=lambda q: q.stat().st_mtime):
    shot = load_shot(p)
    with h5py.File(p, "r") as h5:
        g = h5["globals"].attrs
        off, x = float(g["mod_offset"]), float(g["iris_x_mm"])
    if shot["pda_v"].max() >= 9.95:
        print(f"{p.name}: PDA clipped, skipped")
        continue
    pmt = cse.pmt_cloud_fit(shot, shot["mod_freq"] or 6125.0)      # f only sets the smoothing
    if not pmt["ok"]:
        print(f"{p.name}: no PMT fit, skipped")
        continue
    # 10 kHz blocks in the enabled window, dark subtracted
    k = int(round(shot["fs_pda"] / 10e3))
    m = len(shot["pda_v"]) // k * k
    tb = shot["t_pda"][:m].reshape(-1, k).mean(axis=1)
    vb = shot["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
    win = (tb > shot["enable_start_s"] + 0.002) & (tb < shot["enable_stop_s"] - 0.002)
    tw, vw, tc = tb[win], vb[win], pmt["center_s"]
    # Pass A: common shape
    pa, ca = optimize.curve_fit(lambda t, c0, c1, c2, d: dip_model(t, tc, c0, c1, c2, d, tc + OFFSET, WIDTH),
                                tw, vw, p0=[np.median(vw), 0, 0, 0.01])
    # Pass B: free centre and width (bounded to a sensible range)
    try:
        pb, _ = optimize.curve_fit(lambda t, *a: dip_model(t, tc, *a), tw, vw,
                                   p0=list(pa) + [tc + OFFSET, WIDTH],
                                   bounds=([-np.inf] * 3 + [-0.5, tc - 0.005, 0.3e-3],
                                           [np.inf] * 3 + [0.5, tc + 0.005, 5e-3]))
        free = (pb[3], pb[4] - tc, pb[5])
    except (RuntimeError, ValueError):
        free = (np.nan, np.nan, np.nan)
    try:
        mot = fit_mot_loading(load_mot_trace(p))["atoms_v_at_coils_off"]
    except Exception:
        mot = np.nan
    f_dt, peak_cps = pmt_integral_deadtime(shot, tc)
    f_raw, _ = pmt_integral_deadtime(shot, tc, tau=0.0)          # same sum, no correction
    rows.append({"name": p.stem, "off": off, "x": x, "D": pa[3], "D_err": np.sqrt(ca[3, 3]), "light": pa[0],
                 "F": pmt["amp_cps"] * pmt["sigma_s"], "pmt_sigma": pmt["sigma_s"], "mot": mot,
                 "free": free, "F_sum": f_raw, "F_dt": f_dt, "peak_cps": peak_cps})

print(f"{len(rows)} shots, iris at {sorted({r['x'] for r in rows})} mm\n")
print(" # offset mV  dip %  (fit err)  light V  PMT counts  MOT mV")
for i, r in enumerate(rows):
    print(f"{i:2d}  {1e3*r['off']:6.1f}   {100*r['D']:5.2f}  ({100*r['D_err']:.2f})   {r['light']:6.3f}   "
          f"{r['F']:8.0f}   {1e3*r['mot']:5.1f}")

summary = {}
print("\noffset   shots  dip %  (sem)  rel sd   PMT counts  (sem)  rel sd   MOT mV   free width ms  free centre ms")
for off in sorted({r["off"] for r in rows}):
    g = [r for r in rows if r["off"] == off]
    d = np.array([r["D"] for r in g])
    f = np.array([r["F"] for r in g])
    mo = np.array([r["mot"] for r in g])
    fw = np.array([r["free"][2] for r in g])
    fc = np.array([r["free"][1] for r in g])
    n = len(g)
    summary[off] = (d.mean(), d.std(ddof=1) / np.sqrt(n), f.mean(), f.std(ddof=1) / np.sqrt(n))
    print(f"{1e3*off:5.1f} mV  {n:4d}  {100*d.mean():5.2f} ({100*d.std(ddof=1)/np.sqrt(n):.2f})  "
          f"{100*d.std(ddof=1)/d.mean():4.0f}%   {f.mean():8.0f} ({f.std(ddof=1)/np.sqrt(n):4.0f})  "
          f"{100*f.std(ddof=1)/f.mean():4.0f}%   {1e3*np.nanmean(mo):5.1f}   "
          f"{1e3*np.nanmedian(fw):5.2f}          {1e3*np.nanmedian(fc):+5.2f}")
    print(f"           corr(dip, PMT) {np.corrcoef(d, f)[0,1]:+.2f}, corr(dip, MOT) {np.corrcoef(d, mo)[0,1]:+.2f}, "
          f"PMT transit sigma {1e3*np.median([r['pmt_sigma'] for r in g]):.2f} ms, "
          f"measured peak {np.median([r['peak_cps'] for r in g])/1e6:.1f} Mcps")
    summary[off] += (np.mean([r["F_sum"] for r in g]), np.mean([r["F_dt"] for r in g]))

offs = sorted(summary)
if len(offs) == 2:
    (d0, e0, f0, g0, s0, t0), (d1, e1, f1, g1, s1, t1) = summary[offs[0]], summary[offs[1]]
    print(f"\nPMT summed counts: {s0:.0f} -> {s1:.0f} (ratio {s1/s0:.2f}); dead-time corrected "
          f"({1e9*TAU_DT:.0f} ns): {t0:.0f} -> {t1:.0f} (ratio {t1/t0:.2f})")
    rd = d1 / d0
    rd_err = rd * np.sqrt((e0 / d0) ** 2 + (e1 / d1) ** 2)
    rf = f1 / f0
    rf_err = rf * np.sqrt((g0 / f0) ** 2 + (g1 / f1) ** 2)
    print(f"\n{1e3*offs[1]:.1f} mV / {1e3*offs[0]:.1f} mV: dip ratio {rd:.2f} +- {rd_err:.2f}, "
          f"PMT ratio {rf:.2f} +- {rf_err:.2f}")
    # Lorentzian reading of the dip ratio: L = 1/(1 + x^2), x = 2 delta / Gamma_eff, line centre = 1
    x_eff = np.sqrt(max(rd - 1, 0))
    print(f"if 9.5 mV is line centre: offset 0 sits at 2|delta|/Gamma = {x_eff:.2f} (Lorentzian)")
