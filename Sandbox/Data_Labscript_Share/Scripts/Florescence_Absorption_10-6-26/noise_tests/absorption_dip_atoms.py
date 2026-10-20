"""Atom number from the direct (DC) absorption dip on the upward pass, modulation off.

Usage: python absorption_dip_atoms.py <helper dir> <shot folder>

Per shot:
  1. PMT Gaussian fit -> cloud centre t_c and transit width sigma (independent of the photodiode).
  2. Photodiode block-averaged to 10 kHz inside the modulation-enabled window (the frequency jump
     is on there, so the probe sits at the chosen offset), excluding 2 ms at each edge.
  3. Fit  V(t) = quadratic baseline * (1 - D * exp(-(t - t_c)^2 / (2 sigma^2)))  with t_c, sigma
     fixed from the PMT. D is the peak fractional dip of the detected (double-pass) probe power.
     A second fit with centre and width free checks that the dip sits where the PMT says.
  4. Time integral of the fractional dip: I = D * sqrt(2 pi) * sigma.

Conversion to atoms (optically thin, atoms crossing a horizontal Gaussian sheet vertically):
  single-pass fractional absorption D1(t) = sigma_eff * sum_i I(x_i, y_i) / P, with
  I/P = 2/(pi w_x w_y) exp(-2x^2/w_x^2 - 2y^2/w_y^2). Integrating over the vertical transit at
  speed v removes the beam height w_y:
      integral D1 dt = sigma_eff * sqrt(2/pi) / (w_x v) * sum_i exp(-2 x_i^2 / w_x^2)
  For a cloud of horizontal rms size s_x centred in the sheet, the sum is N * (1 + 4 s_x^2/w_x^2)^-1/2.
  The retroreflected probe passes the atoms twice, so integral D1 dt = I / 2.
  sigma_eff = sigma_0 / (1 + s) on resonance, s = total intensity / I_sat.
"""
import sys
import warnings
from pathlib import Path
import numpy as np
from scipy import optimize

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot, load_mot_trace, fit_mot_loading
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
folder = Path(sys.argv[2])
files = sorted(folder.glob("*.h5"), key=lambda p: int(p.stem.split("_")[-1]))

rows = []
for p in files:
    shot = load_shot(p)
    if shot["pda_v"].max() >= 9.95:
        print(f"{p.name}: clipped, skipped")
        continue
    pmt = cse.pmt_cloud_fit(shot, shot["mod_freq"] or 6125.0)
    if not pmt["ok"]:
        print(f"{p.name}: no PMT fit, skipped")
        continue
    t, v = shot["t_pda"], shot["pda_v"]
    fs = shot["fs_pda"]
    n = int(round(fs / 10e3))                     # block-average to 10 kHz
    m_len = len(v) // n * n
    tb = t[:m_len].reshape(-1, n).mean(axis=1)
    vb = v[:m_len].reshape(-1, n).mean(axis=1)
    win = (tb > shot["enable_start_s"] + 0.002) & (tb < shot["enable_stop_s"] - 0.002)
    tw, vw = tb[win], vb[win]
    tc, sg = pmt["center_s"], pmt["sigma_s"]
    u = tw - tc                                   # centred time for a well-conditioned baseline

    def fixed(tt, c0, c1, c2, depth):
        x = tt - tc
        return (c0 + c1 * x + c2 * x ** 2) * (1 - depth * np.exp(-0.5 * (x / sg) ** 2))

    def free(tt, c0, c1, c2, depth, centre, width):
        x = tt - tc
        return (c0 + c1 * x + c2 * x ** 2) * (1 - depth * np.exp(-0.5 * ((tt - centre) / width) ** 2))

    p0 = [np.median(vw), 0.0, 0.0, 0.005]
    pf, cov = optimize.curve_fit(fixed, tw, vw, p0=p0)
    depth, depth_err = pf[3], np.sqrt(cov[3, 3])
    resid = vw - fixed(tw, *pf)
    try:
        pg, _ = optimize.curve_fit(free, tw, vw, p0=list(pf) + [tc, sg], maxfev=20000)
        centre_free, width_free = pg[4], abs(pg[5])
    except RuntimeError:
        centre_free, width_free = np.nan, np.nan
    mot = load_mot_trace(p)
    mot_mv = 1e3 * fit_mot_loading(mot)["atoms_v_at_coils_off"] if mot is not None else np.nan
    rows.append({"file": p.name, "depth": depth, "depth_err": depth_err, "sigma_s": sg, "tc": tc,
                 "baseline_v": pf[0], "resid_mv": 1e3 * resid.std(),
                 "centre_shift_ms": 1e3 * (centre_free - tc), "width_free_ms": 1e3 * width_free,
                 "pmt_kcps": pmt["amp_cps"] / 1e3, "mot_mv": mot_mv,
                 "integral_s": depth * np.sqrt(2 * np.pi) * sg})

print(f"{'run':>4} {'dip %':>7} {'+/-':>6} {'base V':>7} {'resid mV':>8} {'PMT sig ms':>10} {'free width':>10} "
      f"{'centre shift':>12} {'PMT kcps':>8} {'MOT mV':>6}")
for r in rows:
    print(f"{r['file'].rsplit('_', 1)[1][:-3]:>4} {100*r['depth']:7.3f} {100*r['depth_err']:6.3f} {r['baseline_v']:7.2f} "
          f"{r['resid_mv']:8.2f} {1e3*r['sigma_s']:10.2f} {r['width_free_ms']:10.2f} {r['centre_shift_ms']:12.2f} "
          f"{r['pmt_kcps']:8.0f} {r['mot_mv']:6.1f}")
depth = np.array([r["depth"] for r in rows])
integral = np.array([r["integral_s"] for r in rows])
pmt_k = np.array([r["pmt_kcps"] for r in rows])
print(f"\n{len(rows)} shots: peak dip {100*depth.mean():.3f} +/- {100*depth.std(ddof=1)/np.sqrt(len(depth)):.3f} % "
      f"(shot-to-shot spread {100*depth.std(ddof=1)/depth.mean():.0f}%), corr(dip, PMT amplitude) = {np.corrcoef(depth, pmt_k)[0,1]:+.2f}")
print(f"time-integrated dip I = {1e6*integral.mean():.2f} +/- {1e6*integral.std(ddof=1)/np.sqrt(len(integral)):.2f} us")

# ---- conversion to atoms ------------------------------------------------------------------
v_probe = 2.30          # m/s at the probe, from the up/down crossing times
w_x = 9e-3              # 1/e^2 half-width of the 18 mm sheet (m); user estimate
w_y = 0.5e-3            # 1/e^2 half-height of the 1 mm sheet (m); drops out of the result
s_x = 4.0e-3            # horizontal rms cloud size on the upward pass (m), assumed = vertical
P_atoms = 80e-6         # probe power per pass at the atoms (W)
peak_I_single = 2 * P_atoms / (np.pi * w_x * w_y) * 1e3 / 1e4   # mW/cm^2
I_total = 2 * peak_I_single                                    # two passes
geometry = 1 / np.sqrt(1 + 4 * s_x ** 2 / w_x ** 2)
print(f"\npeak intensity {peak_I_single:.2f} mW/cm^2 per pass, {I_total:.2f} mW/cm^2 both passes; "
      f"geometry factor {geometry:.3f}")
for label, sigma0_cm2, isat in (("cycling (sigma+ on F=4->5')", 3.468e-9, 1.102),
                                ("pi-polarized", 2.307e-9, 1.657),
                                ("isotropic", 1.405e-9, 2.702)):
    s = I_total / isat
    sigma_eff = sigma0_cm2 * 1e-4 / (1 + s)
    n_atoms = (integral.mean() / 2) * w_x * v_probe / (sigma_eff * np.sqrt(2 / np.pi) * geometry)
    print(f"  {label:28s}: s = {s:4.2f}, sigma_eff = {sigma_eff*1e4:.2e} cm^2 -> N_upward = {n_atoms:.2e}")
