"""What distinguishes the low-dip ("dropout") shots in 2026-10-16 0036/0037?

Usage: python dropout_check.py <helper dir> <dark mV> <folder> [<folder> ...]

Per shot, in time order: mod_offset, the dip D (fixed scan shape), PMT integral F, PMT cloud time
t_c (from the analysis-window origin) and transit width, the recorded jump-modulation drive
(RelJumpModulation, mean over the enabled window: does the commanded offset arrive?), the dip's
free-fit centre and width, the probe light and the MOT signal.
A frequency excursion scales D and F alike and leaves the timing and dip shape alone; a launch
problem (fewer or slower atoms) can move t_c and change the widths.
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
DARK = float(sys.argv[2]) * 1e-3
DAY = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\16")
OFFSET, WIDTH = -0.13e-3, 1.48e-3

rows = []
for folder in sys.argv[3:]:
    for p in sorted((DAY / folder).glob("*.h5"), key=lambda q: q.stat().st_mtime):
        s = load_shot(p)
        with h5py.File(p, "r") as h5:
            off = float(h5["globals"].attrs["mod_offset"])
        pmt = cse.pmt_cloud_fit(s, 6125.0)
        tc = pmt["center_s"]
        k = int(round(s["fs_pda"] / 10e3))
        m = len(s["pda_v"]) // k * k
        tb = s["t_pda"][:m].reshape(-1, k).mean(axis=1)
        vb = s["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
        win = (tb > s["enable_start_s"] + 0.002) & (tb < s["enable_stop_s"] - 0.002)
        tw, vw = tb[win], vb[win]

        def model(t, c0, c1, c2, d, cen, wid):
            u = t - tc
            return (c0 + c1 * u + c2 * u ** 2) * (1 - d * np.exp(-0.5 * ((t - cen) / wid) ** 2))

        pa, _ = optimize.curve_fit(lambda t, c0, c1, c2, d: model(t, c0, c1, c2, d, tc + OFFSET, WIDTH),
                                   tw, vw, p0=[np.median(vw), 0, 0, 0.03])
        try:
            pb, _ = optimize.curve_fit(model, tw, vw, p0=list(pa) + [tc + OFFSET, WIDTH],
                                       bounds=([-np.inf] * 3 + [-0.5, tc - 0.005, 0.3e-3],
                                               [np.inf] * 3 + [0.5, tc + 0.005, 5e-3]))
        except (RuntimeError, ValueError):
            pb = [np.nan] * 6
        # Recorded jump drive inside the enabled window
        en = (s["t_pda"] > s["enable_start_s"] + 0.002) & (s["t_pda"] < s["enable_stop_s"] - 0.002)
        ref = s["reference_v"][en]
        try:
            mot = fit_mot_loading(load_mot_trace(p))["atoms_v_at_coils_off"]
        except Exception:
            mot = np.nan
        rows.append((folder, off, pa[3], pmt["amp_cps"] * pmt["sigma_s"], tc, pmt["sigma_s"],
                     ref.mean(), ref.std(), pb[4] - tc, pb[5], pa[0], mot))

print("fold  off mV   dip %   PMT   t_c ms  PMT sig  drive mV (sd)  dip ctr  dip wid  light V  MOT mV")
for r in rows:
    print(f"{r[0]}  {1e3*r[1]:5.1f}  {100*r[2]:6.2f}  {r[3]:6.0f}  {1e3*r[4]:7.2f}  {1e3*r[5]:5.2f}   "
          f"{1e3*r[6]:6.2f} ({1e3*r[7]:.2f})  {1e3*r[8]:+6.2f}  {1e3*r[9]:6.2f}   {r[10]:6.3f}  {1e3*r[11]:5.1f}")

# Dropouts: dip below 60% of the median at their offset (line-centre region only)
a = np.array([[r[1], r[2], r[4], r[5], r[6], r[8], r[9]] for r in rows])
sel = a[:, 0] > 0.006
med = {o: np.median(a[sel & (a[:, 0] == o), 1]) for o in np.unique(a[sel, 0])}
low = np.array([sel[i] and a[i, 1] < 0.6 * med[a[i, 0]] for i in range(len(a))])
ok = sel & ~low
print(f"\n{low.sum()} dropouts of {sel.sum()} shots at offsets >= 6 mV")
for name, j, scale in (("t_c ms", 2, 1e3), ("PMT sigma ms", 3, 1e3), ("drive - offset mV", 4, 1e3),
                       ("dip centre ms", 5, 1e3), ("dip width ms", 6, 1e3)):
    vals = a[:, j] - (a[:, 0] if j == 4 else 0)
    print(f"  {name:18s}: normal {scale*np.nanmean(vals[ok]):+.3f} (sd {scale*np.nanstd(vals[ok]):.3f}), "
          f"dropout {scale*np.nanmean(vals[low]):+.3f} (sd {scale*np.nanstd(vals[low]):.3f})")
