"""PMT fluorescence profile (dead-time corrected) against the PDA dip profile, 10/17 0 mm shots.

Usage: python pmt_vs_pda_profile.py <helper dir> <dark mV>
PMT: smoothed count rate (2 modulation periods), corrected for non-paralysable dead time
n = m / (1 - m tau), tau = 60 ns, minus a constant background (median of u < -18 ms and u > +18 ms),
averaged over shots. Its background is scattered probe light, so a probe-power hump would LOWER the
counts where the PDA shows a 'dip'; atoms raise both the dip and the counts.
PDA: fractional dip with quadratic baselines at |u| > 6 and > 10 ms.
Both normalised to their peak; 1 ms bins. No-atom shot 0000 #7 shown for the PMT too.
"""
import sys
import warnings
from pathlib import Path
import numpy as np
from scipy import signal

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DARK = float(sys.argv[2]) * 1e-3
DAY = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\17")
TAU = 60e-9
V = 2.30
edges = np.arange(-0.0165, 0.0165 + 1e-9, 1e-3)
cen = 0.5 * (edges[1:] + edges[:-1])


def binned(u, y):
    idx = np.digitize(u, edges) - 1
    return np.array([y[idx == i].mean() for i in range(len(cen))])


def pmt_profile(s, tc):
    kern, _ = cse.triangular_kernel(s["fs_pmt"], 6125.0, 2)
    m = signal.fftconvolve(s["pmt_rate"], kern, mode="same")
    n = m / (1 - m * TAU)
    u = s["t_pmt"] - tc
    bg = np.median(n[((u < -0.018) & (u > -0.025)) | ((u > 0.018) & (u < 0.040))])
    return binned(u, n - bg), bg


pda = {6: [], 10: []}
pmt, bgs = [], []
for f in ("0003", "0008", "0013"):
    for p in sorted((DAY / f).glob("*.h5")):
        s = load_shot(p)
        tc = cse.pmt_cloud_fit(s, 6125.0)["center_s"]
        k = int(round(s["fs_pda"] / 10e3))
        mm = len(s["pda_v"]) // k * k
        tb = s["t_pda"][:mm].reshape(-1, k).mean(axis=1)
        vb = s["pda_v"][:mm].reshape(-1, k).mean(axis=1) - DARK
        w = (tb > s["enable_start_s"] + 0.002) & (tb < s["enable_stop_s"] - 0.002)
        u, v = tb[w] - tc, vb[w]
        for lo in pda:
            far = np.abs(u) > lo * 1e-3
            pda[lo].append(binned(u, 1 - v / np.polyval(np.polyfit(u[far], v[far], 2), u)))
        pr, bg = pmt_profile(s, tc)
        pmt.append(pr)
        bgs.append(bg)
pm = np.mean(pmt, axis=0)
pm_sem = np.std(pmt, axis=0, ddof=1) / np.sqrt(len(pmt))
e_prof, e_bg = pmt_profile(load_shot(sorted((DAY / "0000").glob("*_7.h5"))[0]), 0.07017)
print(f"PMT background {np.mean(bgs)/1e6:.2f} Mcps (dead-time corrected); peak {pm.max()/1e6:.2f} Mcps; "
      f"no-atom shot background {e_bg/1e6:.2f} Mcps")
print("u ms        " + " ".join(f"{1e3*c:+6.1f}" for c in cen))
print("z mm        " + " ".join(f"{-1e3*V*c:+6.1f}" for c in cen))
print("PMT /peak   " + " ".join(f"{x/pm.max():+6.3f}" for x in pm))
print("  sem       " + " ".join(f"{x/pm.max():6.3f}" for x in pm_sem))
for lo in pda:
    a = np.mean(pda[lo], axis=0)
    print(f"PDA{lo:>2d} /peak " + " ".join(f"{x/a.max():+6.3f}" for x in a))
print("no-atom PMT " + " ".join(f"{x/pm.max():+6.3f}" for x in e_prof))
