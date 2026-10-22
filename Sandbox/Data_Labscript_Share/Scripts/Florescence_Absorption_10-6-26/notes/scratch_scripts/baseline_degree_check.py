import sys
import numpy as np
sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse
p = r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\17\0003\2026-10-17_0003_exp12_switchArmsAndLaunch_Absorption_6.h5"
shot = load_shot(p)
pmt = cse.pmt_cloud_fit(shot, shot["mod_freq"] or 6125.0)
k = int(round(shot["fs_pda"] / 10e3)); m = len(shot["pda_v"]) // k * k
tb = shot["t_pda"][:m].reshape(-1, k).mean(axis=1)
vb = shot["pda_v"][:m].reshape(-1, k).mean(axis=1) - 0.0527
win = (tb > shot["enable_start_s"] + 0.002) & (tb < shot["enable_stop_s"] - 0.002)
u, v = tb[win] - pmt["center_s"], vb[win]
far = np.abs(u) > 0.010
b = 1e-3 / 2.30
pk = np.abs(u - 0.5e-3) < b / 2      # the peak 1 mm bin (centre +0.5 ms, as in the figure)
for deg in (1, 2):
    c = np.polyfit(u[far], v[far], deg)
    base = np.polyval(c, u)
    rms = np.std(v[far] - base[far])
    print(f"degree {deg}: baseline at cloud {np.polyval(c, 0):.4f} V, peak-bin dip {100*np.mean(1 - v[pk]/base[pk]):.2f}%, "
          f"rms residual in baseline region {1e3*rms:.2f} mV" + (f", curvature term {c[0]:.3g} V/s^2" if deg == 2 else ""))
