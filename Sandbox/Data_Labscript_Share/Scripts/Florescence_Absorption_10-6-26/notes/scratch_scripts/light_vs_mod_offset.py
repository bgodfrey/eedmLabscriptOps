import sys
from pathlib import Path
import h5py
import numpy as np
sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot
DARK = 0.0527
day = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\16")
for folder in ("0036", "0037"):
    rows = {}
    for p in sorted((day / folder).glob("*.h5")):
        s = load_shot(p)
        off = s["mod_offset"]
        t, v = s["t_pda"], s["pda_v"] - DARK
        inw = (t > s["enable_start_s"] + 0.002) & (t < s["enable_stop_s"] - 0.002)
        # keep the window but drop the cloud region (dip): use the first 15 ms and last 25 ms of the window
        early = inw & (t < s["enable_start_s"] + 0.017)
        late = inw & (t > s["enable_stop_s"] - 0.027)
        out = (t > s["enable_stop_s"] + 0.005)          # after the jump is disabled
        rows.setdefault(off, []).append((np.median(v[early | late]), np.median(v[out])))
    print(f"{folder}: mod_offset (mV): light in jump window, light after window (V, mean +- sd), ratio in/out")
    for off in sorted(rows):
        a = np.array(rows[off])
        r = a[:, 0] / a[:, 1]
        print(f"   {1e3*off:5.1f} mV  n={len(a):2d}  in {a[:,0].mean():.3f} +- {a[:,0].std(ddof=1):.3f}   "
              f"after {a[:,1].mean():.3f} +- {a[:,1].std(ddof=1):.3f}   in/after {r.mean():.4f} +- {r.std(ddof=1):.4f}")
