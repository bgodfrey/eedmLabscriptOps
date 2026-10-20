import sys
from pathlib import Path
import h5py
import numpy as np

KEYS = ("mod_amp", "mod_freq", "mod_offset", "pda_gain_db", "probe_power_uW",
        "t_PMT_delay", "v_laser_jump_rel")

for p in map(Path, sys.argv[1:]):
    with h5py.File(p, "r") as h5:
        g = h5["globals"].attrs
        print(f"\n=== {p.name}")
        print("  globals:", {k: g[k] for k in KEYS if k in g})
        # Units/expansion of the test group, as runmanager stored them.
        grp = h5["globals"].get("absorption_tests")
        if grp is not None:
            print("  absorption_tests group attrs:", dict(grp.attrs))
            for sub in grp:
                print(f"    {sub}:", dict(grp[sub].attrs))
        ao = h5["devices/NI6738/ANALOG_OUTS"][()] if "devices/NI6738/ANALOG_OUTS" in h5 else None
        if ao is not None:
            names = ao.dtype.names
            col = [n for n in names if n.endswith("ao3")]
            if col:
                v = ao[col[0]]
                print(f"  NI6738 {col[0]} programmed: min {v.min():.4f}  max {v.max():.4f}  n={len(v)}")
        if "data/traces" not in h5:
            print("  no traces")
            continue
        pda = h5["data/traces/TOF_florescence"][()]
        ref = h5["data/traces/RelJumpModulation"][()]
        pmt = h5["data/traces/photon_counts"][()]["values"]
        t = pda["t"] - pda["t"][0]
        fs = 1 / np.median(np.diff(t))
        v = pda["values"]
        r = ref["values"]
        print(f"  fs {fs/1e3:.0f}k  PDA mean {v.mean():.3f} V std {1e3*v.std():.2f} mV  min {v.min():.3f} max {v.max():.3f}")
        for a, b in ((0.0, 0.04), (0.045, 0.11), (0.115, 0.14), (0.15, 0.2)):
            m = (t >= a) & (t < b)
            if m.any():
                print(f"    {1e3*a:5.0f}-{1e3*b:3.0f} ms: PDA {v[m].mean():.3f} V, ref mean {1e3*r[m].mean():7.2f} mV, "
                      f"ref p2p {1e3*np.ptp(r[m]):7.2f} mV")
        print(f"  PMT total counts {int(pmt[-1]) if len(pmt) else 'empty'}")
