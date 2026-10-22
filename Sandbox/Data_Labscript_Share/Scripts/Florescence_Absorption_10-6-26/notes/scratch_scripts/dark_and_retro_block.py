"""PDA levels for the dark shots (0009, probe blocked, iris closed) and the retro-blocked shots
(0010, mirror blocked, iris 1 mm), 2026-10-15. Prints the relevant globals and the PDA mean and
spread per file, over the whole 0.2 s record and over its last 100 ms (well after the shutter)."""
import sys
from pathlib import Path
import h5py
import numpy as np

DAY = Path(sys.argv[1])
G = ("t_PMT_delay", "MRR_TRIG_BOOL", "mod_amp", "mod_offset", "pda_gain_db", "probe_power_uW", "iris_x_mm")
for folder in sys.argv[2:]:
    print(f"\n{folder}")
    for p in sorted((DAY / folder).glob("*.h5")):
        with h5py.File(p, "r") as h5:
            g = {k: h5["globals"].attrs.get(k) for k in G}
            v = h5["data/traces/TOF_florescence"][()]["values"]
        late = v[len(v) // 2:]
        print(f"  {p.name[-8:]}: mean {1e3*v.mean():8.3f} mV, last 100 ms {1e3*late.mean():8.3f} mV, "
              f"std {1e3*v.std():6.3f} mV, min {1e3*v.min():7.2f}, max {1e3*v.max():7.2f}  "
              + ", ".join(f"{k}={g[k]}" for k in G))

