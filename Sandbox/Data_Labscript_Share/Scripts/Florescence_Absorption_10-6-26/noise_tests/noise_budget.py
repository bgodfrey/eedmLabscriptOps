import sys
from pathlib import Path
import numpy as np
import h5py
from scipy import signal

FS = 100_000.0
FM = 6123.3


def band_asd(f, asd, lo, hi):
    m = (f >= lo) & (f <= hi)
    return float(np.median(asd[m]))


for path in sys.argv[1:]:
    path = Path(path)
    with h5py.File(path, "r") as h5:
        pda = h5["data/traces/TOF_florescence"][()]
        ref = h5["data/traces/RelJumpModulation"][()]
        g = dict(h5["globals"].attrs)
    v = np.asarray(pda["values"], float)
    r = np.asarray(ref["values"], float)
    t = np.arange(len(v)) / FS
    print(f"\n=== {path.name}  N={len(v)}")
    print(f"PDA min/median/max: {v.min():.4f} / {np.median(v):.4f} / {v.max():.4f} V")
    # Segments: pre-modulation (0-42 ms), enabled (44-112 ms), after enable (115-200 ms)
    segs = {"pre (0-42ms, no mod)": (0.001, 0.042),
            "enabled (45-110ms)": (0.045, 0.110),
            "post-enable (115-195ms)": (0.115, 0.195)}
    for name, (a, b) in segs.items():
        m = (t >= a) & (t < b)
        if m.sum() < 1000:
            continue
        x = v[m]
        f, p = signal.welch(x, fs=FS, nperseg=1024, detrend="linear")
        asd = np.sqrt(p)
        rr = r[m]
        print(f"  {name}: mean {x.mean():.4f} V, std {1e3*x.std():.3f} mV, ref std {1e3*rr.std():.2f} mV")
        print(f"     ASD [uV/rtHz]: 0.3-1k {1e6*band_asd(f,asd,300,1000):.2f}, "
              f"near fm (5.5-6.0k) {1e6*band_asd(f,asd,5500,5900):.2f}, "
              f"@fm {1e6*asd[np.argmin(abs(f-FM))]:.2f}, "
              f"@2fm {1e6*asd[np.argmin(abs(f-2*FM))]:.2f}, "
              f"20-45k {1e6*band_asd(f,asd,20000,45000):.2f}")
    # list a few gain-ish globals
    keys = [k for k in g if any(s in k.lower() for s in ("gain", "pda", "probe", "power", "shutter", "detun", "jump", "mod"))]
    for k in keys:
        print(f"  global {k} = {g[k]}")
