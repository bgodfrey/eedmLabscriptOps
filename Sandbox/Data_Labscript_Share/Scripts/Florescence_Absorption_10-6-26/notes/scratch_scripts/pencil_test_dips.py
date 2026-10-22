"""Absorption dip at the PMT cloud time for single test shots (2026-10-15 0007-0009).

Usage: python pencil_test_dips.py <helper dir> <day dir> <dark mV> <folder> [<folder> ...]

Dip depth from a zero-mean Gaussian template (sigma 1.3 ms) centred on the PMT peak, as a voltage
and as a fraction of the light above the dark level. Also the light level itself and the PMT peak.
"""
import sys
from pathlib import Path
import numpy as np
from scipy import signal

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot

DAY, DARK = Path(sys.argv[2]), float(sys.argv[3]) * 1e-3
SIGMA_T = 1.3e-3
for folder in sys.argv[4:]:
    for p in sorted((DAY / folder).glob("*.h5")):
        s = load_shot(p)
        fs, v = s["fs_pda"], s["pda_v"]
        # PMT peak from a 1 ms boxcar of the count rate
        k = max(1, int(1e-3 * s["fs_pmt"]))
        pmt = np.convolve(s["pmt_rate"], np.ones(k) / k, mode="same")
        t0 = s["t_pmt"][np.argmax(pmt)]
        half = int(4 * SIGMA_T * fs)
        i0 = np.argmin(abs(s["t_pda"] - t0))
        tt = np.arange(-half, half + 1) / fs
        g = np.exp(-0.5 * (tt / SIGMA_T) ** 2)
        g0 = g - g.mean()
        seg = v[i0 - half:i0 + half + 1]
        dip = -np.dot(seg, g0) / np.dot(g0, g)        # positive = less light at the cloud
        light = np.median(v) - DARK
        print(f"{folder} {p.name[-8:]}: light above dark {1e3*light:8.1f} mV, PMT peak {1e3*t0:6.1f} ms "
              f"({pmt.max()/1e6:.2f} Mcps), dip {1e3*dip:6.2f} mV = {100*dip/light:5.2f}% of the light")
