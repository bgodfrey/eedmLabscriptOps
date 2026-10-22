"""Noise floor for the pencil-probe dip measurement at 70 dB, from the dark shots.

Usage: python dark_dip_noise_70db.py <day dir> <folder> [<folder> ...]

For each folder: band ASDs of the PDA record (Welch, 4096 points), and the noise on the amplitude
of a Gaussian absorption dip (sigma 1.3 ms, the upward width from 0130) fitted at a known time.
The dip estimate uses a zero-mean Gaussian template over +-4 sigma, so a constant offset in the
window drops out, as the per-shot baseline does in the real fit. Its standard deviation, taken
over template positions stepped by 2 ms along each record and pooled over shots, is the per-shot
amplitude noise in mV.
"""
import sys
from pathlib import Path
import h5py
import numpy as np
from scipy import signal

DAY = Path(sys.argv[1])
SIGMA_T = 1.3e-3
BANDS = [(40, 100), (100, 400), (400, 1000), (1000, 3000)]
for folder in sys.argv[2:]:
    psds, amps = [], []
    for p in sorted((DAY / folder).glob("*.h5")):
        with h5py.File(p, "r") as h5:
            tr = h5["data/traces/TOF_florescence"][()]
        t, v = tr["t"] - tr["t"][0], tr["values"]
        fs = 1 / np.median(np.diff(t))
        f, pxx = signal.welch(v, fs=fs, nperseg=4096, noverlap=2048, detrend="linear")
        psds.append(pxx)
        # zero-mean Gaussian template; amplitude = projection / (template . dip shape)
        half = int(4 * SIGMA_T * fs)
        tt = np.arange(-half, half + 1) / fs
        g = np.exp(-0.5 * (tt / SIGMA_T) ** 2)
        g0 = g - g.mean()
        scale = np.dot(g0, g)            # response of the projection to a unit-depth dip
        for a in range(0, len(v) - len(g), int(2e-3 * fs)):
            amps.append(np.dot(v[a:a + len(g)], g0) / scale)
    asd = np.sqrt(np.mean(psds, axis=0))
    print(f"{folder}: {len(psds)} shots, fs {fs/1e3:.0f} kS/s; band ASD (uV/rtHz): " + ", ".join(
        f"{lo}-{hi} Hz {1e6*np.median(asd[(f>=lo)&(f<hi)]):.0f}" for lo, hi in BANDS))
    print(f"   Gaussian-dip amplitude noise (sigma {1e3*SIGMA_T:.1f} ms): {1e3*np.std(amps):.2f} mV per shot")
