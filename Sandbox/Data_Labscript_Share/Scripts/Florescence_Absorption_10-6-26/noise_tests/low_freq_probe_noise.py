"""How much of the low-frequency (1/f-like) noise comes from the probe beam, and is it steady?

Usage: python low_freq_probe_noise.py <helper dir> <dark folder> <lit day dir>

Compares probe-on, no-launch shots (2026-10-14 0059-0069) with photodetector-blocked shots
(2026-10-15 0005) using 25 ms windows (40 Hz resolution), stepped by 12.5 ms along the record.
For each band the probe's share of the noise power is (lit^2 - dark^2) / lit^2.

Two readings are possible for light-induced low-frequency noise:
  - it decays along the record -> the tail of the probe-shutter transient (mechanical ringing);
  - it is flat along the record -> steady noise carried by the probe beam (laser, pointing, etalons).
Windows that straddle the modulation switches (enable on/off, sine end) are skipped, because a step in
the probe frequency can step the transmitted power and leak into the low-frequency bins.
"""
import sys
import warnings
from pathlib import Path
import numpy as np
from scipy import signal

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot

warnings.filterwarnings("ignore")
DARK, LIT_DAY = Path(sys.argv[2]), Path(sys.argv[3])
WIN_S, STEP_S = 0.025, 0.0125
BANDS = {"40-100 Hz": (40, 100), "120 Hz line": (110, 130), "140-400 Hz": (140, 400), "0.4-1 kHz": (400, 1000)}


def load(folders):
    out = []
    for folder in folders:
        for p in sorted(folder.glob("*.h5")):
            try:
                out.append(load_shot(p))
            except Exception:
                pass
    return out


def band_vs_time(shots):
    fs = shots[0]["fs_pda"]
    n, step = int(WIN_S * fs), int(STEP_S * fs)
    starts = np.arange(0, len(shots[0]["pda_v"]) - n + 1, step)
    hann = signal.windows.hann(n)
    norm = fs * np.sum(hann ** 2)
    f = np.fft.rfftfreq(n, 1 / fs)
    acc = {b: np.zeros(len(starts)) for b in BANDS}
    for s in shots:
        for i, a in enumerate(starts):
            psd = 2 * np.abs(np.fft.rfft(signal.detrend(s["pda_v"][a:a + n]) * hann)) ** 2 / norm
            for b, (lo, hi) in BANDS.items():
                m = (f >= lo) & (f < hi)
                # the line band uses the mean (its power); broadband uses the median (ignores lines)
                acc[b][i] += (psd[m].mean() if "line" in b else np.median(psd[m])) / len(shots)
    sine0 = np.median([s["sine_start_s"] for s in shots])
    switches = (sine0, sine0 + 0.070, sine0 + 0.100)
    ok = np.array([not any(a / fs < sw < (a + n) / fs for sw in switches) for a in starts])
    return (starts + n / 2) / fs, acc, ok


dark = load([DARK])
lit = load([LIT_DAY / f"{i:04d}" for i in range(59, 70)])
tc, acc_l, ok_l = band_vs_time(lit)
_, acc_d, ok_d = band_vs_time(dark)
ok = ok_l & ok_d
print(f"lit {len(lit)} shots, dark {len(dark)} shots; 25 ms windows (40 Hz bins); * = straddles a switch (skipped below)")
print(f"{'centre ms':>9} " + " ".join(f"{b:>22s}" for b in BANDS))
print(f"{'':>9} " + " ".join(f"{'lit / dark (probe %)':>22s}" for _ in BANDS))
for i, t in enumerate(tc):
    cells = []
    for b in BANDS:
        l, d = np.sqrt(acc_l[b][i]), np.sqrt(acc_d[b][i])
        share = max(l ** 2 - d ** 2, 0) / l ** 2
        cells.append(f"{1e6*l:7.1f} /{1e6*d:6.1f} ({100*share:3.0f}%)")
    print(f"{1e3*t:8.1f}{'*' if not ok[i] else ' '} " + " ".join(f"{c:>22s}" for c in cells))

late = ok & (tc > 0.100)
print("\nsteady (windows centred after 100 ms, no switch):")
for b in BANDS:
    l, d = np.sqrt(acc_l[b][late].mean()), np.sqrt(acc_d[b][late].mean())
    print(f"  {b:12s}: lit {1e6*l:6.1f}, blocked {1e6*d:6.1f} uV/rtHz -> probe beam carries "
          f"{100*max(l**2-d**2, 0)/l**2:3.0f}% of the noise power; probe part {1e6*np.sqrt(max(l**2-d**2,0)):6.1f} uV/rtHz "
          f"= {np.sqrt(max(l**2-d**2,0))/np.mean([s['pda_v'].mean() for s in lit]):.1e} /rtHz relative")
