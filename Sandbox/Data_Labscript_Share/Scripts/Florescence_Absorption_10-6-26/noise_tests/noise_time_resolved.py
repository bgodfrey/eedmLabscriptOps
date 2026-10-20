"""Is the low-frequency noise hump steady, or a transient after the probe shutter opens?

Usage: python noise_time_resolved.py <helper dir> <day dir> <figure.png>

The probe shutter is commanded open 7 ms before each 0.2 s PDA record starts. A mechanical transient
(shutter blade ringing, the shock shaking mounts) would be strongest early in the record and die away;
laser or electronic noise would be the same throughout. Uses the no-launch shots 0059-0069 (no atoms
anywhere in the record), so the only changes along the record are the modulation (sine 43-143 ms) and
the frequency jump (enabled 43-113 ms).

Output: band-limited noise in 8 ms windows stepped along the record, plus average spectra of three
stretches: before the modulation, during the enabled window, and after the sine ends.
"""
import sys
import warnings
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot

warnings.filterwarnings("ignore")
DAY, FIG = Path(sys.argv[2]), Path(sys.argv[3])
FOLDERS = [f"{i:04d}" for i in range(59, 70)]
WIN = 2048                 # 8.2 ms windows at 250 kS/s
BANDS = {"0.15-1 kHz": (150, 1000), "1-3 kHz": (1000, 3000), "3.17 kHz line": (3100, 3250),
         "3.5-5.5 kHz": (3500, 5500), "6.5-9 kHz": (6500, 9000)}

shots = []
for folder in FOLDERS:
    for p in sorted((DAY / folder).glob("*.h5")):
        try:
            shots.append(load_shot(p))
        except Exception:
            pass
fs = shots[0]["fs_pda"]
sine0 = np.median([s["sine_start_s"] for s in shots])
switches = [sine0, sine0 + 0.070, sine0 + 0.100]          # enable on, enable off, sine end

# --- band noise vs time: Hann window, per-window linear detrend, one window every 4 ms -------------
starts = np.arange(0, int(0.2 * fs) - WIN, WIN // 2)
centres = (starts + WIN / 2) / fs
hann = signal.windows.hann(WIN)
norm = fs * np.sum(hann ** 2)
f_win = np.fft.rfftfreq(WIN, 1 / fs)
band_t = {b: np.zeros(len(starts)) for b in BANDS}
for s in shots:
    v = s["pda_v"]
    for i, a in enumerate(starts):
        x = signal.detrend(v[a:a + WIN], type="linear")
        psd = 2 * np.abs(np.fft.rfft(x * hann)) ** 2 / norm
        for b, (lo, hi) in BANDS.items():
            m = (f_win >= lo) & (f_win < hi)
            # mean for the line band (captures its power), median for broadband (ignores lines)
            band_t[b][i] += (psd[m].mean() if "line" in b else np.median(psd[m])) / len(shots)
# a window straddling a switch sees the step itself; mark those
straddle = np.array([any(a / fs < sw < (a + WIN) / fs for sw in switches) for a in starts])

print(f"{len(shots)} no-launch shots; sine starts {1e3*sine0:.1f} ms into the record; "
      f"probe shutter commanded open at -7 ms")
print("band ASD (uV/rtHz) vs window centre (ms after record start); * = window straddles a switch")
print(f"{'t ms':>6} " + " ".join(f"{b:>13s}" for b in BANDS))
for i in range(0, len(starts), 2):
    print(f"{1e3*centres[i]:6.1f}{'*' if straddle[i] else ' '}" + " ".join(
        f"{1e6*np.sqrt(band_t[b][i]):13.2f}" for b in BANDS))

# --- average spectra of three stretches -------------------------------------------------------------
stretches = {"before modulation (10 ms to sine start)": (0.010, sine0 - 0.002),
             "enabled window (jump + sine)": (sine0 + 0.003, sine0 + 0.067),
             "after the sine ends": (sine0 + 0.103, 0.198)}
fig, ax = plt.subplots(2, 1, figsize=(10, 9))
for label, (a, b) in stretches.items():
    psds = []
    for s in shots:
        m = (s["t_pda"] >= a) & (s["t_pda"] < b)
        f, p = signal.welch(s["pda_v"][m], fs=fs, nperseg=WIN, noverlap=WIN // 2, detrend="linear")
        psds.append(p)
    asd = np.sqrt(np.mean(psds, axis=0))
    ax[0].loglog(f[1:], 1e6 * asd[1:], lw=0.9, label=label)
dc = np.mean([s["pda_v"].mean() for s in shots])
shot = 1e6 * np.sqrt(2 * 1.602e-19 * 1.51e6 * dc / (1 + (f[1:] / 11e3) ** 2))
ax[0].loglog(f[1:], shot, "k--", lw=0.8, label=f"shot noise at {dc:.1f} V")
ax[0].axvline(6125, color="gray", ls=":", lw=0.8)
ax[0].set(xlabel="Frequency (Hz)", ylabel="PDA ASD (µV/√Hz)", title="No-launch shots 0059–0069, average spectra")
ax[0].legend(fontsize=8)
ax[0].grid(True, which="both", alpha=0.3)
for b in BANDS:
    ax[1].semilogy(1e3 * centres, 1e6 * np.sqrt(band_t[b]), ".-", lw=0.8, label=b)
for sw in switches:
    ax[1].axvline(1e3 * sw, color="gray", ls=":", lw=0.8)
ax[1].set(xlabel="Window centre (ms after record start; probe shutter opened at −7 ms)",
          ylabel="Band ASD (µV/√Hz)", title="Noise along the record (dotted: enable on, enable off, sine end)")
ax[1].legend(fontsize=8)
ax[1].grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig(FIG, dpi=110)
print(f"figure: {FIG}")
