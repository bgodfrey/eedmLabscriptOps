"""Photodetector blocked vs probe on: which noise is electronic and which needs light?

Usage: python dark_vs_lit.py <helper dir> <dark folder> <lit day dir> <figure.png>

Dark: PDA100A2 blocked, everything else running as usual (exp12, 2026-10-15 0005).
Lit:  no-launch shots 0059-0069 from 2026-10-14 (probe on, no atoms), same gain and sample rate.

Noise that is the same with the detector blocked is electronic (detector amplifier, cable, NI,
pickup from the track or the shutter driver). Noise that disappears when blocked needs light: laser
intensity noise, shot noise, beam motion. For each band the light-induced part is
sqrt(lit^2 - dark^2), since independent noises add in power.

Also checks the modulation carrier: with no light, any signal at f_m on the PDA channel is electrical
crosstalk from the modulation drive, and it would sit inside the lock-in background.
"""
import sys
import warnings
from pathlib import Path
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot, fit_modulation_reference
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DARK, LIT_DAY, FIG = Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
LIT_FOLDERS = [f"{i:04d}" for i in range(59, 70)]
WIN = 2048                                   # 8.2 ms Welch segments at 250 kS/s
GAIN, BW = 1.51e6, 11e3                      # PDA100A2 at 60 dB (Hi-Z), -3 dB bandwidth
BANDS = [(150, 300), (300, 1000), (1000, 3000), (3000, 5500), (6500, 9000), (15e3, 25e3), (40e3, 60e3), (90e3, 120e3)]
TIME_BANDS = {"0.15-1 kHz": (150, 1000), "1-3 kHz": (1000, 3000), "3.17 kHz line": (3100, 3250), "6.5-9 kHz": (6500, 9000)}
GLOBALS = ("t_PMT_delay", "MRR_TRIG_BOOL", "mod_amp", "mod_freq", "mod_offset", "pda_gain_db", "probe_power_uW")


def load_folder(folders):
    shots = []
    for folder in folders:
        for p in sorted(folder.glob("*.h5"), key=lambda q: int(q.stem.split("_")[-1])):
            try:
                s = load_shot(p)
            except Exception as exc:
                print(f"  skip {p.name}: {exc}")
                continue
            with h5py.File(p, "r") as h5:
                s["globals"] = {k: h5["globals"].attrs.get(k) for k in GLOBALS}
            shots.append(s)
    return shots


def stretch_psd(shots, a_s, b_s, key="pda_v"):
    """Average Welch PSD of one time stretch of the record, per-segment linear detrend."""
    out = []
    for s in shots:
        m = (s["t_pda"] >= a_s) & (s["t_pda"] < b_s)
        f, p = signal.welch(s[key][m], fs=s["fs_pda"], nperseg=WIN, noverlap=WIN // 2, detrend="linear")
        out.append(p)
    return f, np.mean(out, axis=0)


def band(f, psd, lo, hi):
    m = (f >= lo) & (f < hi)
    return float(np.sqrt(np.median(psd[m])))       # band-median ASD: the broadband floor, not the lines


def band_vs_time(shots):
    """Band ASD in 8 ms Hann windows stepped by 4 ms along the record, averaged over shots."""
    fs = shots[0]["fs_pda"]
    starts = np.arange(0, int(0.2 * fs) - WIN, WIN // 2)
    hann = signal.windows.hann(WIN)
    norm = fs * np.sum(hann ** 2)
    fw = np.fft.rfftfreq(WIN, 1 / fs)
    acc = {b: np.zeros(len(starts)) for b in TIME_BANDS}
    for s in shots:
        for i, a in enumerate(starts):
            psd = 2 * np.abs(np.fft.rfft(signal.detrend(s["pda_v"][a:a + WIN]) * hann)) ** 2 / norm
            for b, (lo, hi) in TIME_BANDS.items():
                m = (fw >= lo) & (fw < hi)
                acc[b][i] += (psd[m].mean() if "line" in b else np.median(psd[m])) / len(shots)
    return (starts + WIN / 2) / fs, {b: np.sqrt(v) for b, v in acc.items()}


def carrier(shots):
    """Lock-in amplitude at f_m (1f and 2f) averaged over the enabled window, as the analysis computes it."""
    rows = []
    for s in shots:
        ref = fit_modulation_reference(s, s["mod_freq"] or 6125.0)
        m = (s["t_pda"] >= s["enable_start_s"] + 0.005) & (s["t_pda"] < s["enable_stop_s"] - 0.005)
        row = []
        for h in (1, 2):
            t, z = cse.demodulate(s, ref, harmonic=h)
            row.append(np.mean(z[m]))
        rows.append(row)
    return np.array(rows)


dark = load_folder([DARK])
lit = load_folder([LIT_DAY / f for f in LIT_FOLDERS])
print(f"dark shots: {len(dark)}   lit shots: {len(lit)}")
print("dark-run globals (first shot):", {k: (round(float(v), 4) if v is not None else None)
                                        for k, v in dark[0]["globals"].items()})
varying = {k for k in GLOBALS if len({str(s['globals'][k]) for s in dark}) > 1}
if varying:
    print("  globals that vary across the dark run:", varying)
dc_d = np.array([s["pda_v"].mean() for s in dark])
dc_l = np.array([s["pda_v"].mean() for s in lit])
print(f"PDA DC level: dark {1e3*dc_d.mean():.2f} mV (shot-to-shot spread {1e3*dc_d.std():.2f} mV), "
      f"lit {dc_l.mean():.2f} V")
passes = {("return" if s["pmt_delay_s"] > 0.2 else "upward") for s in dark}
print(f"dark passes recorded: {sorted(passes)}")

sine0 = np.median([s["sine_start_s"] for s in dark])
stretches = {"before modulation": (0.010, sine0 - 0.002),
             "enabled window": (sine0 + 0.003, sine0 + 0.067),
             "after the sine": (sine0 + 0.103, 0.198)}
spec = {}
for name, (a, b) in stretches.items():
    spec[("dark", name)] = stretch_psd(dark, a, b)
    spec[("lit", name)] = stretch_psd(lit, a, b)

f = spec[("lit", "after the sine")][0]
shot_asd = np.sqrt(2 * 1.602e-19 * GAIN * dc_l.mean() / (1 + (f / BW) ** 2))
print("\nBand-median ASD (uV/rtHz), steady stretch after the sine ends:")
print(f"  {'band (Hz)':>14} {'lit':>7} {'dark':>7} {'light-induced':>14} {'shot (lit DC)':>14} {'light-ind./shot':>15}")
for lo, hi in BANDS:
    l, d = band(f, spec[("lit", "after the sine")][1], lo, hi), band(f, spec[("dark", "after the sine")][1], lo, hi)
    li = np.sqrt(max(l ** 2 - d ** 2, 0))
    sh = band(f, shot_asd ** 2, lo, hi)
    print(f"  {lo:>6.0f}-{hi:<7.0f} {1e6*l:7.2f} {1e6*d:7.2f} {1e6*li:14.2f} {1e6*sh:14.2f} {li/sh:15.2f}")

print("\nSame, stretch before the modulation (shutter transient region, 10 ms to sine start):")
for lo, hi in BANDS[:5]:
    l, d = band(f, spec[("lit", "before modulation")][1], lo, hi), band(f, spec[("dark", "before modulation")][1], lo, hi)
    print(f"  {lo:>6.0f}-{hi:<7.0f} lit {1e6*l:7.2f}   dark {1e6*d:7.2f}")

# Noise next to f_m, the number the lock-in sees (excluding the carrier itself).
for label, key in (("lit", "lit"), ("dark", "dark")):
    fe, pe = spec[(key, "enabled window")]
    m = ((fe > 5300) & (fe < 5900)) | ((fe > 6350) & (fe < 6950))
    print(f"\nnext to f_m in the enabled window, {label}: {1e6*np.sqrt(np.median(pe[m])):.2f} uV/rtHz", end="")
print()

# Narrow lines: present when dark -> electrical.
for label, key in (("dark", "dark"), ("lit", "lit")):
    fa, pa = stretch_psd(dark if key == "dark" else lit, 0.010, 0.198)
    ratio = pa / signal.medfilt(pa, 31)
    pk, _ = signal.find_peaks(ratio, height=4.0)
    top = pk[np.argsort(ratio[pk])[::-1][:14]]
    print(f"lines, whole record, {label:4s}: " + ", ".join(f"{fa[k]/1e3:.2f}k (x{np.sqrt(ratio[k]):.1f}, "
          f"{1e6*np.sqrt(pa[k]):.0f})" for k in sorted(top)))

# Electrical crosstalk of the modulation into the PDA channel.
cz_d, cz_l = carrier(dark), carrier(lit)
print("\nLock-in carrier over the enabled window (complex mean, |Z| in uV):")
for h in (0, 1):
    print(f"  {h+1}f: dark |Z| = {1e6*np.abs(cz_d[:, h]).mean():8.1f} (phase {np.degrees(np.angle(cz_d[:, h].mean())):6.1f} deg), "
          f"lit |Z| = {1e6*np.abs(cz_l[:, h]).mean():8.1f} (phase {np.degrees(np.angle(cz_l[:, h].mean())):6.1f} deg)")

# Coherence with the reference channel ai1 while dark.
fc, coh = [], []
for s in dark:
    m = (s["t_pda"] >= sine0 + 0.103) & (s["t_pda"] < 0.198)
    fc, c = signal.coherence(s["pda_v"][m], s["reference_v"][m], fs=s["fs_pda"], nperseg=WIN)
    coh.append(c)
coh = np.mean(coh, axis=0)
print("\ndark coherence ai0-ai1 (after the sine), band medians: " + ", ".join(
    f"{lo:.0f}-{hi:.0f} Hz {np.median(coh[(fc >= lo) & (fc < hi)]):.2f}" for lo, hi in BANDS))

# ---- figure ------------------------------------------------------------------------------------
LIT_C, DARK_C, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, ax = plt.subplots(2, 1, figsize=(10, 9), facecolor="#fcfcfb")
for a in ax:
    a.set_facecolor("#fcfcfb")
    a.grid(True, which="major", color="#e4e3df", lw=0.6)
for key, col, name in (("lit", LIT_C, f"Probe on, no atoms (0059–0069, {dc_l.mean():.1f} V)"),
                       ("dark", DARK_C, "Photodetector blocked (2026-10-15 0005)")):
    fs_, ps_ = spec[(key, "after the sine")]
    ax[0].loglog(fs_[1:], 1e6 * np.sqrt(ps_[1:]), color=col, lw=1.2, label=name)
ax[0].loglog(f[1:], 1e6 * shot_asd[1:], color=MUTED, lw=1.0, ls="--", label="Shot noise at the lit DC level")
ax[0].axvline(6125, color=MUTED, lw=0.8, ls=":")
ax[0].text(6125 * 1.04, ax[0].get_ylim()[1] * 0.6, "f_m", color=MUTED)
ax[0].set(xlabel="Frequency (Hz)", ylabel="PDA noise (µV/√Hz)",
          title="Steady noise, 146–198 ms into the record (after the modulation ends)")
ax[0].legend(frameon=False, loc="upper right")
for key, col, ls in (("lit", LIT_C, "-"), ("dark", DARK_C, "-")):
    tc, bands_t = band_vs_time(lit if key == "lit" else dark)
    for (b, v), mk in zip(bands_t.items(), ("o", "s", "^", "D")):
        ax[1].semilogy(1e3 * tc, 1e6 * v, color=col, lw=1.2, marker=mk, ms=4, ls=ls,
                       label=f"{'Probe on' if key == 'lit' else 'Blocked'}: {b}")
for sw in (sine0, sine0 + 0.07, sine0 + 0.1):
    ax[1].axvline(1e3 * sw, color=MUTED, lw=0.8, ls=":")
ax[1].set(xlabel="Time in the record (ms); probe shutter commanded open at −7 ms",
          ylabel="Band noise (µV/√Hz)", title="Noise along the record (dotted: enable on, enable off, sine end)")
ax[1].legend(frameon=False, fontsize=7.5, ncol=2, loc="upper right")
fig.tight_layout()
fig.savefig(FIG, dpi=110, facecolor=fig.get_facecolor())
print(f"\nfigure: {FIG}")
