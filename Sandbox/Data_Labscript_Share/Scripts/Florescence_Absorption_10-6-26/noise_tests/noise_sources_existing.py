"""Where does the PDA noise come from? First pass using shots already on disk (no new data).

Usage: python noise_sources_existing.py <helper dir> <day dir>

Each 0.2 s PDA record has a stretch with the probe on but no modulation and no atoms: from the start of
the acquisition until the sine starts (~43 ms in). That "pre" stretch is a clean noise measurement.
This script asks four questions of it:

  1. Spectrum shape: how steep is the low-frequency (1/f-like) part, and where does it meet the white
     floor relative to f_m = 6.1 kHz? That sets the range worth scanning for the modulation frequency.
  2. Scaling with light level: the probe power drifts shot to shot (thermal swing), so the DC level
     differs between shots. For each band, fit  ASD ~ DC^s  across shots:
         s ~ 1   noise proportional to the light   -> optical (laser intensity, pointing, etalons, polarization)
         s ~ 0.5 shot noise
         s ~ 0   independent of the light           -> electronic (detector, cabling, ADC, pickup)
  3. Other channels: the reference channel ai1 carries no light. Lines or 1/f that also appear there,
     coherent with ai0, are electrical pickup. (ai0 leaks into ai1 at about -3 mV/V; the script
     subtracts that expectation before calling anything pickup.) The MOT PDA100A on ai2 is dark for
     the first ~8 ms of its record: a second, independent detector in the same electrical environment.
  4. Track motion: launched shots (stage moves) vs no-launch shots (MRR_TRIG_BOOL = 0, stage still),
     same settings, same session. Extra noise in launched shots comes from the stage moving.
"""
import sys
import warnings
from pathlib import Path
import h5py
import numpy as np
from scipy import signal

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot

warnings.filterwarnings("ignore")
DAY = Path(sys.argv[2])
LAUNCHED = [f"{i:04d}" for i in range(49, 59)]          # 0049-0058: launched, mod_offset 6 mV, mod_amp 10 mV
NO_LAUNCH = [f"{i:04d}" for i in range(59, 70)]         # 0059-0068 + 0069: same settings, stage not triggered
GAIN_V_PER_A = 1.51e6     # PDA100A2 at 60 dB, high-impedance load (the NI input)
PDA_BW_HZ = 11e3          # PDA100A2 -3 dB bandwidth at 60 dB
NPERSEG = 2048            # 8.2 ms Welch segments at 250 kS/s -> 122 Hz resolution
CROSSTALK = -3e-3         # measured ai0 -> ai1 leak (V per V)
BANDS = [(150, 300), (300, 1000), (1000, 3000), (3000, 5500), (6500, 9000), (15e3, 25e3), (40e3, 60e3), (90e3, 120e3)]


def welch(x, fs, detrend="linear"):
    """One-sided PSD (V^2/Hz). Each 8 ms segment gets its own straight-line detrend, so a slow drift
    across the record does not leak into the low-frequency bins through the window sidelobes."""
    return signal.welch(x, fs=fs, nperseg=min(NPERSEG, len(x)), noverlap=min(NPERSEG, len(x)) // 2,
                        detrend=detrend)


def band_median(f, psd, lo, hi):
    """Median PSD in a band: robust against narrow lines, so it reports the broadband floor."""
    m = (f >= lo) & (f < hi)
    return float(np.median(psd[m]))


def collect(folders):
    """Pre-stretch PSDs of ai0 and ai1, their cross-spectrum, the DC level, and the dark ai2 PSD."""
    out = []
    for folder in folders:
        for p in sorted((DAY / folder).glob("*.h5")):
            try:
                shot = load_shot(p)
            except Exception as exc:             # skip anything the validated loader rejects
                print(f"  skip {p.name}: {exc}")
                continue
            fs, t = shot["fs_pda"], shot["t_pda"]
            # Probe shutter is commanded open 7 ms before the record starts; skip 10 ms more for it to settle.
            pre = (t >= 0.010) & (t < shot["sine_start_s"] - 0.002)
            x, r = shot["pda_v"][pre], shot["reference_v"][pre]
            f, pxx = welch(x, fs)
            _, prr = welch(r, fs)
            _, pxr = signal.csd(x, r, fs=fs, nperseg=NPERSEG, noverlap=NPERSEG // 2, detrend="linear")
            _, pxx_const = welch(x, fs, detrend="constant")
            # ai2: MOT PDA100A (70 dB) record starts 10 ms before the MOT shutter opens -> dark at first.
            with h5py.File(p, "r") as h5:
                mot = h5["data/traces/MOT_fluorescence"][()] if "MOT_fluorescence" in h5["data/traces"] else None
            p_dark = None
            if mot is not None:
                dark = np.asarray(mot["values"], float)[: int(0.008 * fs)]
                _, p_dark = welch(dark, fs)
            out.append({"name": p.name, "acq": "return" if shot["pmt_delay_s"] > 0.2 else "upward",
                        "dc": float(x.mean()), "f": f, "pxx": pxx, "pxx_const": pxx_const, "prr": prr,
                        "pxr": pxr, "p_dark": p_dark, "fs": fs})
    return out


def summarize(label, shots):
    """Average PSD over shots; print band ASDs in uV/rtHz and relative to DC."""
    f = shots[0]["f"]
    pxx = np.mean([s["pxx"] for s in shots], axis=0)
    dc = np.mean([s["dc"] for s in shots])
    cells = [f"{1e6*np.sqrt(band_median(f, pxx, lo, hi)):6.2f}" for lo, hi in BANDS]
    print(f"  {label:26s} n={len(shots):2d} DC {dc:4.1f} V | " + " ".join(cells))
    return f, pxx, dc


launched, no_launch = collect(LAUNCHED), collect(NO_LAUNCH)
f = no_launch[0]["f"]
hdr = " ".join(f"{(lo/1e3 if lo >= 1e3 else lo):>5g}{'k' if lo >= 1e3 else ''}".rjust(6) for lo, _ in BANDS)
print("\n1) Probe on, no modulation, no atoms ('pre' stretch). Band-median ASD in uV/rtHz; band lower edges (Hz):")
print(f"  {'':26s} {'':16s}   {hdr}")
groups = {}
for acq in ("upward", "return"):
    for label, shots in (("launched", launched), ("no launch", no_launch)):
        sel = [s for s in shots if s["acq"] == acq]
        if sel:
            groups[(acq, label)] = summarize(f"{acq}, {label}", sel)
all_nl = no_launch
_, pxx_nl, dc_nl = summarize("all no-launch", all_nl)
shot_asd = np.sqrt(2 * 1.602e-19 * GAIN_V_PER_A * dc_nl / (1 + (f / PDA_BW_HZ) ** 2))
print(f"  {'shot noise at that DC':26s} {'':16s} | " + " ".join(
    f"{1e6*np.sqrt(band_median(f, shot_asd**2, lo, hi)):6.2f}" for lo, hi in BANDS))
pxx_c = np.mean([s["pxx_const"] for s in all_nl], axis=0)
print(f"  {'no-launch, mean removed only':26s} {'':16s} | " + " ".join(
    f"{1e6*np.sqrt(band_median(f, pxx_c, lo, hi)):6.2f}" for lo, hi in BANDS)
      + "   <- without per-segment detrend (drift leakage)")

# Spectral slope of the excess over shot noise, 150 Hz - 3 kHz.
excess = pxx_nl - shot_asd ** 2
m = (f >= 150) & (f <= 3000) & (excess > 0)
alpha = -np.polyfit(np.log(f[m]), np.log(excess[m]), 1)[0]
print(f"\n   excess PSD over shot noise falls as 1/f^{alpha:.2f} between 150 Hz and 3 kHz")
for target in (1000, 3000, 6125, 9000, 12000):
    k = np.argmin(abs(f - target))
    lo = max(k - 3, 0)
    tot = np.median(pxx_nl[lo:k + 4])
    print(f"   at {target/1e3:5.2f} kHz: total {1e6*np.sqrt(tot):5.2f}, shot {1e6*shot_asd[k]:4.2f} uV/rtHz "
          f"-> excess/shot (power) = {(tot - shot_asd[k]**2)/shot_asd[k]**2:5.2f}")

print("\n2) Scaling with light level across no-launch shots: ASD ~ DC^s (1 optical, 0.5 shot, 0 electronic)")
dcs = np.array([s["dc"] for s in all_nl])
print(f"   DC range {dcs.min():.2f}-{dcs.max():.2f} V over {len(dcs)} shots")
for lo, hi in BANDS:
    a = np.array([np.sqrt(band_median(s["f"], s["pxx"], lo, hi)) for s in all_nl])
    s_fit, cov = np.polyfit(np.log(dcs), np.log(a), 1, cov=True)
    print(f"   {lo:>7.0f}-{hi:<7.0f} Hz: s = {s_fit[0]:+.2f} +/- {np.sqrt(cov[0,0]):.2f}")

print("\n3) Electrical check: reference channel ai1 (no light) and its coherence with ai0, no-launch pre stretch")
prr = np.mean([s["prr"] for s in all_nl], axis=0)
pxr = np.mean([s["pxr"] for s in all_nl], axis=0)
coh = np.abs(pxr) ** 2 / (pxx_nl * prr)
for lo, hi in BANDS:
    mm = (f >= lo) & (f < hi)
    leak = abs(CROSSTALK) * np.sqrt(band_median(f, pxx_nl, lo, hi))
    print(f"   {lo:>7.0f}-{hi:<7.0f} Hz: ai1 {1e6*np.sqrt(band_median(f, prr, lo, hi)):6.2f} uV/rtHz "
          f"(crosstalk would give {1e6*leak:5.3f}), coherence with ai0 {np.median(coh[mm]):.2f}")
darks = [s["p_dark"] for s in all_nl if s["p_dark"] is not None]
if darks:
    fd = signal.welch(np.zeros(int(0.008 * all_nl[0]["fs"])), fs=all_nl[0]["fs"], nperseg=NPERSEG)[0]
    pd_ = np.mean(darks, axis=0)
    print("   ai2 MOT PDA100A, dark (first 8 ms, 70 dB): " + ", ".join(
        f"{lo:.0f}-{hi:.0f} Hz {1e6*np.sqrt(band_median(fd, pd_, lo, hi)):.2f}" for lo, hi in BANDS[:5]) + " uV/rtHz")

print("\n4) Narrow lines in the no-launch pre stretch (ratio to the local median, top 12, ai0 and ai1):")
for name, spec in (("ai0 PDA", pxx_nl), ("ai1 reference", prr)):
    base = signal.medfilt(spec, 31)
    ratio = spec / base
    peaks, _ = signal.find_peaks(ratio, height=4.0)
    top = peaks[np.argsort(ratio[peaks])[::-1][:12]]
    print(f"   {name:14s}: " + ", ".join(f"{f[k]/1e3:.2f} kHz (x{np.sqrt(ratio[k]):.1f})" for k in sorted(top)))
