"""Compare launched shots with no-launch (MRR_TRIG_BOOL = False) shots taken at the same settings.

Usage: python no_atom_reference.py <helper dir> <shots dir> <launched range> <no-launch range> [acq]
  e.g. ... 0005-0014 0015-0019 1
Checks: PMT fluorescence with no launch, per-shot noise, and whether subtracting the averaged
no-launch lock-in trace removes the repeatable background from the averaged launched trace.
"""
import sys
import warnings
from pathlib import Path
import numpy as np
from scipy import signal

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot, fit_modulation_reference, analyze_shot

warnings.filterwarnings("ignore")
SHOTS_DIR = Path(sys.argv[2])
ACQ = sys.argv[5] if len(sys.argv) > 5 else "1"
HALF = 0.008


def shot_range(text):
    lo, hi = (int(x) for x in text.split("-"))
    return [f"{n:04d}" for n in range(lo, hi + 1)]


def analyze(number):
    path = next((SHOTS_DIR / number).glob(f"*_{ACQ}.h5"))
    shot = load_shot(path)
    ref = fit_modulation_reference(shot, frequency_guess_hz=shot["mod_freq"] or 6123.0)
    res = analyze_shot(shot, ref, n_periods=2, cloud_half_width_s=HALF,
                       background_guard_s=0.002, switch_guard_s=0.001, harmonics=(1,))
    r = res["results"][("PDA100A2", 1)]
    pmt = res["results"][("PMT", 1)]
    t, v = shot["t_pda"], shot["pda_v"]
    m = (t >= shot["enable_start_s"] + 0.002) & (t <= shot["enable_stop_s"] - 0.002)
    fs, fm = shot["fs_pda"], shot["mod_freq"] or 6123.0
    f, psd = signal.welch(v[m], fs=fs, nperseg=int(round(fs / 97.66)), detrend="linear")
    asd_fm = float(np.median(np.sqrt(psd)[(f >= fm - 600) & (f <= fm - 200)]))
    # Remove each shot's own mean so only the time-dependent pattern is averaged.
    z0 = r["z"] - r["z"][r["valid"]].mean()
    # PMT fluorescence contrast: smoothed peak above the median, in units of the smoothed std.
    ps = res["pmt_smooth"]
    return dict(n=number, t=r["t"], z0=z0, valid=r["valid"], cloud_s=res["cloud_center_s"],
                after_drop=shot["acquisition_after_drop_s"], dc=v[m].mean(), asd=asd_fm,
                snr=r["peak_amp"] / r["noise_rms"], noise=r["noise_rms"],
                pmt_snr=pmt["peak_amp"] / pmt["noise_rms"],
                pmt_counts=int(shot["cumulative_counts"][-1]),
                pmt_peak_rate=float(np.nanmax(ps)), pmt_median_rate=float(np.nanmedian(ps)),
                phase_peak=np.angle(r["z"][r["peak_index"]] - r["background"]))


launched = [analyze(n) for n in shot_range(sys.argv[3])]
nolaunch = [analyze(n) for n in shot_range(sys.argv[4])]

print(f"acquisition _{ACQ}")
print(f"{'shot':>5} {'type':>9} {'PDA V':>6} {'ASD fm':>7} {'noise mV':>9} {'1f SNR':>7} "
      f"{'PMT SNR':>8} {'PMT counts':>11} {'PMT peak/median':>16}")
for kind, group in (("launched", launched), ("no-launch", nolaunch)):
    for s in group:
        print(f"{s['n']:>5} {kind:>9} {s['dc']:6.2f} {s['asd']*1e6:7.2f} {s['noise']*1e3:9.3f} "
              f"{s['snr']:7.2f} {s['pmt_snr']:8.2f} {s['pmt_counts']:11d} "
              f"{s['pmt_peak_rate'] / s['pmt_median_rate']:16.3f}")
for kind, group in (("launched", launched), ("no-launch", nolaunch)):
    print(f"mean {kind:>9}: 1f SNR {np.mean([s['snr'] for s in group]):.2f} ± {np.std([s['snr'] for s in group], ddof=1):.2f}, "
          f"noise {1e3*np.mean([s['noise'] for s in group]):.3f} mV, PMT SNR {np.mean([s['pmt_snr'] for s in group]):.2f}, "
          f"ASD {1e6*np.mean([s['asd'] for s in group]):.2f} uV/rtHz")

# Repeatable background: average each group, compare outside the launched cloud window.
t = launched[0]["t"]
valid = np.logical_and.reduce([s["valid"] for s in launched + nolaunch])
center = np.mean([s["cloud_s"] for s in launched])
outside = valid & (np.abs(t - center) > HALF + 0.002)
avg_l = np.mean([s["z0"] for s in launched], axis=0)
avg_n = np.mean([s["z0"] for s in nolaunch], axis=0)


def crms(z):
    return float(np.sqrt(np.mean(np.abs(z[outside]) ** 2)))


single = np.mean([crms(s["z0"]) for s in launched])
print(f"\nOff-cloud complex RMS: single launched shot {1e3*single:.3f} mV")
print(f"  {len(launched)}-shot launched average {1e3*crms(avg_l):.3f} mV (sqrt(N) expects {1e3*single/np.sqrt(len(launched)):.3f})")
print(f"  {len(nolaunch)}-shot no-launch average {1e3*crms(avg_n):.3f} mV (sqrt(N) expects {1e3*single/np.sqrt(len(nolaunch)):.3f})")
diff = avg_l - avg_n
expect_diff = single * np.sqrt(1 / len(launched) + 1 / len(nolaunch))
print(f"  launched avg - no-launch avg {1e3*crms(diff):.3f} mV (pure noise expects {1e3*expect_diff:.3f})")
corr = np.corrcoef(np.r_[avg_l[outside].real, avg_l[outside].imag],
                   np.r_[avg_n[outside].real, avg_n[outside].imag])[0, 1]
print(f"  correlation of the two averaged off-cloud patterns: {corr:.2f}")

# Signal after subtraction, projected on the launched shots' common 1f phase.
phi0 = np.angle(np.mean([np.exp(1j * s["phase_peak"]) for s in launched]))
cloud = valid & (np.abs(t - center) <= HALF)
for label, trace in (("launched average", avg_l), ("minus no-launch average", diff)):
    x = (trace * np.exp(-1j * phi0)).real
    print(f"  {label:>24}: projected peak {1e3*x[cloud].max():.3f} mV, off-cloud std {1e3*x[outside].std():.4f} mV, "
          f"SNR {x[cloud].max() / x[outside].std():.1f}")
