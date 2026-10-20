"""Fixed-phase projection and coherent multi-shot averaging of the 1f PDA lock-in."""
import sys
import warnings
from pathlib import Path
import numpy as np

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot, fit_modulation_reference, analyze_shot

warnings.filterwarnings("ignore")
HALF = 0.008

shots = []
for p in map(Path, sys.argv[2:]):
    shot = load_shot(p)
    ref = fit_modulation_reference(shot, frequency_guess_hz=shot["mod_freq"] or 6123.0)
    res = analyze_shot(shot, ref, n_periods=2, cloud_half_width_s=HALF,
                       background_guard_s=0.002, switch_guard_s=0.001, harmonics=(1,))
    r = res["results"][("PDA100A2", 1)]
    w = r["z"] - r["background"]
    cloud = np.abs(r["t"] - res["cloud_center_s"]) <= HALF
    shots.append(dict(name=p.stem.split("_")[1], w=w, valid=r["valid"], bg=r["background_mask"] & r["valid"],
                      cloud=cloud & r["valid"], phase=np.angle(w[r["peak_index"]]),
                      r_snr=r["peak_amp"] / r["noise_rms"]))

# Common signal phase from all shots (circular mean of the peak phases).
phi0 = np.angle(np.mean([np.exp(1j * s["phase"]) for s in shots]))
print(f"Common 1f signal phase: {np.rad2deg(phi0):.1f} deg")

print(f"{'shot':>5} {'R SNR':>7} {'proj SNR':>9}")
for s in shots:
    x = (s["w"] * np.exp(-1j * phi0)).real
    noise = x[s["bg"]].std()
    peak = x[s["cloud"]].max()
    s["x"] = x
    print(f"{s['name']:>5} {s['r_snr']:7.2f} {peak / noise:9.2f}")
r_mean = np.mean([s["r_snr"] for s in shots])
p_mean = np.mean([s["x"][s["cloud"]].max() / s["x"][s["bg"]].std() for s in shots])
print(f"mean  {r_mean:7.2f} {p_mean:9.2f}   (projection gain {p_mean / r_mean:.2f}x)")

# Coherent average: identical sample grids, so average the projected traces directly.
valid = np.logical_and.reduce([s["valid"] for s in shots])
bg = np.logical_and.reduce([s["bg"] for s in shots])
cloud = np.logical_or.reduce([s["cloud"] for s in shots]) & valid
xavg = np.mean([s["x"] for s in shots], axis=0)
noise_avg = xavg[bg].std()
peak_avg = xavg[cloud].max()
print(f"\n{len(shots)}-shot coherent average, projected: peak {1e3*peak_avg:.3f} mV, noise {1e3*noise_avg:.4f} mV, "
      f"SNR {peak_avg / noise_avg:.1f}  (single-shot projected mean {p_mean:.1f}; sqrt(N) would give "
      f"{p_mean * np.sqrt(len(shots)):.1f})")
