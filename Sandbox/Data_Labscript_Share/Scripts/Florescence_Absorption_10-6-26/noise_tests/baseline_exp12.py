"""Per-shot lock-in and MOT metrics for exp12 shots, plus mean/std per acquisition (_0 upward, _1 returning).

Usage: python baseline_exp12.py <helper dir> <shot files...>
SNR uses a background more than 20 ms from the cloud (cloud ±8 ms + 12 ms guard), where the
cloud's own wings no longer leak in; SNR1_near repeats the original ±8 ms + 2 ms guard for
comparison with older log entries. MOT columns are empty for shots without the ai2 trace.
"""
import sys
import warnings
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import signal

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import (load_shot, fit_modulation_reference, analyze_shot,
                                            load_mot_trace, fit_mot_loading)

warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
COMMON = dict(n_periods=2, cloud_half_width_s=0.008, switch_guard_s=0.001, harmonics=(1, 2))

rows = []
for p in map(Path, sys.argv[2:]):
    try:
        shot = load_shot(p)
        fm = shot["mod_freq"] or 6123.0
        ref = fit_modulation_reference(shot, frequency_guess_hz=fm)
        res = analyze_shot(shot, ref, background_guard_s=0.012, **COMMON)
        near_res = analyze_shot(shot, ref, background_guard_s=0.002, **COMMON)
    except Exception as e:
        print(f"{p.name}: {type(e).__name__}: {e}")
        continue
    with h5py.File(p, "r") as h5:
        launch = h5["globals"].attrs.get("MRR_TRIG_BOOL")
    mot = load_mot_trace(p)
    mot_fit = fit_mot_loading(mot) if mot is not None else None
    fs = shot["fs_pda"]
    t, v = shot["t_pda"], shot["pda_v"]
    m = (t >= shot["enable_start_s"] + 0.002) & (t <= shot["enable_stop_s"] - 0.002)
    f, psd = signal.welch(v[m], fs=fs, nperseg=int(round(fs / 97.66)), detrend="linear")
    asd = np.sqrt(psd)
    near = (f >= fm - 600) & (f <= fm - 200)
    floor = (f >= 60e3) & (f <= min(200e3, 0.8 * fs / 2))
    p1, p2 = res["results"][("PDA100A2", 1)], res["results"][("PDA100A2", 2)]
    pmt = res["results"][("PMT", 1)]
    # Signed phase at the 1f peak after background subtraction: tells which side of the line.
    z_peak = p1["z"][p1["peak_index"]]
    rows.append({
        "shot": p.stem.split("_")[1], "acq": p.stem[-1],
        "launch": None if launch is None else int(launch),
        "fs_kHz": fs / 1e3,
        "offset_mV": 1e3 * shot["mod_offset"] if shot["mod_offset"] is not None else np.nan,
        "amp_mV": 1e3 * shot["mod_amp"] if shot["mod_amp"] is not None else np.nan,
        "gain_dB": shot["pda_gain_db"], "P_uW": shot["probe_power_uW"],
        "PDA_V": v[m].mean(),
        "ASD_fm_uV": 1e6 * np.median(asd[near]),
        "floor_uV": 1e6 * np.median(asd[floor]) if floor.any() else np.nan,
        "RIN_fm_dBc": 20 * np.log10(np.median(asd[near]) / v[m].mean()),
        "bg1_mV": 1e3 * p1["noise_rms"], "pk1_mV": 1e3 * p1["peak_amp"],
        "SNR1": p1["peak_amp"] / p1["noise_rms"],
        "SNR1_near": near_res["results"][("PDA100A2", 1)]["peak_amp"]
                     / near_res["results"][("PDA100A2", 1)]["noise_rms"],
        "frac1_ppm": 1e6 * p1["peak_amp"] / v[m].mean(),
        "phase1_deg": np.rad2deg(np.angle(z_peak)),
        "SNR2": p2["peak_amp"] / p2["noise_rms"],
        "PMT_SNR": pmt["peak_amp"] / pmt["noise_rms"],
        "MOT_atoms_mV": 1e3 * mot_fit["atoms_v_at_coils_off"] if mot_fit else np.nan,
        "MOT_tau_s": mot_fit["tau_s"] if mot_fit else np.nan,
        "MOT_scatter_mV": 1e3 * mot_fit["scatter_v"] if mot_fit else np.nan,
        "PMT_peak_t_ms": 1e3 * (res["pmt_peak_s"] + shot["acquisition_after_drop_s"]),
        "PDA_peak_t_ms": 1e3 * (p1["peak_time_s"] + shot["acquisition_after_drop_s"]),
    })

df = pd.DataFrame(rows)
# Peak per unit MOT atom signal removes shot-to-shot atom-number swings from the comparison.
df["pk1_per_MOT"] = df["pk1_mV"] / df["MOT_atoms_mV"]
print(df.round(2).to_string(index=False))
num = ["PDA_V", "ASD_fm_uV", "floor_uV", "RIN_fm_dBc", "bg1_mV", "pk1_mV", "SNR1", "SNR1_near",
       "frac1_ppm", "SNR2", "PMT_SNR", "MOT_atoms_mV", "MOT_tau_s", "MOT_scatter_mV", "pk1_per_MOT"]
for acq, g in df.groupby("acq"):
    print(f"\n--- acquisition _{acq}: {len(g)} shots ---")
    summary = pd.DataFrame({"mean": g[num].mean(), "std": g[num].std(),
                            "cv_%": 100 * g[num].std() / g[num].mean().abs()})
    print(summary.round(3).to_string())
    for col in ("pk1_mV", "SNR1", "PMT_SNR"):
        if g["MOT_atoms_mV"].notna().sum() > 2:
            print(f"corr({col}, MOT_atoms_mV) = {g[col].corr(g['MOT_atoms_mV']):.2f}")
    print(f"corr(SNR1, PMT_SNR) = {g['SNR1'].corr(g['PMT_SNR']):.2f}")
