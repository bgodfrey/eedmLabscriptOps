import sys
import warnings
from pathlib import Path
import numpy as np
from scipy import signal

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot, fit_modulation_reference, analyze_shot

warnings.filterwarnings("ignore")
FM = 6123.0
rows = []
for p in sys.argv[2:]:
    p = Path(p)
    try:
        shot = load_shot(p)
        ref = fit_modulation_reference(shot, frequency_guess_hz=FM)
        res = analyze_shot(shot, ref, n_periods=2, cloud_half_width_s=0.008,
                           background_guard_s=0.002, switch_guard_s=0.001, harmonics=(1, 2))
    except Exception as e:
        print(f"{p.name}: {type(e).__name__}: {e}")
        continue
    fs = shot["fs_pda"]
    t, v = shot["t_pda"], shot["pda_v"]
    m = (t >= shot["enable_start_s"] + 0.002) & (t <= shot["enable_stop_s"] - 0.002)
    # Fixed ~98 Hz resolution regardless of sample rate, so ASDs compare directly.
    nper = int(round(fs / 97.66))
    f, psd = signal.welch(v[m], fs=fs, nperseg=nper, detrend="linear")
    asd = np.sqrt(psd)
    band = lambda lo, hi: 1e6 * float(np.median(asd[(f >= lo) & (f <= hi)]))
    pda1 = res["results"][("PDA100A2", 1)]
    pmt1 = res["results"][("PMT", 1)]
    print(f"{p.name}: fs {fs/1e3:.0f}k | DC {v[m].mean():.2f} V | ASD uV/rtHz: "
          f"5.5-5.9k {band(5500, 5900):.2f}, 6.4-7k {band(6400, 7000):.2f}, "
          f"20-45k {band(20e3, 45e3):.2f}, 60-200k {band(60e3, 200e3) if fs > 150e3 else float('nan'):.2f} | "
          f"lockin 1f: bg RMS {1e3*pda1['noise_rms']:.3f} mV, peak {1e3*pda1['peak_amp']:.3f} mV, "
          f"ratio {pda1['peak_amp']/pda1['noise_rms']:.2f} | PMT ratio {pmt1['peak_amp']/pmt1['noise_rms']:.2f}, "
          f"ref amp {1e3*ref['amp_v']:.2f} mV resid {1e3*ref['residual_rms_v']:.3f} mV")
