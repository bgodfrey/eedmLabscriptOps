import sys, warnings
from pathlib import Path
import numpy as np
sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot, load_mot_trace, fit_mot_loading
import cloud_signal_estimation as cse
warnings.filterwarnings("ignore")
root = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10")
for day, folder in (("14", "0130"), ("16", "0023"), ("16", "0028"), ("16", "0033")):
    F, M = [], []
    for p in sorted((root / day / folder).glob("*.h5")):
        s = load_shot(p)
        pm = cse.pmt_cloud_fit(s, s["mod_freq"] or 6125.0)
        if not pm["ok"]:
            continue
        F.append(pm["amp_cps"] * pm["sigma_s"])
        try:
            M.append(fit_mot_loading(load_mot_trace(p))["atoms_v_at_coils_off"])
        except Exception:
            M.append(np.nan)
    F, M = np.array(F), np.array(M)
    print(f"10/{day} {folder}: {len(F)} shots, mod_offset {1e3*s['mod_offset']:.1f} mV, PMT integral {F.mean():.0f} counts, "
          f"sd {100*F.std()/F.mean():.0f}%, PMT/MOT sd {100*np.nanstd(F/M)/np.nanmean(F/M):.0f}%, MOT {1e3*np.nanmean(M):.1f} mV")
