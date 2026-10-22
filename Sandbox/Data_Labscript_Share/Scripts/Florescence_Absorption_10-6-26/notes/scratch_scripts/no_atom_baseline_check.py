"""Does the early 'tail' in the vertical dip profile appear in shots with no atoms?

Usage: python no_atom_baseline_check.py <helper dir> <dark mV>
10/17 0000/0001 include shots with a full lock loss (MOT ~ 0, no cloud). Those traces go through the
same baseline procedure as the 0 mm atom shots (0003/0008/0013), aligned to the median cloud time
of the atom shots (the sequence timing is fixed). If the no-atom 'profile' shows the +z tail, it is a
probe-light transient; if it is flat, the tail is atoms.
"""
import sys
import warnings
from pathlib import Path
import numpy as np

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot, load_mot_trace, fit_mot_loading
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DARK = float(sys.argv[2]) * 1e-3
DAY = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\17")
V = 2.30
edges = np.arange(-0.022, 0.012 + 1e-9, 1e-3)          # 1 ms bins, -22..+12 ms
cen = 0.5 * (edges[1:] + edges[:-1])


def trace(s):
    k = int(round(s["fs_pda"] / 10e3))
    m = len(s["pda_v"]) // k * k
    tb = s["t_pda"][:m].reshape(-1, k).mean(axis=1)
    vb = s["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
    w = (tb > s["enable_start_s"] + 0.002) & (tb < s["enable_stop_s"] - 0.002)
    return tb[w], vb[w]


def profile(tb, vb, tc, lo):
    u = tb - tc
    far = np.abs(u) > lo
    fr = 1 - vb / np.polyval(np.polyfit(u[far], vb[far], 2), u)
    idx = np.digitize(u, edges) - 1
    return np.array([fr[idx == i].mean() for i in range(len(cen))])


atom, tcs = [], []
for f in ("0003", "0008", "0013"):
    for p in sorted((DAY / f).glob("*.h5")):
        s = load_shot(p)
        pm = cse.pmt_cloud_fit(s, 6125.0)
        tcs.append(pm["center_s"])
        atom.append(trace(s) + (pm["center_s"],))
tc_med = np.median(tcs)
print(f"atom shots: t_c median {1e3*tc_med:.2f} ms, spread {1e3*np.std(tcs):.2f} ms")
empty = []
for f in ("0000", "0001"):
    for p in sorted((DAY / f).glob("*.h5")):
        s = load_shot(p)
        try:
            mot = fit_mot_loading(load_mot_trace(p))["atoms_v_at_coils_off"]
        except Exception:
            mot = np.nan
        pm = cse.pmt_cloud_fit(s, 6125.0)
        tb, vb = trace(s)
        amp = pm["amp_cps"] if pm["ok"] else 0.0
        print(f"  {f}/{p.name}: MOT {1e3*mot:6.1f} mV, PMT amp {amp/1e6:5.2f} Mcps, light {np.median(vb):.3f} V")
        if (np.isfinite(mot) and mot < 0.005) or amp < 0.3e6:
            empty.append((tb, vb))
print(f"{len(empty)} no-atom shots")
for lo in (0.006, 0.010):
    a = np.mean([profile(tb, vb, tc, lo) for tb, vb, tc in atom], axis=0)
    e = np.array([profile(tb, vb, tc_med, lo) for tb, vb in empty])
    em, es = e.mean(axis=0), e.std(axis=0, ddof=1) / np.sqrt(len(e))
    print(f"\nbaseline |u| > {1e3*lo:.0f} ms        u ms:  " + " ".join(f"{1e3*c:+6.1f}" for c in cen[::2]))
    print("  z mm:                         " + " ".join(f"{-1e3*V*c:+6.1f}" for c in cen[::2]))
    print("  atoms (dip %):                " + " ".join(f"{100*x:+6.2f}" for x in a[::2]))
    print("  no atoms (dip %):             " + " ".join(f"{100*x:+6.2f}" for x in em[::2]))
    print("  no atoms sem:                 " + " ".join(f"{100*x:6.2f}" for x in es[::2]))
