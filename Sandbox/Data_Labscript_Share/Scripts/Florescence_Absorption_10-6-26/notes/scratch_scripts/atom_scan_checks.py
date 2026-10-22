"""Checks on the 2026-10-16 iris atom scan (0023-0035).

Usage: python atom_scan_checks.py <helper dir> <dark mV> <figure.png>

1. Per-shot scatter: does the dip follow the PMT fluorescence (detuning and atom number act on both)
   or the MOT signal (atom number only)? Ratio D / (PMT integral / light) cancels both.
2. Shape of the dip at each position: mean fractional deviation from a quadratic baseline fitted
   outside +-6 ms of the cloud, aligned on the PMT centre. A real absorption dip has the PMT-like
   Gaussian shape; a mechanical (launch) transient at a beam edge would not.
3. Same fixed-shape dip fit on the 2026-10-15 line-centre shots (0007, 0008: mod_offset 0.0095,
   iris at 0 mm) with their MOT signal, for comparison with today's 0 mm (mod_offset 0).
"""
import sys
import warnings
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import optimize

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot, load_mot_trace, fit_mot_loading
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DARK, FIG = float(sys.argv[2]) * 1e-3, Path(sys.argv[3])
ROOT = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10")
OFFSET, WIDTH = -0.13e-3, 1.48e-3          # common dip shape from iris_atom_scan.py


def analyse(p):
    shot = load_shot(p)
    with h5py.File(p, "r") as h5:
        g = h5["globals"].attrs
        x, off = float(g["iris_x_mm"]), float(g["mod_offset"])
    pmt = cse.pmt_cloud_fit(shot, shot["mod_freq"] or 6125.0)
    k = int(round(shot["fs_pda"] / 10e3))
    m = len(shot["pda_v"]) // k * k
    tb = shot["t_pda"][:m].reshape(-1, k).mean(axis=1)
    vb = shot["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
    win = (tb > shot["enable_start_s"] + 0.002) & (tb < shot["enable_stop_s"] - 0.002)
    tw, vw, tc = tb[win], vb[win], pmt["center_s"]
    centre = tc + OFFSET

    def model(t, c0, c1, c2, d):
        u = t - tc
        return (c0 + c1 * u + c2 * u ** 2) * (1 - d * np.exp(-0.5 * ((t - centre) / WIDTH) ** 2))

    pf, _ = optimize.curve_fit(model, tw, vw, p0=[np.median(vw), 0, 0, 0.01])
    # baseline outside the cloud, for the shape check
    out = np.abs(tw - tc) > 0.006
    cb = np.polyfit(tw[out] - tc, vw[out], 2)
    frac = vw / np.polyval(cb, tw - tc) - 1
    try:
        mot = fit_mot_loading(load_mot_trace(p))["atoms_v_at_coils_off"]
    except Exception:
        mot = np.nan
    return {"x": x, "off": off, "D": pf[3], "light": pf[0], "pmt": pmt["amp_cps"] * pmt["sigma_s"],
            "pmt_amp": pmt["amp_cps"], "mot": mot, "u": tw - tc, "frac": frac}


rows = [analyse(p) for n in range(23, 36) for p in sorted((ROOT / "16" / f"{n:04d}").glob("*.h5"))]
xs = np.array([r["x"] for r in rows])
D = np.array([r["D"] for r in rows])
F = np.array([r["pmt"] for r in rows])
L = np.array([r["light"] for r in rows])
M = np.array([r["mot"] for r in rows])

# --- 1. scatter and correlations at the bright positions ---------------------------------------
print("position: corr(D, PMT), corr(D, MOT), corr(D, light); rel sd of D, D/MOT, D/(PMT/light)")
pooled = {"D": [], "D/MOT": [], "D/(PMT/L)": [], "PMT": []}
for x in (-4, -2, 0, 2, 4):
    s = xs == x
    d, f, l, mo = D[s], F[s], L[s], M[s]
    ratio = d / (f / l)
    ratio *= d.mean() / ratio.mean()
    print(f"{x:+3.0f} mm ({s.sum():2d}): {np.corrcoef(d, f)[0,1]:+.2f}, {np.corrcoef(d, mo)[0,1]:+.2f}, "
          f"{np.corrcoef(d, l)[0,1]:+.2f};  {100*d.std()/d.mean():4.0f}%, "
          f"{100*(d/mo).std()/(d/mo).mean():4.0f}%, {100*ratio.std()/ratio.mean():4.0f}%   "
          f"(PMT integral sd {100*f.std()/f.mean():.0f}%)")
    pooled["D"] += list(d / d.mean())
    pooled["D/MOT"] += list((d / mo) / (d / mo).mean())
    pooled["D/(PMT/L)"] += list(ratio / ratio.mean())
    pooled["PMT"] += list(f / f.mean())
print("pooled relative sd: " + ", ".join(f"{k} {100*np.std(v):.0f}%" for k, v in pooled.items()))
print(f"PMT integral vs MOT over all shots: corr {np.corrcoef(F, M)[0,1]:+.2f}; PMT sd {100*F.std()/F.mean():.0f}%, "
      f"MOT sd {100*np.nanstd(M)/np.nanmean(M):.1f}%")

# --- 2. mean dip shape per position -------------------------------------------------------------
grid = np.arange(-0.02, 0.02, 1e-4)
fig, ax = plt.subplots(figsize=(8, 5))
print("\nmean fractional change at the cloud (|u| < 1 ms) and rms away from it (6-20 ms):")
for x in (-10, -8, -6, 0, 6, 8, 10):
    tr = np.array([np.interp(grid, r["u"], r["frac"]) for r in rows if r["x"] == x])
    mean = tr.mean(axis=0)
    centre = np.abs(grid) < 1e-3
    far = (np.abs(grid) > 6e-3)
    print(f"  {x:+3.0f} mm: {100*mean[centre].mean():+.3f}%   (far rms {100*mean[far].std():.3f}%)")
    ax.plot(1e3 * grid, 100 * mean, lw=1.5, label=f"{x:+.0f} mm")
ax.set_xlabel("Time from PMT centre (ms)")
ax.set_ylabel("Light change (%)")
ax.axhline(0, color="grey", lw=1)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG, dpi=110)

# --- 3. line-centre shots from 10/15 --------------------------------------------------------------
print("\n2026-10-15 line-centre shots vs today's 0 mm:")
for n in (7, 8, 9):
    for p in sorted((ROOT / "15" / f"{n:04d}").glob("*.h5")):
        r = analyse(p)
        print(f"  10/15 {n:04d}: mod_offset {1e3*r['off']:.1f} mV, iris {r['x']:+.0f} mm, light {r['light']:.3f} V, "
              f"D {100*r['D']:.2f}%, MOT {1e3*r['mot']:.1f} mV, PMT amp {r['pmt_amp']/1e6:.2f} Mcps, "
              f"PMT integral {1e3*r['pmt']:.1f} kcounts")
s = xs == 0
print(f"  10/16 0 mm (24 shots): mod_offset 0, light {L[s].mean():.3f} V, D {100*D[s].mean():.2f}%, "
      f"MOT {1e3*np.nanmean(M[s]):.1f} mV, PMT amp {np.mean([r['pmt_amp'] for r in rows if r['x']==0])/1e6:.2f} Mcps, "
      f"PMT integral {1e3*F[s].mean():.1f} kcounts")
