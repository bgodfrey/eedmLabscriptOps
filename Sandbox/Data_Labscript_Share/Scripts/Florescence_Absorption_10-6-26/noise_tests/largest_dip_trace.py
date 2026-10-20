"""Find the single shot with the deepest absorption dip (%) in the 10/16 and 10/17 iris scans, and plot it.

Usage: python largest_dip_trace.py <helper dir> <figure.png> [--show]
(Runs: 2026-10-16 0023-0035 and 2026-10-17 0003-0013, as in beam_profile_softedge_fit.py.)

Dip per shot, model-free (no PMT normalisation: this is the dip the trace itself shows):
- PDA trace minus the 52.7 mV dark level, averaged to 10 kHz, inside the jump-enabled window
  (2 ms guards), as in dip_vs_offset_hv.py.
- Baseline: quadratic in time, fitted only more than BASE_FAR_S from the PMT cloud centre, so the
  atoms (including the leading tail) do not pull it down much.
- Fractional dip = 1 - light / baseline (positive = absorption), averaged in 1 mm bins
  (0.43 ms at 2.30 m/s, the sheet height). The shot's peak dip is the largest bin within
  +-PEAK_SEARCH_S of the cloud centre.
The fixed-shape D of dip_vs_offset_hv.py (quadratic x (1 - D Gaussian), 10/17 shape) is printed
alongside as a cross-check of the ranking. Shots with a saturated PDA (>= 9.95 V) or a failed PMT
fit are skipped, and so are shots whose baseline light at the cloud is below LIGHT_MIN: outside the
lit range (e.g. 10/16 +10 mm, 10-27 mV) a few mV of noise or dark-level error over almost no light
gives fake "dips" of 6-29%.

Figure: top = the deepest shot's light (V) with its baseline; bottom = its fractional dip (%),
10 kHz (faint) and 1 mm bins (dots), with the peak bin marked.
"""
import sys
import warnings
from pathlib import Path
import h5py
import numpy as np
import matplotlib
if "--show" not in sys.argv:
    matplotlib.use("Agg")                      # file only, unless asked to open a window
import matplotlib.pyplot as plt
from scipy import optimize

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
FIG = Path(sys.argv[2])
ROOT = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10")
RUNS = [("10/16", ROOT / "16", 23, 35), ("10/17", ROOT / "17", 3, 13)]
DARK = 0.0527                                  # V at 70 dB (10/15 0010, probe blocked)
V_PROBE = 2.30                                 # m/s through the probe sheet
BIN_S = 1e-3 / V_PROBE                         # 1 mm of cloud height per bin
BASE_FAR_S = 0.010                             # baseline fitted only at |t - t_c| > this
PEAK_SEARCH_S = 0.005                          # look for the peak bin within this of t_c
LIGHT_MIN = 0.10                               # V; need ~20x the 5 mV dark uncertainty for a % dip
OFFSET, WIDTH = 0.05e-3, 1.51e-3               # fixed dip shape (10/17, as in dip_vs_offset_hv.py)
BLUE, ORANGE, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"

# --- every shot: model-free peak dip and the fixed-shape D --------------------------------------
rows = []
for name, day, first, last in RUNS:
    for n in range(first, last + 1):
        for p in sorted((day / f"{n:04d}").glob("*.h5")):
            shot = load_shot(p)
            with h5py.File(p, "r") as h5:
                x = float(h5["globals"].attrs["iris_x_mm"])
            pmt = cse.pmt_cloud_fit(shot, shot["mod_freq"] or 6125.0)
            if shot["pda_v"].max() >= 9.95 or not pmt["ok"]:
                continue
            # 10 kHz averages inside the enabled window, dark subtracted
            k = int(round(shot["fs_pda"] / 10e3))
            m = len(shot["pda_v"]) // k * k
            tb = shot["t_pda"][:m].reshape(-1, k).mean(axis=1)
            vb = shot["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
            win = (tb > shot["enable_start_s"] + 0.002) & (tb < shot["enable_stop_s"] - 0.002)
            tw, vw, tc = tb[win], vb[win], pmt["center_s"]
            u = tw - tc                                        # time from the cloud centre
            # Baseline from far outside the cloud, then the fractional dip
            far = np.abs(u) > BASE_FAR_S
            bc = np.polyfit(u[far], vw[far], 2)
            base = np.polyval(bc, u)
            frac = 1 - vw / base
            # 1 mm bins, centred on the cloud centre so one bin sits on the peak
            edges = np.arange(-PEAK_SEARCH_S - BIN_S / 2, PEAK_SEARCH_S + BIN_S, BIN_S)
            idx = np.digitize(u, edges) - 1
            nb = len(edges) - 1
            bins = np.array([frac[idx == i].mean() if np.any(idx == i) else np.nan for i in range(nb)])
            ipk = np.nanargmax(bins)

            # Fixed-shape D for comparison (quadratic baseline x (1 - D Gaussian))
            def fixed(t, c0, c1, c2, d):
                uu = t - tc
                return (c0 + c1 * uu + c2 * uu ** 2) * (1 - d * np.exp(-0.5 * ((t - tc - OFFSET) / WIDTH) ** 2))
            pf, _ = optimize.curve_fit(fixed, tw, vw, p0=[np.median(vw), 0, 0, 0.01])
            rows.append({"run": name, "file": p, "folder": n, "x": x, "peak": bins[ipk],
                         "t_peak": 0.5 * (edges[ipk] + edges[ipk + 1]), "D": pf[3],
                         "F": pmt["amp_cps"] * pmt["sigma_s"], "light": np.polyval(bc, 0),
                         "tw": tw, "vw": vw, "base": base, "tc": tc, "frac": frac, "u": u,
                         "bins": bins, "edges": edges})

# --- ranking ------------------------------------------------------------------------------------
dim = [r for r in rows if r["light"] < LIGHT_MIN]
rows = [r for r in rows if r["light"] >= LIGHT_MIN]
print(f"{len(dim)} shots below {1e3*LIGHT_MIN:.0f} mV of light dropped (iris "
      f"{', '.join(f'{x:+.0f}' for x in sorted({r['x'] for r in dim}))} mm)")
rows.sort(key=lambda r: r["peak"], reverse=True)
print(f"{len(rows)} shots. Top 10 by model-free peak dip (1 mm bin):")
print("  run    folder  file                                           iris mm  peak %   D %   "
      "light V  PMT F (rel)")
f_med = np.median([r["F"] for r in rows])
for r in rows[:10]:
    print(f"  {r['run']}  {r['folder']:04d}   {r['file'].name:46s} {r['x']:+5.1f}  {100*r['peak']:5.2f}  "
          f"{100*r['D']:5.2f}   {r['light']:6.3f}   {r['F']/f_med:5.2f}")
for name, *_ in RUNS:                          # deepest shot of each run, and the run's typical dip
    rr = [r for r in rows if r["run"] == name]
    r0 = [r for r in rr if r["x"] == 0]
    print(f"{name}: deepest {rr[0]['file'].name} at {rr[0]['x']:+.0f} mm, peak {100*rr[0]['peak']:.2f}%; "
          f"0 mm shots: median peak {100*np.median([r['peak'] for r in r0]):.2f}%, "
          f"median PMT F {np.median([r['F'] for r in r0])/f_med:.2f}")
by_d = max(rows, key=lambda r: r["D"])
print(f"Largest fixed-shape D: {by_d['run']} {by_d['file'].name} at {by_d['x']:+.1f} mm, "
      f"D {100*by_d['D']:.2f}% (peak bin {100*by_d['peak']:.2f}%)")

# --- figure: the deepest shot -------------------------------------------------------------------
r = rows[0]
t_ms = 1e3 * r["tw"]
fig, (top, bot) = plt.subplots(2, 1, figsize=(10, 6.4), sharex=True, height_ratios=(1, 1),
                               facecolor=SURFACE)
for a_ in (top, bot):
    a_.set_facecolor(SURFACE)
    a_.tick_params(colors=MUTED)
    for spine in a_.spines.values():
        spine.set_color(MUTED)
    a_.grid(True, color=MUTED, alpha=0.15)
    a_.axvline(1e3 * r["tc"], color=MUTED, lw=1, ls=":")       # PMT cloud centre
    # Shade the region the baseline was fitted on (|t - t_c| > BASE_FAR_S)
    a_.axvspan(t_ms.min(), 1e3 * (r["tc"] - BASE_FAR_S), color=MUTED, alpha=0.07, lw=0)
    a_.axvspan(1e3 * (r["tc"] + BASE_FAR_S), t_ms.max(), color=MUTED, alpha=0.07, lw=0)
top.plot(t_ms, r["vw"], color=BLUE, lw=1.2, label="Probe light (10 kHz)")
top.plot(t_ms, r["base"], color=ORANGE, lw=1.8, ls="--", label="Baseline (fitted in shaded region)")
top.set_ylabel("Probe light (V)", color=INK)
top.set_title(f"Deepest dip: {r['run']} {r['file'].name}, iris {r['x']:+.0f} mm", color=INK)
top.legend(frameon=False, fontsize=9, loc="lower left", labelcolor=INK)
bot.axhline(0, color=MUTED, lw=1)
bot.plot(1e3 * (r["u"] + r["tc"]), 100 * r["frac"], color=BLUE, lw=0.8, alpha=0.35, label="10 kHz")
bc_ms = 1e3 * (0.5 * (r["edges"][1:] + r["edges"][:-1]) + r["tc"])
bot.plot(bc_ms, 100 * r["bins"], "o", ms=6, color=BLUE, mec=SURFACE, mew=1.2, label="1 mm bins")
bot.plot(1e3 * (r["t_peak"] + r["tc"]), 100 * r["peak"], "o", ms=11, mfc="none", mec=ORANGE, mew=2,
         label=f"Peak {100*r['peak']:.2f}%")
bot.set_ylabel("Dip (%)", color=INK)
bot.set_xlabel("Time in record (ms)", color=INK)
bot.legend(frameon=False, fontsize=9, loc="upper right", labelcolor=INK)
fig.tight_layout()
fig.savefig(FIG, dpi=130, facecolor=SURFACE)
print(f"figure: {FIG}")
if "--show" in sys.argv:
    plt.show()
