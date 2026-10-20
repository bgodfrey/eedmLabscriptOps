"""Upward-pass dip versus probe offset, every shot, with the PMT for comparison (2026-10-16 0036, 0037).

Usage: python offset_scan_plot.py <helper dir> <day dir> <dark mV> <figure.png> <folder> [<folder> ...]

Per shot: the dip D from the fixed scan shape (quadratic baseline x (1 - D Gaussian), centre
-0.13 ms from the PMT, sigma 1.48 ms; as in iris_atom_scan.py) and the PMT integral F (amplitude x
transit width, not corrected for dead time). Iris at 0 mm, modulation off.

Left: D against mod_offset, one dot per shot, folders in different markers, the median per offset
and folder as a bar. Right: D against F. If the low shots ("dropouts") are atoms that absorb and
scatter less (fewer atoms, or the probe off resonance), they fall on the same curve as the normal
shots; the curve bends because the PMT loses counts at high rates.
Shots below 60% of their offset's median at >= 6 mV are drawn hollow as dropouts.
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
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DAY, DARK, FIG = Path(sys.argv[2]), float(sys.argv[3]) * 1e-3, Path(sys.argv[4])
FOLDERS = sys.argv[5:]
OFFSET, WIDTH = -0.13e-3, 1.48e-3
BLUE, ORANGE, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"
COLORS = {FOLDERS[0]: ORANGE, FOLDERS[-1]: BLUE} if len(FOLDERS) > 1 else {FOLDERS[0]: BLUE}

rows = []
for folder in FOLDERS:
    for p in sorted((DAY / folder).glob("*.h5")):
        s = load_shot(p)
        with h5py.File(p, "r") as h5:
            off = 1e3 * float(h5["globals"].attrs["mod_offset"])
        pmt = cse.pmt_cloud_fit(s, 6125.0)
        if s["pda_v"].max() >= 9.95 or not pmt["ok"]:
            continue
        tc = pmt["center_s"]
        k = int(round(s["fs_pda"] / 10e3))
        m = len(s["pda_v"]) // k * k
        tb = s["t_pda"][:m].reshape(-1, k).mean(axis=1)
        vb = s["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
        win = (tb > s["enable_start_s"] + 0.002) & (tb < s["enable_stop_s"] - 0.002)

        def model(t, c0, c1, c2, d):
            u = t - tc
            return (c0 + c1 * u + c2 * u ** 2) * (1 - d * np.exp(-0.5 * ((t - tc - OFFSET) / WIDTH) ** 2))

        pf, _ = optimize.curve_fit(model, tb[win], vb[win], p0=[np.median(vb[win]), 0, 0, 0.03])
        rows.append((folder, off, 100 * pf[3], pmt["amp_cps"] * pmt["sigma_s"] / 1e3))

rows_a = np.array([(r[1], r[2], r[3]) for r in rows])
folder_of = np.array([r[0] for r in rows])
# Dropouts: below 60% of the median at that offset (line-centre region, all folders pooled)
med = {o: np.median(rows_a[rows_a[:, 0] == o, 1]) for o in np.unique(rows_a[:, 0])}
drop = np.array([o >= 6 and d < 0.6 * med[o] for o, d, _ in rows_a])
for f in FOLDERS:
    sel = folder_of == f
    for o in np.unique(rows_a[sel, 0]):
        s2 = sel & (rows_a[:, 0] == o)
        print(f"{f} {o:5.1f} mV: {s2.sum():2d} shots, median {np.median(rows_a[s2, 1]):.2f}%, "
              f"dropouts {np.sum(drop & s2)}")

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5), width_ratios=(1.3, 1), facecolor=SURFACE)
for a_ in (ax, ax2):
    a_.set_facecolor(SURFACE)
    a_.tick_params(colors=MUTED)
    for spine in a_.spines.values():
        spine.set_color(MUTED)
    a_.grid(True, color=MUTED, alpha=0.15)
rng = np.random.default_rng(1)
for f in FOLDERS:
    sel = folder_of == f
    c = COLORS.get(f, BLUE)
    shift = -0.12 if (len(FOLDERS) > 1 and f == FOLDERS[0]) else 0.12
    jit = shift + rng.uniform(-0.08, 0.08, sel.sum())
    x, d, fl = rows_a[sel, 0] + jit, rows_a[sel, 1], rows_a[sel, 2]
    dr = drop[sel]
    ax.plot(x[~dr], d[~dr], "o", ms=6, color=c, alpha=0.8, mec="none", label=f"{f}")
    ax.plot(x[dr], d[dr], "o", ms=6, mfc=SURFACE, mec=c, mew=1.5,
            label=f"{f}, dropout" if dr.any() else None)          # no legend entry when none
    for o in np.unique(rows_a[sel, 0]):
        mm = np.median(rows_a[sel & (rows_a[:, 0] == o), 1])
        ax.plot([o + shift - 0.25, o + shift + 0.25], [mm, mm], color=c, lw=2.5)
    ax2.plot(fl[~dr], d[~dr], "o", ms=6, color=c, alpha=0.8, mec="none")
    ax2.plot(fl[dr], d[dr], "o", ms=6, mfc=SURFACE, mec=c, mew=1.5)
ax.set_xlabel("mod_offset (mV)", color=INK)
ax.set_ylabel("Peak dip at 0 mm (%)", color=INK)
ax.set_title("Dip vs probe offset (bars = median)", color=INK)
ax.set_ylim(0, None)
ax.legend(frameon=False, fontsize=9, loc="lower right", labelcolor=INK)
ax2.set_xlabel("PMT integral (k counts, no dead-time correction)", color=INK)
ax2.set_ylabel("Peak dip at 0 mm (%)", color=INK)
ax2.set_title("Dip vs PMT, every shot", color=INK)
ax2.set_xlim(0, None)
ax2.set_ylim(0, None)
fig.tight_layout()
fig.savefig(FIG, dpi=130, facecolor=SURFACE)
print(f"figure: {FIG}")
