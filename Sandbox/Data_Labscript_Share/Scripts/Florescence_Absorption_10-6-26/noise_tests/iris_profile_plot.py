"""Coarse probe-beam profile through a ~2 mm iris (no launch), with a Gaussian through the middle.

Usage: python iris_profile_plot.py <day dir> <dark mV> <figure.png> <folder>:<offset mm> [...]

Example (2026-10-15):
    python iris_profile_plot.py C:/Experiments/.../2026/10/15 52.7 iris_profile_0012-0016.png
        0012:-12 0013:-6 0014:0 0015:6 0016:12

Each shot's level is the median PDA voltage over its 0.2 s record minus the dark level (blocked-probe
run 0010). Per-shot levels are drawn as small dots, position means with their shot-to-shot spread
as larger markers. The dashed curve is the Gaussian that passes exactly through the three middle
means (a quadratic in ln I), extended to the outer positions: where the measured light falls far
below it, something clips the beam.
"""
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DAY, DARK, FIG = Path(sys.argv[1]), float(sys.argv[2]) * 1e-3, Path(sys.argv[3])
# Colors: data blue, model orange, ink and muted grey, warm-white surface (same as dark_vs_lit.py)
BLUE, ORANGE, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"

xs, shots = [], []
for arg in sys.argv[4:]:
    folder, offset = arg.split(":")
    levels = []
    for p in sorted((DAY / folder).glob("*.h5")):
        with h5py.File(p, "r") as h5:
            levels.append(np.median(h5["data/traces/TOF_florescence"][()]["values"]) - DARK)
    xs.append(float(offset))
    shots.append(np.array(levels))
xs = np.array(xs)
means = np.array([s.mean() for s in shots])
spreads = np.array([s.std() for s in shots])

# Gaussian through the three middle points: ln I = c0 + c1 x + c2 x^2, exactly determined
mid = np.argsort(np.abs(xs))[:3]
c2, c1, c0 = np.polyfit(xs[mid], np.log(means[mid]), 2)
x0, w = -c1 / (2 * c2), np.sqrt(-2 / c2)            # centre and 1/e^2 radius
xx = np.linspace(xs.min() - 2, xs.max() + 2, 400)
model = np.exp(c0 + c1 * xx + c2 * xx ** 2)

fig, ax = plt.subplots(figsize=(8, 5), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
ax.plot(xx, model, ls="--", color=ORANGE, lw=2,
        label=f"Gaussian through middle 3 (centre {x0:+.1f} mm, w {w:.0f} mm)")
for x, s in zip(xs, shots):
    ax.plot(np.full(len(s), x), s, "o", ms=4, color=BLUE, alpha=0.45, mec="none")
ax.errorbar(xs, means, yerr=spreads, fmt="o", ms=9, color=BLUE, mec=SURFACE, mew=2,
            capsize=4, lw=2, label="Mean of 4 shots (± spread)")
ax.axhline(0, color=MUTED, lw=1)
# Point out the clipped edge: what the Gaussian predicts there against what was measured
edge = np.argmax(xs)
pred = np.exp(c0 + c1 * xs[edge] + c2 * xs[edge] ** 2)
ax.annotate(f"predicted {pred:.1f} V,\nmeasured {1e3*means[edge]:.0f} mV",
            xy=(xs[edge], means[edge]), xytext=(xs[edge] - 7.5, 0.35 * pred),
            color=INK, fontsize=10, arrowprops=dict(arrowstyle="->", color=MUTED))
ax.set_xlabel("Stage offset (mm; + = left, facing N)", color=INK)
ax.set_ylabel("PDA above dark (V, 70 dB)", color=INK)
ax.set_title("Iris profile, no launch (0012–0016)", color=INK)
ax.tick_params(colors=MUTED)
for spine in ax.spines.values():
    spine.set_color(MUTED)
ax.grid(True, color=MUTED, alpha=0.15)
ax.legend(frameon=False, fontsize=9, loc="upper left", labelcolor=INK)
fig.tight_layout()
fig.savefig(FIG, dpi=130, facecolor=SURFACE)
print(f"centre {x0:+.2f} mm, w {w:.1f} mm; figure: {FIG}")
