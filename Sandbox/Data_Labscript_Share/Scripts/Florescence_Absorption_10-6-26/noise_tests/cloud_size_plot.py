"""Size of the upward cloud at the probe: horizontal and vertical absorption profiles on one axis.

Usage: python cloud_size_plot.py <helper dir> <day dir> <dark mV> <scan first> <scan last>
                                 <line-centre folder> <figure.png>
Example (2026-10-16): ... 2026/10/16 52.7 23 35 36 cloud_size_0023-0036.png

Left panel, both profiles normalised to their own peak and centred on the cloud:
- Horizontal (x): the iris scan (0023-0035). Per shot the fixed-shape dip D (as in
  iris_atom_scan.py: quadratic baseline x (1 - D Gaussian), shape -0.13 ms / 1.48 ms from the
  scan), normalised by the shot's PMT integral over its relative light level, then averaged per
  position and fitted with a Gaussian averaged over the 2 mm iris (positions |x| <= 8 mm only).
  The scan ran at mod_offset 0; detuning scales every position alike, so the shape is unaffected.
- Vertical (z): the line-centre shots of the interleaved check (mod_offset 0.0095, iris at 0 mm).
  Each trace is divided by a quadratic baseline fitted more than 6 ms from the PMT cloud time,
  giving the fractional dip versus time; the 12 traces are averaged on a common grid aligned to
  the PMT centre. Time becomes height through the speed at the probe: z = v (t_c - t), v = 2.30 m/s
  (the top of the cloud arrives first, so earlier time is higher). The 1 mm sheet height adds
  about 0.08 mm^2 to sigma_z^2, which is negligible.
  Its error comes from the spread of single-shot Gaussian fits (the averaged trace has correlated
  noise, so its own fit error is too small).

Right panel, to scale: the cloud's 1 sigma and 2 sigma ellipses in the x-z plane when it is centred
in the probe sheet, over the sheet drawn 1 mm tall with its brightness set by the light the iris
saw on 2026-10-16 (that is the light reaching the PDA; whether an aperture clips it before or after
the atoms is still open).
"""
import sys
import warnings
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.colors import LinearSegmentedColormap
from scipy import optimize
from scipy.special import erf

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse

warnings.filterwarnings("ignore")
DAY, DARK = Path(sys.argv[2]), float(sys.argv[3]) * 1e-3
FIRST, LAST, LC_FOLDER, FIG = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), Path(sys.argv[7])
D_IRIS_MM, FIT_MAX_MM = 2.0, 8.0
OFFSET, WIDTH = -0.13e-3, 1.48e-3          # common dip shape from iris_atom_scan.py
V_PROBE = 2.30                             # m/s through the probe (0130 crossing times)
T_RELEASE_MS = 87                          # upward crossing, ms after release (0130)
LINE_CENTRE = 0.0095                       # mod_offset at line centre (V)
# Colors as in the other noise_tests figures: horizontal blue, vertical orange, ink, muted, surface
BLUE, ORANGE, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"


def pda_window(shot):
    """PDA at 10 kHz inside the jump-enabled window (2 ms guards), dark subtracted."""
    k = int(round(shot["fs_pda"] / 10e3))
    m = len(shot["pda_v"]) // k * k
    tb = shot["t_pda"][:m].reshape(-1, k).mean(axis=1)
    vb = shot["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
    win = (tb > shot["enable_start_s"] + 0.002) & (tb < shot["enable_stop_s"] - 0.002)
    return tb[win], vb[win]


# --- horizontal: iris scan ---------------------------------------------------------------------------
rows = []
for n in range(FIRST, LAST + 1):
    for p in sorted((DAY / f"{n:04d}").glob("*.h5")):
        shot = load_shot(p)
        with h5py.File(p, "r") as h5:
            x = float(h5["globals"].attrs["iris_x_mm"])
        pmt = cse.pmt_cloud_fit(shot, shot["mod_freq"] or 6125.0)
        if shot["pda_v"].max() >= 9.95 or not pmt["ok"]:
            continue
        tw, vw = pda_window(shot)
        tc = pmt["center_s"]

        def fixed(t, c0, c1, c2, d):
            u = t - tc
            return (c0 + c1 * u + c2 * u ** 2) * (1 - d * np.exp(-0.5 * ((t - tc - OFFSET) / WIDTH) ** 2))

        pf, _ = optimize.curve_fit(fixed, tw, vw, p0=[np.median(vw), 0, 0, 0.01])
        rows.append({"x": x, "D": pf[3], "light": pf[0], "F": pmt["amp_cps"] * pmt["sigma_s"]})
xs = np.array([r["x"] for r in rows])
positions = np.unique(xs)
# PMT normalisation, keeping each position's mean (see iris_atom_scan.py)
light_mean = {x: np.mean([r["light"] for r in rows if r["x"] == x]) for x in positions}
q = np.array([r["F"] / (r["light"] / light_mean[r["x"]]) for r in rows])
dn = np.array([r["D"] for r in rows]) * q.mean() / q
mean = np.array([dn[xs == x].mean() for x in positions])
sem = np.array([dn[xs == x].std(ddof=1) / np.sqrt((xs == x).sum()) for x in positions])
light = np.array([light_mean[x] for x in positions])


def profile(x, a, xc, sx):
    """Gaussian averaged over the iris width (analytic)."""
    h = D_IRIS_MM / 2
    return a * sx * np.sqrt(np.pi / 2) / D_IRIS_MM * (erf((x + h - xc) / (np.sqrt(2) * sx))
                                                     - erf((x - h - xc) / (np.sqrt(2) * sx)))


inside = np.abs(positions) <= FIT_MAX_MM
ph, ch = optimize.curve_fit(profile, positions[inside], mean[inside], p0=(mean.max(), 0, 3),
                            sigma=sem[inside], absolute_sigma=True)
chi2 = np.sum(((mean[inside] - profile(positions[inside], *ph)) / sem[inside]) ** 2)
eh = np.sqrt(np.diag(ch) * max(1, chi2 / (inside.sum() - 3)))
sx, sx_err, xc = abs(ph[2]), eh[2], ph[1]
print(f"horizontal: centre {xc:+.2f} mm, sigma_x {sx:.2f} +- {sx_err:.2f} mm")

# --- vertical: line-centre traces --------------------------------------------------------------------
grid = np.arange(-0.008, 0.008 + 1e-9, 1e-4)        # s from the PMT centre
traces, widths = [], []
for p in sorted((DAY / f"{LC_FOLDER:04d}").glob("*.h5")):
    with h5py.File(p, "r") as h5:
        if not np.isclose(float(h5["globals"].attrs["mod_offset"]), LINE_CENTRE):
            continue
    shot = load_shot(p)
    pmt = cse.pmt_cloud_fit(shot, shot["mod_freq"] or 6125.0)
    if shot["pda_v"].max() >= 9.95 or not pmt["ok"]:
        continue
    tw, vw = pda_window(shot)
    u = tw - pmt["center_s"]
    far = np.abs(u) > 0.006                          # baseline from outside the cloud
    frac = 1 - vw / np.polyval(np.polyfit(u[far], vw[far], 2), u)   # positive = absorption
    tr = np.interp(grid, u, frac)
    traces.append(tr)
    try:                                             # single-shot width, for the error bar
        g1, _ = optimize.curve_fit(lambda t, a, t0, s: a * np.exp(-0.5 * ((t - t0) / s) ** 2),
                                   grid, tr, p0=(tr.max(), 0, 1.5e-3))
        widths.append(abs(g1[2]))
    except RuntimeError:
        pass
traces = np.array(traces)
v_mean = traces.mean(axis=0)
v_sem = traces.std(axis=0, ddof=1) / np.sqrt(len(traces))
pv, _ = optimize.curve_fit(lambda t, a, t0, s, c: c + a * np.exp(-0.5 * ((t - t0) / s) ** 2),
                           grid, v_mean, p0=(v_mean.max(), 0, 1.5e-3, 0))
st = abs(pv[2])
st_err = np.std(widths, ddof=1) / np.sqrt(len(widths))
sz, sz_err = 1e3 * V_PROBE * st, 1e3 * V_PROBE * st_err
print(f"vertical: {len(traces)} line-centre shots, peak {100*pv[0]:.2f}%, sigma_t {1e3*st:.2f} +- "
      f"{1e3*st_err:.2f} ms -> sigma_z {sz:.2f} +- {sz_err:.2f} mm (centre {1e3*pv[1]:+.2f} ms from PMT)")
print(f"aspect sigma_z / sigma_x = {sz/sx:.2f} +- {np.hypot(sz_err/sx, sz*sx_err/sx**2):.2f}")

# --- figure ------------------------------------------------------------------------------------------
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5), width_ratios=(1.15, 1), facecolor=SURFACE)
for a_ in (ax, ax2):
    a_.set_facecolor(SURFACE)
    a_.tick_params(colors=MUTED)
    for spine in a_.spines.values():
        spine.set_color(MUTED)

# Left: normalised profiles against distance from the cloud centre
ax.grid(True, color=MUTED, alpha=0.15)
z = -1e3 * V_PROBE * (grid - pv[1])                 # mm, positive = up
vn = (v_mean - pv[3]) / pv[0]
ax.fill_between(z, vn - v_sem / pv[0], vn + v_sem / pv[0], color=ORANGE, alpha=0.2, lw=0)
ax.plot(z, vn, color=ORANGE, lw=1.6, label=f"Vertical (z): σ {sz:.1f} ± {sz_err:.1f} mm")
zz = np.linspace(-12, 12, 400)
ax.plot(zz, np.exp(-0.5 * (zz / sz) ** 2), color=ORANGE, lw=1.2, ls="--")
ax.errorbar(positions[inside] - xc, mean[inside] / ph[0], yerr=sem[inside] / ph[0], fmt="o", ms=7,
            color=BLUE, mec=SURFACE, mew=1.5, capsize=3, lw=1.5,
            label=f"Horizontal (x): σ {sx:.1f} ± {sx_err:.1f} mm")
ax.plot(zz, profile(zz + xc, *ph) / ph[0], color=BLUE, lw=1.2, ls="--")
ax.axhline(0, color=MUTED, lw=1)
ax.set_xlim(-11, 11)
ax.set_xlabel("Distance from cloud centre (mm)", color=INK)
ax.set_ylabel("Absorption (normalised)", color=INK)
ax.set_title(f"Cloud profiles, {T_RELEASE_MS} ms after release", color=INK)
ax.legend(frameon=False, fontsize=9, loc="upper left", labelcolor=INK)

# Right: to-scale cross-section, cloud ellipses over the probe sheet
xf = np.linspace(-10, 10, 400)
sheet = np.interp(xf, positions, light / light.max())[None, :]
cmap = LinearSegmentedColormap.from_list("sheet", [SURFACE, MUTED])
ax2.imshow(sheet, extent=(-10, 10, -0.5, 0.5), cmap=cmap, vmin=0, vmax=1, aspect="auto", zorder=1)
for k, alpha in ((2, 0.18), (1, 0.38)):
    ax2.add_patch(Ellipse((xc, 0), 2 * k * sx, 2 * k * sz, color=BLUE, alpha=alpha, lw=0, zorder=2))
ax2.plot(xc, 0, "+", color=INK, ms=10, zorder=3)
ax2.annotate("1σ", (xc + sx * 0.72, sz * 0.72), color=INK, fontsize=9, zorder=3)
ax2.annotate("2σ", (xc + 2 * sx * 0.72, 2 * sz * 0.72), color=INK, fontsize=9, zorder=3)
ax2.annotate("Probe sheet, 1 mm tall\n(shade = light at the iris)", xy=(-6, 0.5), xytext=(-10.5, 6.8),
             color=INK, fontsize=9, arrowprops=dict(arrowstyle="->", color=MUTED))
ax2.set_xlim(-11, 11)
ax2.set_ylim(-9, 9)
ax2.set_aspect("equal")
ax2.grid(True, color=MUTED, alpha=0.15)
ax2.set_xlabel("Stage offset x (mm; + = left, facing N)", color=INK)
ax2.set_ylabel("Height z (mm)", color=INK)
ax2.set_title("To scale, cloud centred in the sheet", color=INK)
fig.tight_layout()
fig.savefig(FIG, dpi=130, facecolor=SURFACE)
print(f"figure: {FIG}")
