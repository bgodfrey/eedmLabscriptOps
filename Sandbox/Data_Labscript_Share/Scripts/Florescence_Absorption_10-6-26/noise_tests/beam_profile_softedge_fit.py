"""Probe light versus iris position, fitted as a Gaussian beam chopped by an opening with soft edges.

Usage: python beam_profile_softedge_fit.py <figure.png>
(The runs are listed in RUNS below: 2026-10-15 0017-0033, 2026-10-16 0023-0035, 2026-10-17 0003-0013.)

Why soft edges: iris_profile_fit.py cuts the beam with perfectly sharp edges. That fitted 10/15
(chi2/dof 50/9) but failed on 10/16 (208/6): there the iris still saw 0.16 V at -10 mm, beyond the
- edge, and a sharp edge allows none. The opening sits upstream of the iris, so its shadow arrives
blurred by diffraction (about sqrt(lambda L) ~ 1 mm over a metre) and beam divergence. Each edge is
therefore an erf step with its own rms blur (s_lo, s_hi): a first version with one shared blur
could not fit both dates, which want a sharp + edge (10/15, bright there) and a soft - edge
(10/16-17, bright there). An edge's blur is only measured on the date the beam is bright at it.

Model at stage offset x (u = horizontal position at the iris, mm):
    I(u) = A exp(-2 (u - x0)^2 / w^2) * T(u),   w = 1/e^2 radius
    T(u) = 1/2 [1 + erf((u - lo) / (sqrt2 s_lo))] * 1/2 [1 + erf((hi - u) / (sqrt2 s_hi))]
    level(x) = mean of I(u) over the iris, u in [x - 1, x + 1] (2 mm top-hat, 81 samples; smooth
               now that the edges are soft)

Each date shows only one flank of the beam (its peak sits at or behind an edge), so a single date
cannot separate the centre from the width. The two dates show opposite flanks, so a joint fit with
one beam width constrains both.

Data per run: per-shot PDA median over the 0.2 s record minus the dark level (52.7 mV, 0010 on
10/15), drift-corrected by the 0 mm visits (linear interpolation in time, as in
iris_profile_fit.py). Weights: each position's standard error combined with a relative floor equal
to the rms spread of that run's 0 mm visits (the drift the interpolation cannot follow) and 5 mV.

Fits (amplitudes always per run, since the probe power differed):
  A. 10/15 alone.
  B. 10/16 + 10/17 (same alignment).
  C. All three: one opening (lo, hi, s_lo, s_hi) and one beam width w, a beam centre per alignment.
     This is the picture "same beam, same opening, only the beam moved". Its chi-square against
     A + B says whether that picture holds.
Several starting centres are tried for each fit (the centre can sit behind an edge) and the best
kept. Errors come from the covariance, scaled by sqrt(chi2/dof) when that exceeds 1.

What to trust (2026-10-19): the opening (lo, hi) and its blurs. The beam centre x0 and width w are
NOT physical: they come out at 15-21 mm radius with the centre 9-21 mm from the opening, i.e. a
30-40 mm beam that moved ~25 mm after an input nudge. The light inside the opening is not a clean
clipped Gaussian (its slow fall-off over 14 mm is the wrong shape), so the Gaussian part of the model
is just a smooth curve that lets the edges be fitted. Do not quote x0, w, fit C's "beam moved" or
the "% inside the opening" lines below.
"""
import os
import sys
import warnings
from pathlib import Path
import h5py
import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.special import erf

warnings.filterwarnings("ignore")
FIG = Path(sys.argv[1])
ROOT = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10")
DARK = 0.0527                                   # V at 70 dB (10/15 0010, probe blocked)
RUNS = [("10/15", ROOT / "15", 17, 33), ("10/16", ROOT / "16", 23, 35), ("10/17", ROOT / "17", 3, 13)]
U = np.linspace(-1.0, 1.0, 81)                  # sample points across the 2 mm iris
BLUE, ORANGE, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"
INF = np.inf


def load_run(day, first, last):
    """Drift-corrected mean light per position, its error, and the run's relative drift floor."""
    x, t, v, f = [], [], [], []
    for n in range(first, last + 1):
        for p in sorted((day / f"{n:04d}").glob("*.h5")):
            with h5py.File(p, "r") as h5:
                x.append(float(h5["globals"].attrs["iris_x_mm"]))
                v.append(np.median(h5["data/traces/TOF_florescence"][()]["values"]) - DARK)
            t.append(os.path.getmtime(p))
            f.append(n)
    x, t, v, f = map(np.array, (x, t, v, f))
    ref_f = sorted(set(f[x == 0]))
    rt = np.array([t[f == k].mean() for k in ref_f])
    rv = np.array([v[f == k].mean() for k in ref_f])
    vc = v / (np.interp(t, rt, rv) / rv.mean())
    floor = rv.std() / rv.mean()
    pos = np.unique(x)
    mean = np.array([vc[x == p].mean() for p in pos])
    sem = np.array([vc[x == p].std(ddof=1) / np.sqrt((x == p).sum()) for p in pos])
    err = np.sqrt(sem ** 2 + (floor * mean) ** 2 + 0.005 ** 2)
    return pos, mean, err, floor


def level(x, a, x0, w, lo, hi, s_lo, s_hi):
    """Soft-edged chopped Gaussian averaged over the iris."""
    u = np.asarray(x, dtype=float)[:, None] + U[None, :]
    beam = np.exp(-2 * (u - x0) ** 2 / w ** 2)
    edges = 0.25 * (1 + erf((u - lo) / (np.sqrt(2) * s_lo))) * (1 + erf((hi - u) / (np.sqrt(2) * s_hi)))
    return a * (beam * edges).mean(axis=1)


def fit(blocks, p0, lo_b, hi_b, x0_index, x0_starts):
    """Least squares over (x, y, err, model) blocks sharing one parameter vector; best of several
    starting values for the beam centre(s) at x0_index."""
    def resid(p):
        return np.concatenate([(y - m(x, p)) / e for x, y, e, m in blocks])
    best = None
    for start in x0_starts:
        q0 = np.array(p0, dtype=float)
        for i, s in zip(x0_index, start):
            q0[i] = s
        r = optimize.least_squares(resid, q0, bounds=(lo_b, hi_b), x_scale="jac", max_nfev=20000)
        if best is None or r.cost < best.cost:
            best = r
    chi2 = np.sum(best.fun ** 2)
    dof = sum(len(b[0]) for b in blocks) - len(p0)
    cov = np.linalg.pinv(best.jac.T @ best.jac)
    err = np.sqrt(np.diag(cov) * max(1.0, chi2 / max(dof, 1)))
    return best.x, err, chi2, dof


def show(label, names, p, e, chi2, dof):
    print(f"\n{label}: chi2/dof {chi2:.1f}/{dof}")
    for n_, v_, e_ in zip(names, p, e):
        print(f"  {n_:18s} {v_:8.2f} +- {e_:.2f}")


data = {name: load_run(day, a, b) for name, day, a, b in RUNS}
for name, (pos, mean, err, floor) in data.items():
    print(f"{name}: {len(pos)} positions {pos.min():+.0f}..{pos.max():+.0f} mm, drift floor {100*floor:.1f}%")
x15, y15, e15, _ = data["10/15"]
x16, y16, e16, _ = data["10/16"]
x17, y17, e17, _ = data["10/17"]
EDGE_LO = [-14, 4, 0.05, 0.05]                  # bounds for lo, hi, s_lo, s_hi
EDGE_HI = [-4, 14, 3, 3]

# --- A. 10/15 alone: p = [A, x0, w, lo, hi, s_lo, s_hi] --------------------------------------------
pA, eA, cA, dA = fit([(x15, y15, e15, lambda x, p: level(x, *p))],
                     [y15.max(), 5, 12, -8, 9, 0.5, 0.5], [0, -40, 3] + EDGE_LO, [INF, 40, 40] + EDGE_HI,
                     [1], [(2,), (5,), (8,), (12,)])
show("A. 10/15", ["A (V)", "x0 (mm)", "w (mm)", "lo (mm)", "hi (mm)", "s_lo (mm)", "s_hi (mm)"],
     pA, eA, cA, dA)

# --- B. 10/16 + 10/17: p = [A16, A17, x0, w, lo, hi, s_lo, s_hi] ------------------------------------
pB, eB, cB, dB = fit([(x16, y16, e16, lambda x, p: level(x, p[0], *p[2:])),
                      (x17, y17, e17, lambda x, p: level(x, p[1], *p[2:]))],
                     [y16.max(), y17.max(), -8, 12, -8, 9, 0.5, 0.5],
                     [0, 0, -40, 3] + EDGE_LO, [INF, INF, 40, 40] + EDGE_HI,
                     [2], [(-5,), (-8,), (-12,), (-20,)])
show("B. 10/16 + 10/17", ["A 10/16 (V)", "A 10/17 (V)", "x0 (mm)", "w (mm)", "lo (mm)", "hi (mm)",
                          "s_lo (mm)", "s_hi (mm)"], pB, eB, cB, dB)

# --- C. one opening, one width, centre per alignment ------------------------------------------------
# p = [A15, A16, A17, x0_15, x0_16, w, lo, hi, s_lo, s_hi]
mC = lambda a_i, x0_i: (lambda x, p: level(x, p[a_i], p[x0_i], *p[5:]))
pC, eC, cC, dC = fit([(x15, y15, e15, mC(0, 3)), (x16, y16, e16, mC(1, 4)), (x17, y17, e17, mC(2, 4))],
                     [pA[0], pB[0], pB[1], 5, -8, 13, -7.7, 8.7, 0.9, 0.2],
                     [0, 0, 0, -40, -40, 3] + EDGE_LO, [INF, INF, INF, 40, 40, 40] + EDGE_HI,
                     [3, 4], [(5, -8), (8, -10), (10, -14), (3, -6)])
show("C. one opening, one beam width", ["A 10/15 (V)", "A 10/16 (V)", "A 10/17 (V)", "x0 10/15 (mm)",
                                        "x0 10/16-17 (mm)", "w (mm)", "lo (mm)", "hi (mm)", "s_lo (mm)",
                                        "s_hi (mm)"], pC, eC, cC, dC)
print(f"\nA + B separately: chi2 {cA + cB:.1f} for {dA + dB} dof;  C: chi2 {cC:.1f} for {dC} dof "
      f"(extra chi2 {cC - cA - cB:+.1f} for {dC - dA - dB} more dof)")
w, lo, hi = pC[5], pC[6], pC[7]
print("\n(x0, w and everything derived from them below are NOT physical; see the docstring. "
      "Trust the opening and its blurs only.)")
print(f"C: opening {hi - lo:.1f} mm wide, centred {0.5 * (lo + hi):+.2f} mm; beam moved "
      f"{pC[4] - pC[3]:+.1f} mm between 10/15 and 10/16")
for label, x0 in (("10/15", pC[3]), ("10/16-17", pC[4])):
    s_ = np.sqrt(2) / w
    frac = 0.5 * (erf(s_ * (hi - x0)) - erf(s_ * (lo - x0)))
    print(f"  {label}: beam centre {x0:+.1f} mm; {100*frac:.0f}% of the horizontal profile inside the opening")
s_ = np.sqrt(2) / w
print(f"  centred on the opening: {100*erf(s_ * (hi - lo) / 2):.0f}% inside")

# --- figure: each date's own fit (A, B), in volts -------------------------------------------------
# Plain volts (drift-corrected to the run's average 0 mm level): the point is where the light sits
# along the stage, not its size. Each run gets its own curve (its own amplitude, the shared shape),
# so in the right panel 10/16 and 10/17 sit at different heights (different probe power).
# The beam centre and width are not constrained (one flank per date), so no "beam alone" curve.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), facecolor=SURFACE)
xx = np.linspace(-13, 13, 600)
shapeA, shapeB = tuple(pA[1:]), tuple(pB[2:])     # (x0, w, lo, hi, s_lo, s_hi)
panels = [("Before realignment (10/15)", shapeA, cA, dA, [("10/15", x15, y15, e15, pA[0], BLUE)]),
          ("After realignment (10/16, 10/17)", shapeB, cB, dB,
           [("10/16", x16, y16, e16, pB[0], BLUE), ("10/17", x17, y17, e17, pB[1], ORANGE)])]
for ax, (title, shape, chi2, dof, sets) in zip(axes, panels):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=MUTED)
    for spine in ax.spines.values():
        spine.set_color(MUTED)
    ax.grid(True, color=MUTED, alpha=0.15)
    _, _, lo_, hi_, slo_, shi_ = shape
    ax.axvspan(lo_, hi_, color=MUTED, alpha=0.07, lw=0)
    # Legend line for the shared fit shape (edges, blurs, chi2); the curves themselves are per run
    ax.plot([], [], color=INK, lw=2,
            label=f"Fit: edges {lo_:+.1f} (blur {slo_:.1f}) / {hi_:+.1f} (blur {shi_:.1f}) mm, "
                  f"χ²/dof {chi2:.0f}/{dof}")
    for name, x, y, e, a, c in sets:            # each run with its own fitted amplitude
        ax.plot(xx, level(xx, a, *shape), color=c, lw=2)
        ax.errorbar(x, y, yerr=e, fmt="o", ms=7, color=c, mec=SURFACE, mew=1.5, capsize=3, lw=1.5,
                    label=name)
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xlim(-13, 13)
    ax.set_ylim(bottom=-0.05 * ax.get_ylim()[1])
    ax.set_xlabel("Stage offset (mm; + = left, facing N)", color=INK)
    ax.set_ylabel("Probe light (V)", color=INK)
    ax.set_title(title, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper left", labelcolor=INK)
fig.tight_layout()
fig.savefig(FIG, dpi=130, facecolor=SURFACE)
plt.show()
print(f"\nfigure: {FIG}")
