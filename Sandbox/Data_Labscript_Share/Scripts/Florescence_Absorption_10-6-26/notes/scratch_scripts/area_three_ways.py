"""10/17 0003-0013: the per-position absorption area three ways, PMT-normalised identically.
  old   = D_fixed x sqrt(2 pi) x 1.51 ms   (template: centre and width fixed)
  free  = D_free  x sqrt(2 pi) x w_free    (Gaussian with its own centre and width per shot, no cut)
  integ = sum of (1 - V/baseline) dt over |u| <= 10 ms, baseline quadratic at |u| > 10 ms
The free Gaussian adapts to the dip's width but not to its tail, and averages noise better than the
raw sum. If it agrees with the integral in the wings, the template is the biased one.
Usage: python area_three_ways.py <helper dir>
"""
import sys
from pathlib import Path
import h5py
import numpy as np
from scipy import optimize
from scipy.special import erf

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse

DAY = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\17")
DARK, OFFSET, WIDTH, HALF = 0.0527, 0.05e-3, 1.51e-3, 0.010
V, SIG_ISO = 2.30, 1.405e-9 * 1e-4 / 1.2
R2P = np.sqrt(2 * np.pi)


def model(t, tc, c0, c1, c2, d, centre, width):
    u = t - tc
    return (c0 + c1 * u + c2 * u ** 2) * (1 - d * np.exp(-0.5 * ((t - centre) / width) ** 2))


rows = []
for n in range(3, 14):
    for p in sorted((DAY / f"{n:04d}").glob("*.h5")):
        shot = load_shot(p)
        with h5py.File(p, "r") as h5:
            x = float(h5["globals"].attrs["iris_x_mm"])
        pmt = cse.pmt_cloud_fit(shot, shot["mod_freq"] or 6125.0)
        if shot["pda_v"].max() >= 9.95 or not pmt["ok"]:
            continue
        k = int(round(shot["fs_pda"] / 10e3)); m = len(shot["pda_v"]) // k * k
        tb = shot["t_pda"][:m].reshape(-1, k).mean(axis=1)
        vb = shot["pda_v"][:m].reshape(-1, k).mean(axis=1) - DARK
        win = (tb > shot["enable_start_s"] + 0.002) & (tb < shot["enable_stop_s"] - 0.002)
        tw, vw, tc = tb[win], vb[win], pmt["center_s"]
        u = tw - tc
        # template (old)
        pf, _ = optimize.curve_fit(lambda t, c0, c1, c2, d: model(t, tc, c0, c1, c2, d, tc + OFFSET, WIDTH),
                                   tw, vw, p0=[np.median(vw), 0, 0, 0.01])
        # free Gaussian, same bounds as iris_atom_scan.py pass 1, no quality cut
        try:
            pg, _ = optimize.curve_fit(lambda t, *a: model(t, tc, *a), tw, vw,
                                       p0=[np.median(vw), 0, 0, 0.02, tc, 1.3e-3],
                                       bounds=([-np.inf] * 3 + [-0.5, tc - 0.005, 0.3e-3],
                                               [np.inf] * 3 + [0.5, tc + 0.005, 5e-3]))
            free = pg[3] * R2P * pg[5]
        except (RuntimeError, ValueError):
            free = np.nan
        # model-free integral
        far = np.abs(u) > HALF
        base = np.polyval(np.polyfit(u[far], vw[far], 2), u)
        integ = np.sum((1 - vw / base)[~far]) * np.median(np.diff(u))
        rows.append((x, pf[3] * R2P * WIDTH, free, integ, pf[0], pmt["amp_cps"] * pmt["sigma_s"]))

a = np.array(rows)
x, old, free, integ, light, F = a.T
pos = np.unique(x)
lm = {p: light[x == p].mean() for p in pos}
q = F / (light / np.array([lm[v] for v in x]))
sc = q.mean() / q
old, free, integ = old * sc, free * sc, integ * sc
print("  x mm   old (% ms)        free Gaussian      integrated        free/old  integ/old")
sums = {"old": 0, "free": 0, "integ": 0}
for p in pos:
    s = x == p
    out = []
    for name, v in (("old", old), ("free", free), ("integ", integ)):
        vv = v[s][np.isfinite(v[s])]
        out.append((vv.mean(), vv.std(ddof=1) / np.sqrt(len(vv))))
        sums[name] += vv.mean() * 2.0
    print(f"  {p:+5.1f}   {1e5*out[0][0]:5.2f} ({1e5*out[0][1]:.2f})     {1e5*out[1][0]:5.2f} ({1e5*out[1][1]:.2f})"
          f"      {1e5*out[2][0]:5.2f} ({1e5*out[2][1]:.2f})      {out[1][0]/out[0][0]:5.2f}     {out[2][0]/out[0][0]:5.2f}")
for name, sm in sums.items():
    print(f"  {name:6s}: sum {1e5*sm:6.2f} % ms mm -> N isotropic {V*sm*1e-3/(2*SIG_ISO):.2e}")


def profile(xx, a_, xc, sx):
    h = 1.0
    return a_ * sx * np.sqrt(np.pi / 2) / 2 * (erf((xx + h - xc) / (np.sqrt(2) * sx)) - erf((xx - h - xc) / (np.sqrt(2) * sx)))


for name, v in (("old", old), ("free", free), ("integ", integ)):
    m_ = np.array([np.nanmean(v[x == p]) for p in pos])
    e_ = np.array([np.nanstd(v[x == p], ddof=1) / np.sqrt(np.isfinite(v[x == p]).sum()) for p in pos])
    pp, cov = optimize.curve_fit(profile, pos, m_, p0=(m_.max(), 0.5, 3), sigma=e_, absolute_sigma=True)
    chi2 = np.sum(((m_ - profile(pos, *pp)) / e_) ** 2)
    err = np.sqrt(np.diag(cov) * max(1, chi2 / (len(pos) - 3)))
    print(f"  {name:6s}: sigma_x {abs(pp[2]):.2f} +- {err[2]:.2f} mm (chi2/dof {chi2:.1f}/{len(pos)-3})")
