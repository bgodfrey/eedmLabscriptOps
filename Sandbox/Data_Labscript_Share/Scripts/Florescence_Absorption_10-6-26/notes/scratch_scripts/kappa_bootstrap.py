"""Uncertainty of the matched-filter noise calibration kappa, respecting correlations.

Readings within one no-launch shot are correlated (Gaussian-filtered noise has correlation
exp(-d^2 / 4 sigma^2) at spacing d), so the reading count overstates the information. Resampling
whole shots (bootstrap over shots) gives an honest confidence interval. Also compares readings
spaced 1 sigma vs 3 sigma apart, and the single reading per shot at the typical cloud time.
"""
import json
import os
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import numpy as np

NB = Path(r"C:\Users\EEDM\labscript-suite\userlib\labscriptlib\eedmLabscriptOps\Sandbox\Data_Labscript_Share\Scripts\Florescence_Absorption_10-6-26\Matched_Filter_Shot_Averaging.ipynb")
sys.path.insert(0, str(NB.parent))
os.chdir(NB.parent)
nb = json.loads(NB.read_text(encoding="utf-8"))
ns = {}
for c in nb["cells"][:10]:                      # through Step 4
    if c["cell_type"] != "code":
        continue
    src = "".join(c["source"]).replace("%matplotlib widget", "").replace("SAVE_OUTPUTS = True", "SAVE_OUTPUTS = False")
    exec(src, ns)
    if "def show_figure" in src:
        ns["show_figure"] = lambda fig, name: None
        ns["display"] = lambda x: None
        ns["Markdown"] = lambda s: s

cse, PHASE, SIGMA, TC = ns["cse"], ns["PHASE"], ns["SIGMA_MATCH"], ns["T_C_TYPICAL"]
rng = np.random.default_rng(0)
for acq, name in (("1", "returning"), ("0", "upward")):
    per_shot = {1: [], 3: [], "tc": []}
    for d in (n for n in ns["noise_shots"] if n["acq"] == acq):
        tr = cse.shot_traces(d["shot"], d["ref"], TC[acq], SIGMA[acq])
        t, ok = tr["t"], tr["valid_match"]
        if not ok.any():
            continue
        lo, hi = t[ok][0], t[ok][-1]
        for spacing in (1, 3):
            pos = np.arange(lo, hi, spacing * SIGMA[acq])
            per_shot[spacing].append(np.array([cse.project(cse.value_at(t, tr["z_match"], x), PHASE[acq])
                                               / tr["noise_match"] for x in pos]))
        per_shot["tc"].append(np.array([cse.project(cse.value_at(t, tr["z_match"], TC[acq]), PHASE[acq])
                                        / tr["noise_match"]]))
    print(f"\n=== {name}: {len(per_shot[1])} no-launch shots, matched sigma {1e3*SIGMA[acq]:.2f} ms")
    for key, label in ((1, "readings 1 sigma apart"), (3, "readings 3 sigma apart"), ("tc", "one reading per shot at t_c")):
        shots = per_shot[key]
        kappa = np.std(np.concatenate(shots))
        boot = []
        for _ in range(4000):
            pick = rng.integers(0, len(shots), len(shots))
            boot.append(np.std(np.concatenate([shots[i] for i in pick])))
        lo, hi = np.percentile(boot, [16, 84])
        print(f"  {label:30s}: {sum(len(s) for s in shots):4d} readings, kappa = {kappa:.2f}  "
              f"(68% interval {lo:.2f}-{hi:.2f}, i.e. +/-{50*(hi-lo)/kappa:.0f}%)")
