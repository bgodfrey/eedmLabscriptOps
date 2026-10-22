"""Does the 5 mV absolute error term in beam_profile_softedge_fit.py matter?
Runs that script's source with the 0.005 V term replaced by each value in FLOORS (figure to a scratch
file, Agg backend so plt.show() does not block), and prints the fit lines it reports.
Also prints, per run, how big each error term is at the dimmest and brightest positions.
Usage: python abs_floor_sensitivity.py <softedge script> <scratch dir>
"""
import io
import os
import sys
import contextlib
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
import matplotlib
matplotlib.use("Agg")

SRC = Path(sys.argv[1]).read_text(encoding="utf-8")
OUT = Path(sys.argv[2])
FLOORS = [0.0, 0.002, 0.005, 0.010, 0.020]
assert SRC.count("0.005 ** 2") == 1

for fl in FLOORS:
    code = SRC.replace("0.005 ** 2", f"{fl!r} ** 2")
    buf = io.StringIO()
    sys.argv = ["softedge", str(OUT / f"softedge_floor_{int(fl*1e3)}mV.png")]
    g = {"__name__": "__main__"}
    with contextlib.redirect_stdout(buf):
        exec(compile(code, "softedge", "exec"), g)
    print(f"\n===== absolute term {1e3*fl:.0f} mV =====")
    for line in buf.getvalue().splitlines():
        if any(k in line for k in ("lo ", "hi ", "s_lo", "s_hi", "chi2", "edge", "blur")):
            print(line)
    if fl == 0.005:
        # Size of each error term per position (the data are the same for every floor)
        for name, (pos, mean, err, floor) in g["data"].items():
            sem = (err ** 2 - (floor * mean) ** 2 - fl ** 2).clip(0) ** 0.5
            print(f"  {name}: x, mean mV, sem mV, drift-term mV")
            for x, m, s in zip(pos, mean, sem):
                print(f"     {x:+5.1f}  {1e3*m:7.1f}  {1e3*s:6.1f}  {1e3*floor*m:6.1f}")
