"""Inspect the end of the MOT_fluorescence trace: coils off, PGC jump, MOT shutter close.

Usage: python mot_shutter_close.py <out dir for PNG> <shot files...>
Prints the commanded event times (time markers and the digital/analog outputs that change in
the window) and the PDA level in short windows, then saves one zoomed plot per file.
"""
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(sys.argv[1])


def decode(value):
    return value.decode() if isinstance(value, (bytes, np.bytes_)) else str(value)


for path in map(Path, sys.argv[2:]):
    with h5py.File(path, "r") as h5:
        tr = h5["data/traces/MOT_fluorescence"][()]
        markers = {decode(r["label"]): float(r["time"]) for r in h5["time_markers"][()]}
        g = h5["globals"].attrs
        t_redshift = float(g["t_redshift"]) if "t_redshift" in g else None
    t, v = tr["t"].astype(float), tr["values"].astype(float)
    t0 = markers["Loading MOT"]
    coils = markers["MOT coils off"] - t0
    pgc = markers["PGC on MOT arm"] - t0
    close = markers["Close MOT shutter -- dark time"] - t0
    switch = markers["Switch to MRR arm"] - t0
    tt = t - t0
    print(f"\n=== {path.name}  (t_redshift = {t_redshift})")
    print(f"  commanded after shutter-open command: coils off {coils*1e3:.1f} ms, PGC {pgc*1e3:.1f} ms, "
          f"shutter close {close*1e3:.1f} ms, LCR switch to MRR {switch*1e3:.1f} ms, trace ends {tt[-1]*1e3:.1f} ms")
    # Level in 1 ms bins from 5 ms before coils off to the end.
    edges = np.arange(coils - 0.005, tt[-1], 0.001)
    print("  t (ms rel. coils off) : mean V")
    for a in edges:
        m = (tt >= a) & (tt < a + 0.001)
        if m.any():
            print(f"    {1e3*(a - coils):6.1f} : {v[m].mean():.4f}")
    # Zoomed plot at full rate with event lines.
    fig, ax = plt.subplots(figsize=(11, 5), layout="constrained")
    show = tt >= coils - 0.010
    ax.plot(1e3 * (tt[show] - coils), 1e3 * v[show], lw=0.5)
    for x, label, c in ((0, "coils off", "C1"), (pgc - coils, "PGC jump", "C2"),
                        (close - coils, "shutter close (commanded)", "C3"),
                        (switch - coils, "LCR switch to MRR arm", "C4")):
        ax.axvline(1e3 * x, color=c, ls="--", label=label)
    ax.set_xlabel("Time after MOT coils off (ms)")
    ax.set_ylabel("MOT PDA (mV)")
    ax.set_title(path.name)
    ax.legend()
    fig.savefig(OUT / f"{path.stem}_shutter_close.png", dpi=90)
    plt.close(fig)
