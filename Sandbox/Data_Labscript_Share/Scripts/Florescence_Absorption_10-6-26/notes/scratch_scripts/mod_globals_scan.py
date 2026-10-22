from pathlib import Path
import h5py
root = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10")
for day, a, b in (("16", 23, 37), ("17", 3, 13)):
    for n in range(a, b + 1):
        fs = sorted((root / day / f"{n:04d}").glob("*.h5"))
        if not fs:
            continue
        vals = set()
        for p in fs:
            with h5py.File(p, "r") as h5:
                g = h5["globals"].attrs
                vals.add((float(g["mod_offset"]), float(g["mod_amp"]), float(g["mod_freq"]), float(g["iris_x_mm"])))
        print(day, f"{n:04d}", len(fs), sorted(vals))
