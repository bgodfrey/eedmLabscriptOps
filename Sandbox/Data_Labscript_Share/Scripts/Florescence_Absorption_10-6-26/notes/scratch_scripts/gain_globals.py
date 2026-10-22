from pathlib import Path
import h5py
root = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10")
for day, a, b in (("15", 17, 33), ("16", 23, 37), ("17", 3, 13)):
    vals = {}
    for n in range(a, b + 1):
        for p in sorted((root / day / f"{n:04d}").glob("*.h5")):
            with h5py.File(p, "r") as h5:
                g = h5["globals"].attrs
                key = (g.get("pda_gain_db"), g.get("probe_power_uW"))
            vals.setdefault(tuple(str(k) for k in key), []).append(n)
    for (gain, power), ns in vals.items():
        print(f"10/{day} {min(ns):04d}-{max(ns):04d} ({len(ns)} shots): pda_gain_db {gain}, probe_power_uW {power}")
