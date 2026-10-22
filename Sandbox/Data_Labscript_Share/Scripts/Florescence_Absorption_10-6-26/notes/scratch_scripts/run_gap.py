from datetime import datetime
from pathlib import Path
import h5py
root = Path(r"C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10")
ts = {}
for name, day, n in (("10/16", "16", 23), ("10/17", "17", 3)):
    first = sorted((root / day / f"{n:04d}").glob("*.h5"))[0]
    with h5py.File(first, "r") as h5:
        ts[name] = datetime.strptime(h5.attrs["run time"], "%Y%m%dT%H%M%S.%f")
    print(name, first.name, ts[name])
d = ts["10/17"] - ts["10/16"]
print("difference:", d, f"= {d.total_seconds()/3600:.2f} h")
