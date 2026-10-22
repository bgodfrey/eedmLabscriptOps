import sys, warnings
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, sys.argv[1])
sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2], "52.7"]
from fluorescence_absorption_lockin import load_shot
import cloud_signal_estimation as cse
warnings.filterwarnings("ignore")
# Dead time that makes the PMT line-centre/offset-0 ratio equal the dip ratio (3.61 +- 0.43)
folder = Path(sys.argv[2])
data = []
for p in sorted(folder.glob("*.h5")):
    s = load_shot(p)
    with h5py.File(p, "r") as h5:
        off = float(h5["globals"].attrs["mod_offset"])
    tc = cse.pmt_cloud_fit(s, 6125.0)["center_s"]
    k = int(round(0.2e-3 * s["fs_pmt"]))
    ms = np.convolve(s["pmt_rate"], np.ones(k) / k, mode="same")
    data.append((off, s["t_pmt"] - tc, ms, s["fs_pmt"]))
def ratio(tau):
    tot = {}
    for off, u, ms, fs in data:
        n = ms / np.clip(1 - ms * tau, 0.05, None)
        bg = np.median(n[(np.abs(u) > 0.008) & (np.abs(u) < 0.03)])
        tot.setdefault(off, []).append(np.sum(n[np.abs(u) <= 0.008] - bg) / fs)
    a, b = sorted(tot)
    return np.mean(tot[b]) / np.mean(tot[a])
for tau in (0, 30e-9, 40e-9, 50e-9, 55e-9, 60e-9, 67e-9):
    print(f"tau {1e9*tau:4.0f} ns: PMT ratio {ratio(tau):.2f}")
