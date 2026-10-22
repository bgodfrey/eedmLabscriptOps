"""Light levels on the MOT-window PDA100A (ai2), to judge whether a C12703-01 APD module could take its place.

The C12703-01 saturates at 50-60 nW, so the question is how much light the window sees: dark (shutter closed),
scattered MOT-beam light (shutter open, no atoms yet) and the atom fluorescence at coils-off.
Usage: python mot_light_levels.py <helper dir> <shot folder>
"""
import sys
import warnings
from pathlib import Path
import numpy as np

sys.path.insert(0, sys.argv[1])
from fluorescence_absorption_lockin import load_mot_trace, fit_mot_loading

warnings.filterwarnings("ignore")
GAIN_70DB = 4.75e6        # PDA100A Hi-Z gain at 70 dB (V/A)
RESP_852 = 0.5            # approximate PDA100A responsivity at 852 nm (A/W)
AREA_RATIO = (np.pi * 1.5 ** 2) / 100.0   # phi 3.0 mm APD vs 10 x 10 mm PDA100A active area

rows = []
for p in sorted(Path(sys.argv[2]).glob("*.h5")):
    mot = load_mot_trace(p)
    if mot is None:
        continue
    fit = fit_mot_loading(mot)
    rows.append((fit["dark_v"], fit["scatter_v"], fit["atoms_v_at_coils_off"], mot["v"].max()))
rows = np.array(rows)
print(f"{len(rows)} shots")
for i, name in enumerate(("dark (shutter closed)", "scatter level (shutter open, no atoms)", "atoms at coils-off", "trace maximum")):
    v = np.median(rows[:, i])
    print(f"  {name:40s}: {1e3*v:8.1f} mV")
light_v = np.median(rows[:, 3]) - np.median(rows[:, 0])   # everything above the dark level
p_pda = light_v / GAIN_70DB / RESP_852
print(f"\nlight above dark at the trace maximum: {1e9*p_pda:.1f} nW on the PDA100A "
      f"-> about {1e9*p_pda*AREA_RATIO:.2f} nW on a phi 3 mm APD at the same spot (uniform illumination assumed)")
atoms_w = np.median(rows[:, 2]) / GAIN_70DB / RESP_852
print(f"atom fluorescence: {1e9*atoms_w:.1f} nW on the PDA100A -> {1e9*atoms_w*AREA_RATIO:.2f} nW on the APD "
      f"-> {1.4e8*atoms_w*AREA_RATIO:.2f} V at 1.4e8 V/W")
