# Absorption detection: handoff, as of 2026-10-21

Cs fountain (MOT → PGC → launch), sequence `exp12_switchArmsAndLaunch_Absorption`. The goals are better absorption SNR, an understanding of the noise, and atom number, cloud size and temperature from the probe-sheet absorption dip. This folder hands the work to a new Claude session that may not have access to the data-acquisition PC.

| File | What it is |
| --- | --- |
| `README.md` | This handoff: where things stand, how the analysis works, conventions, next steps |
| [`lab_notebook.md`](lab_notebook.md) / [`lab_notebook.pdf`](lab_notebook.pdf) | Full lab notebook, exported 2026-10-21. Log is newest first. **Its "Current status", "Files" and "Open questions" sections date from 10/14–15 and are stale; this README supersedes them.** |
| [`scratch_scripts/`](scratch_scripts/) | One-off checks the notebook cites as "scratch". Not maintained; paths are hard-coded for the DAQ PC |

The live notebook is a Claude Doc: https://claude.ai/code/artifact/846313a4-da34-4556-ac8b-b91686b395b0. The user asks that every new test or analysis be logged there.

## Where things stand

All numbers are from the 10/17 line-centre iris scan (0003–0013, 88 shots) unless stated otherwise.

1. **Cloud shape: a round core plus a broad component.**
    - Core: σx = 3.05 ± 0.07 mm (iris scan, centre +0.58 mm) and σz = 3.10 mm (vertical, from transit time at iris 0 mm). Peak dip 4.0% at 0 mm.
    - Broad component, vertical: σ 6.9 mm, centred 7.5 mm above the core, ≥16% of the atoms. It is a leading tail, and it is atoms (the PMT shows it too; a no-atom shot is flat).
    - The time-integrated horizontal profile is wider (σx 3.71 ± 0.28 mm, χ²/dof 15/6). No single σx fits; report a core plus a broad component.
2. **Atom number** (isotropic σ, s = 0.2, as measured):

    | Estimator | N | Status |
    | --- | --- | --- |
    | Integrated dip | 1.45 ± 0.07 × 10⁷ (5.9 × 10⁶ cycling) | Use this; still a lower bound (far tail excluded) |
    | Free Gaussian per shot | 1.36 × 10⁷ | Cross-check; reproduces the wing excess |
    | Fixed template (1.51 ms) | 1.16 × 10⁷ | Floor; under-reads where the dip is wider than the template |

    - Referred to the 10/16 0036 lock point (÷ 0.79): N ≈ 1.8 × 10⁷ isotropic.
    - The MOT holds about 1.25 × 10⁷ by `MOT_CAL`, so the transfer is ≥ 50% (> 100% for isotropic). This tension is unresolved.
    - About 60% of the integral-vs-template gap comes from the ±6–8 mm points, from only 8 shots each.
3. **Which dip estimator to use.**
    - **Template fit:** for per-shot comparisons at a fixed iris position (SNR, alignment, detuning, drift). At 0 mm its shot-to-shot sd is half the integral's (1.2 vs 2.5 %·ms), and it is within 3–6% of the integral at |x| ≤ 2 mm.
    - **Integrated dip:** for totals (atom number, the horizontal profile). The template reads 25–30% low at ±4 mm and 2–4× low at ±6–8 mm, because the dip lasts longer there (free width 1.41 ms at 0 mm vs 2.8–3.1 ms at ±8 mm).
4. **Shot-to-shot scatter is probe-frequency jitter** from the laser lock, about 1 MHz rms (10/16 0036: line centre vs `mod_offset` 0 gives a 3.6 ± 0.4× dip ratio). The lock also drifts within a run: on 10/17 the 0 mm dip fell 4.51 → 3.34% over 22 min.
5. **Detector noise:** in the lock-in era (60 dB), the PDA100A2's own electronics dominated at f_m (3.5 of 4.1 µV/√Hz). The best 1f lock-in SNR was 11.6 ± 1.7 per shot (10/14). Since 10/16 the work uses direct DC absorption (`mod_amp` = 0) at 70 dB.
6. **Probe light at the detector:**
    - It is cut by an opening fixed in the detection path. The **lit range** is −7.5 to +8.9 mm (16.4 mm wide), the same on 10/15, 10/16 and 10/17.
    - Inside it the light is lopsided: ×9–10 brighter at −6 than at +6 mm after the 10/16 input move.
    - Gaussian beam centre and width fits gave unphysical numbers (a 30–37 mm beam centred outside the opening). Don't quote them. To align, balance the light at ±6 mm.
7. **`mod_offset`** is the DC level of the sine on `MAIN_JUMP_AMP_ao`. It sets the probe frequency only during the 70 ms jump window; the probe power changes by < 0.6%. The inferred scale is about 0.5 MHz/mV, with line centre between 7 and 12 mV (9.5 mV is used). There is no direct V→MHz calibration.
8. **Temperature:** not yet measured. It, and whether the tail is a velocity or a position tail, needs the downward pass with a longer acquisition window (notebook, 10/19 entry, "Cause, open").

## Setup and current settings

| Item | Value |
| --- | --- |
| Sequence | `controlsAndMotDiagnosticAndLaunchTest_3_10_26/exp12_switchArmsAndLaunch_Absorption.py`. **Edit only this one, never exp11.** |
| Probe | FM sheet through a 1 mm slit, double-passed (PBS + retro mirror); the return goes through a 2 mm iris on horizontal/vertical stages to a PDA100A2 |
| PDA gain | 70 dB. Dark level 52.7 mV, measured once (10/15 0010) |
| Probe power | 150 µW since 10/16 (80 µW on 10/15) |
| Modulation | `mod_amp` 0, `mod_offset` 0.0095 V (line centre) |
| Jump window | `MAIN_REL_JUMP_do` high for 70 ms: 43–113 ms into the PDA record. Cloud at about 70 ms |
| PMT | H7421 photon counting on ctr3, 10 µs bins. Dead time ≈ 60 ns; compressed on the upward pass |
| Cloud speed in the sheet | 2.30 m/s, so 1 mm = 0.43 ms |

## Data (on the DAQ PC only)

`C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\<DD>\<NNNN>\*.h5`. The exp12 data from 10/14–10/17 is about 15 GB, so it is not in the repo. To rerun the scripts elsewhere, copy the folders below.

| Day, folders | What | Settings |
| --- | --- | --- |
| 10/17 0003–0013 | Line-centre iris scan: 11 folders × 8 shots, iris 0, +6, −2, +4, −8, 0, +2, −6, +8, −4, 0 mm | 9.5 mV, 150 µW |
| 10/17 0000–0002 | Lock check, relock, offset readjust | |
| 10/16 0037 | Offset scan 7–12 mV | 150 µW |
| 10/16 0036 | Line centre vs offset 0, interleaved | 0 / 9.5 mV |
| 10/16 0023–0035 | Iris atom scan, off resonance (0034/0035 = ±10 mm, after the last 0 mm visit) | `mod_offset` 0 |
| 10/15 0017–0033 | Probe beam profile through the iris (light only) | 80 µW |
| 10/15 0010 | Dark level, probe blocked | 70 dB |
| 10/14 | Lock-in era: baselines, offset/amplitude scans, no-launch references | 60 dB |

- **h5 contents:** the PDA trace is `data/traces/TOF_florescence` (`values`). The globals include `iris_x_mm`, `mod_offset`, `mod_amp`, `mod_freq`, `pda_gain_db` and `probe_power_uW` (the last two are recorded only, not used by the sequence). The root attribute `run time` has the format `%Y%m%dT%H%M%S.%f`.

## Analysis recipe

Shared by the scripts in [`../noise_tests/`](../noise_tests/). Each script takes the helper folder (`..`, which holds `fluorescence_absorption_lockin.py` and `cloud_signal_estimation.py`) as its first argument. On the DAQ PC they run with `C:\Users\EEDM\miniconda3\python.exe -I`.

1. `load_shot(path)`. Subtract the dark level (52.7 mV), average to 10 kHz, and keep the jump window with 2 ms guards (45–111 ms).
2. Cloud time t_c = `cse.pmt_cloud_fit(shot, mod_freq or 6125)["center_s"]`. Skip shots that are saturated (PDA ≥ 9.95 V) or whose PMT fit fails. Dip-percent measures also need ≥ 0.10 V of light (outside the lit range, a few mV of error fakes 6–29% dips).
3. **Template fit:** (c0 + c1·u + c2·u²)·[1 − D·exp(−(u − u0)²/(2w²))], with u = t − t_c and four free parameters. The shape is fixed per run as the median of clean pass-1 free fits: 10/17 u0 = +0.05 ms, w = 1.51 ms; 10/16 −0.13 ms, 1.48 ms. Area = D·√(2π)·w.
4. **Integrated dip:** fit a quadratic baseline at |u| > 10 ms, then I = Σ(1 − V/baseline)·dt over |u| ≤ 10 ms. Linear and quadratic baselines agree to 0.2%.
5. **PMT normalisation per shot:** q = PMT integral / (light / mean light at that x). Scale each dip by mean(q)/q.
6. **Profile vs iris x:** a Gaussian averaged over the 2 mm iris (erf form), fitted over |x| ≤ 8 mm (the "fit range"; the "lit range" is −7.5 to +8.9 mm).
7. **Atom number:** N = v/(2σ_eff)·∫I(x)dx, with v = 2.30 m/s, σ_eff = σ0/(1 + s) and s = 0.2. σ0 = 1.405 × 10⁻⁹ cm² (isotropic) or 3.468 × 10⁻⁹ cm² (cycling). The factor 2 is the double pass.
8. **Beam-profile drift correction:** raw × mean(visits) / interp(t; 0 mm visits). Visit times come from file mtimes, which run a steady 11–12 s after `run time`. The per-point error is √(sem² + (f·y)² + (5 mV)²), with f = rms/mean of the 0 mm visits (5.8%, 5.0% and 1.8% on 10/15, 10/16 and 10/17). The 5 mV is a guess; don't quote the 10/15 low-edge blur, which depends on it.

### Main scripts (`../noise_tests/`)

| Script | Arguments | Output |
| --- | --- | --- |
| `integrated_dip_profile.py` | helper, day dir, dark mV, first, last, fig, [dropout], [--show] | Integrated vs template dip per position, profile fits, N, half-window check |
| `iris_atom_scan.py` | helper, day dir, dark mV, first, last, fig, [dropout] | Template D per shot, PMT-normalised profile, σx, N |
| `dip_vs_offset_hv.py` | helper, day dir, dark mV, first, last, fig, [dropout] | Horizontal (iris) and vertical (transit time) dip profiles, core + broad fit |
| `largest_dip_trace.py` | helper, fig, [--show] | Deepest single-shot dip in 10/16–17 (10/17 0003 shot 6, 5.80%) |
| `beam_profile_softedge_fit.py` | fig | Light vs iris with soft-edged opening fit (always calls `plt.show()`, so it blocks until the window closes) |
| `beam_profile_gauss_fit.py`, `ref_level_vs_time.py` | fig | Plain Gaussian check; 0 mm reference level vs time |
| `offset_scan_plot.py`, `detuning_ratio_check.py`, `cloud_size_plot.py` | see docstrings | 10/16 offset scan, detuning ratio, σz from traces |

Example: `python integrated_dip_profile.py .. C:\Experiments\...\2026\10\17 52.7 3 13 integrated_dip_1017.png 0.6`. Figures from each run are saved next to the scripts.

## Working conventions (from the user)

- Comment code heavily. Use short plot labels. Plot beam profiles in plain volts.
- The user likes interactive figures: new scripts take `--show`. Write the PNG either way.
- Don't present physically implausible fit parameters as results; say the model is wrong instead.
- Terms: "lit range" (−7.5 to +8.9 mm, where there is light) vs "fit range" (what a fit uses).
- After editing an `.ipynb` that is open in VS Code, tell the user to close and reopen it (disk edits break the ipympl display).
- Log each test or analysis in the lab notebook (Claude Doc above).

## Next steps (open, roughly in priority order)

- [ ] Firm up the wings: 3–4× more shots at ±6 and ±8 mm, plus an atom-free run (`MRR_TRIG_BOOL` False) through `integrated_dip_profile.py` as a null test.
- [ ] Take blocked-probe dark shots at the start and end of each session (the dark level was measured only once, 10/15 0010). Start and end each iris scan with a 0 mm visit.
- [ ] Downward pass with a longer acquisition: temperature, and velocity vs position tail.
- [ ] Log the laser lock error signal on ai3 to tie the dip scatter to frequency jitter. Calibrate `mod_offset` in MHz.
- [ ] Probe beam: get the collimated beam size and look at the spot with an IR card, iris open. Balance the light at ±6 mm.
- [ ] Optional: refit the template width from |x| ≤ 2 mm shots only (1.41 ms instead of 1.51 ms; changes D by a few %).
- [ ] Earlier offers not yet done: an ND filter on the PMT (dead time), and the `MOT_CAL` load/aperture check behind the MOT-vs-probe tension.
