# Absorption lock-in SNR lab notebook

Oct 14, 2026 · @Ben

## Current status

Best returning-cloud result so far: **1f lock-in SNR 11.6 ± 1.7 per shot** (exp12, 2026-10-14 shots 0049–0058; the morning baseline 0005–0014 gave 10.3 ± 1.5). That uses a background taken >20 ms from the cloud; the original ±8 ms window gave 8.8 ± 0.9. Either way it is about 1.6–2× the exp11 returning shots, mostly from more atoms. Without atoms the lock-in noise is 0.258 mV at 250 kS/s (0.235 mV at 500 kS/s), the MOT atom signal is stable to 4%, and the shot-to-shot scatter enters after the MOT, and pure noise alone reaches an "SNR" of about 2.4.

| Setting | Current value |
| --- | --- |
| Script | `exp12_switchArmsAndLaunch_Absorption.py` |
| NI6363 AI rate | 250 kS/s per channel (500 kS/s until ai2 was added) |
| PDA100A2 gain | 60 dB (recorded in `pda_gain_db`) |
| Probe power | 80–90 µW before a 1 mm slit, about 8 µW in the sheet (recorded in `probe_power_uW`) |
| Modulation | `mod_amp` 0.01 V, `mod_freq` 6125 Hz, `mod_offset` 0.006 V |
| PDA DC level (returning) | 5.7–8.5 V |

- **What limits the noise:** the PDA100A2's own electronics at 60 dB: 3.5 of the 4.1 µV/√Hz next to f\_m (74% of the noise power, measured with the detector blocked on 2026-10-15). Shot noise and laser intensity noise are about 13% each (RIN about −130 dBc/Hz). No-launch shots confirm there's no repeatable pickup or background to subtract.
- **Shots needed per scan setting:** 5–6 to resolve a 20% change, about 8 for 15% (per-shot SNR scatter is 10–15%).
- **Next (plan of 2026-10-15, after the slit finding):** (1) centre the beam on the slit and profile it before the slit in both directions; (2) blocked-detector shots at 50 dB, about 5; (3) about 3× the light in the sheet at 50 dB (more power, or a cylindrical telescope), bracketed by 1× at 60 dB, about 10 shots each: predicted 1.4–1.6× SNR; (4) `mod_freq` 6125 vs 10500 Hz at 50 dB, low priority; (5) open the probe shutter about 300 ms early on return shots only; (6) the per-shot amplitude scatter (16–34%, most likely laser detuning following the thermal swing) is about 10× the noise (1/49 ≈ 2%): 1f+2f amplitude from existing data, and log the laser lock error signal on ai3 (2026-10-16: probe-frequency jitter confirmed, about 1 MHz rms; see 0036).

## Setup

A Cs cloud is launched with the moving-retroreflector (MRR) stage. A frequency-modulated probe sheet detects it twice: on the way up (`_0` file, recording starts 16 ms after drop) and on the way back (`_1` file, recording starts 484 ms after drop). Each recording lasts 200 ms.

- **Probe:** horizontal sheet about 1 mm tall × 18 mm wide, retroreflected. The returning beam is focused onto the PDA100A2 with a lens. A mechanical shutter (`PROBE_SHUTTER_do`) opens it about 40 ms before the modulation starts.
- **Detectors:**
  - PDA100A2 photodiode → NI6363 `ai0` (trace `TOF_florescence`)
  - modulation readback → NI6363 `ai1` (trace `RelJumpModulation`)
  - PMT → photon counter `ctr3`, 10 µs bins (trace `photon_counts`)
- **Modulation:** `MAIN_JUMP_AMP_ao` (NI6738 `ao3`) plays a 100 ms sine, amplitude `mod_amp`, frequency `mod_freq`, DC offset `mod_offset`, generated at 100 kS/s. `MAIN_REL_JUMP_do` enables the jump for the first 70 ms, which is the analyzed window.
- **Lock-in analysis:** complex demodulation against the phase fitted to the recorded reference, then a triangular low-pass filter (two boxcars, `AVERAGING_PERIODS = 2`, reported as a 0.25 ms boxcar). The cloud window is ±8 ms around the PMT fluorescence peak; the background is the rest of the enabled window.
- **SNR definition used throughout:** peak of R = |Z − background| inside the cloud window, divided by the complex RMS of Z over the background.

Test settings live in the runmanager globals group `absorption_tests`, enter values in SI units (the Units column is only a label):

| Global | Meaning | Used by script? |
| --- | --- | --- |
| `mod_amp` | Sine amplitude (V) | Yes |
| `mod_freq` | Modulation frequency (Hz) | Yes |
| `mod_offset` | Sine DC offset (V), sets mean detuning | Yes |
| `pda_gain_db` | PDA100A2 gain knob (dB) | No, record only |
| `probe_power_uW` | Measured probe power (µW) | No, record only |
| `MRR_TRIG_BOOL` | Stage trigger on/off; `False` = no launch, everything else identical | Yes |

## Log

Newest entries first. Shot numbers are `YYYY-MM-DD_NNNN`; data lives under `C:\Experiments\eedmLabscriptOps\<script name>\YYYY\MM\DD\NNNN\`.

### 2026-10-20 · Time-integrated dip replaces the fixed template: profile wider, N ×1.25

**Method:** For each shot, a quadratic baseline is fitted at |t − t\_c| > 10 ms, and 1 − V/baseline is summed over |t − t\_c| ≤ 10 ms. No dip shape is assumed. The same PMT normalisation as iris\_atom\_scan.py is applied. Script `noise_tests\integrated_dip_profile.py`, figures `integrated_dip_1017.png` and `integrated_dip_1016.png`. The old estimate for comparison is D × √(2π) × 1.51 ms.

**10/17 (line centre):**

| Iris (mm) | Integrated (%·ms) | Old (%·ms) | Ratio |
| --- | --- | --- | --- |
| −8 | 2.7 ± 0.8 | 0.9 ± 0.2 | 2.9 |
| −6 | 1.3 ± 1.7 (median 2.6) | 1.6 ± 0.4 | 0.8 (one shot at −9.7) |
| −4 | 6.4 ± 0.7 | 4.9 ± 0.2 | 1.30 |
| −2 | 10.8 ± 1.2 | 10.5 ± 0.4 | 1.03 |
| 0 | 15.7 ± 0.5 | 14.9 ± 0.2 | 1.06 |
| +2 | 13.9 ± 1.2 | 13.1 ± 0.6 | 1.06 |
| +4 | 9.8 ± 0.7 | 7.9 ± 0.2 | 1.25 |
| +6 | 6.4 ± 0.9 | 3.5 ± 0.3 | 1.8 |
| +8 | 6.4 ± 1.7 | 1.6 ± 0.5 | 4.1 |

- Near the centre the template captures almost everything (an extra 3–6%, from the tail within ±10 ms). At ±4 mm the template misses 25–30%. At ±6–8 mm the integrated dip is 2–4× the template value but noisy. The one exception, the −6 mm mean, comes from a single bad shot.
- Horizontal width: the iris-averaged Gaussian fitted to the integrated dip gives σx = **3.71 ± 0.28 mm** (χ²/dof 15/6), vs 3.05 ± 0.07 mm from D. The model-free rms over the lit points is 3.88 mm. So the earlier "round cloud" comparison (σx 3.05 vs vertical core σz 3.10) compared core widths only. The column-integrated horizontal profile is wider, consistent with a broad component in x as well as in z.
- Atom number (s = 0.2, as measured): **1.45 × 10⁷ isotropic / 5.9 × 10⁶ cycling**, vs 1.16 × 10⁷ / 4.7 × 10⁶ from D (×1.25). The sum of medians and the Gaussian area agree (1.46 and 1.43 × 10⁷). Integration half-windows of 8 / 10 / 12 ms give 1.34 / 1.45 / 1.52 × 10⁷, and the far leading tail is still excluded, so this is a lower bound. With the same 0.79 lock-point factor as before, N ≈ 1.8 × 10⁷ isotropic / 7.4 × 10⁶ cycling, which sharpens the tension with the MOT's \~1.25 × 10⁷.
- Cost: the integral picks up slow baseline errors. The shot-to-shot sd at 0 mm is 2.5 %·ms vs 1.2 %·ms for the template. D remains the better per-shot SNR estimator; the integral is the right total.

**10/16 (mod\_offset 0):** too noisy to add anything. The sems are comparable to the signals, so σx = 5.5 ± 1.1 mm is not meaningful. The ratio of the sums is 1.39, the same direction.

**Which to believe (10/20): three estimates of the area.** All three use the same 10/17 shots and PMT normalisation (scratch area\_three\_ways.py): the template (D × √(2π) × 1.51 ms), a free Gaussian per shot with no cut (D × √(2π) × its own width), and the model-free integral.

|  | Template | Free Gaussian | Integral |
| --- | --- | --- | --- |
| N (isotropic, as measured) | 1.16 × 10⁷ | 1.36 × 10⁷ | 1.45 × 10⁷ |
| σx of Gaussian fit (χ²/dof) | 3.05 (17/6) | 3.27 (28/6) | 3.71 (16/6) |
| ±4 mm, ÷ template | 1 | 1.08–1.09 | 1.25–1.30 |
| +6 / +8 mm (%·ms) | 3.5 / 1.6 | 6.8 / 6.7 | 6.4 / 6.4 |

- The free Gaussian reproduces the wing excess, so the excess isn't an artefact of the raw sum. The template is low by construction wherever the dip is wider than 1.51 ms.
- The integral sits about 7% above the free Gaussian. That is the non-Gaussian leading tail within ±10 ms, which the vertical profile already showed is atoms.
- Statistically, the integral's total is 1.45 ± 0.07 × 10⁷ against 1.16 ± 0.02 × 10⁷, a gap of about 4σ.
- Soft spot: ±6–8 mm supply about 60% of the gap, and +8 mm alone about a third, from 8 shots at the edge of the lit range (0.44 V of light).
- Independent consistency check: the vertical profile at x = 0 has ≥16% in a 6.9 mm broad component. If that component is as wide in x, then within the lit range it is about 1.3× the core, close to the 1.25 seen. This assumes the halo is isotropic.
- Conclusion: for N, use the integral (1.45 × 10⁷, still a lower bound); the template value is a floor. For the width, no single σx describes the profile: all three Gaussian fits have χ²/dof 2.5–4.6. 3.05 mm describes the core and about 3.7 mm the whole column. Report a core plus a broad component rather than one number. To firm up the wings, take 3–4× more shots at ±6 and ±8 mm, plus an atom-free run through the same analysis as a null test.

### 2026-10-19 · Dip vs horizontal and vertical offset: the cloud is skewed vertically

**Data:** 10/17 0003–0013 (line centre, 9.5 mV; 88 shots, none below 0.6 × the median PMT). Script `noise_tests\dip_vs_offset_hv.py`, figure `noise_tests\dip_vs_offset_hv_1017.png` (both panels in dip %, same mm scale).

**There is no vertical stage scan.** Only `iris_x_mm` was ever varied (10/14–10/17); the vertical stage stayed put. The vertical profile comes from time of flight: the dip versus time at iris 0 mm (24 shots), with z = v·(t\_c − t) and v = 2.30 m/s, so + is the top of the cloud, which arrives first. Each trace is divided by a quadratic baseline fitted more than 6 ms from the cloud and scaled by the shot's PMT factor.

|  | Horizontal (iris scan) | Vertical (time of flight) |
| --- | --- | --- |
| Peak dip | 4.01 ± 0.10% (fit); 3.93 ± 0.06% at 0 mm | 4.1% (data) |
| Shape | Gaussian, σx 3.05 ± 0.07 mm (χ²/dof 16.7/6) | skewed: σ 2.5 mm below the peak, 4.4 mm above; rms 3.6 mm |
| Centre | +0.58 ± 0.08 mm | maximum about 1 mm below the PMT-fit centre |

The panels agree at the centre (4.0% at z = 0 against 3.9% at the 0 mm iris point).

*Update, later 10/19 (supersedes the two-sided fit and "What it changes" below):* the two-sided Gaussian was the wrong shape: it undershot the peak and missed both the shoulder and the tail. The tail is also longer than first seen; it runs past the ±6 ms baseline window, which was eating part of it. The script now fits the baseline only beyond ±10 ms, uses 1 mm bins (the sheet height) and draws both panels on the same mm scale. Model comparison over z ≤ 23 mm (χ² with per-bin sem, which ignores the bin-to-bin correlation, so compare between models only): Gaussian 1224/93, two-sided Gaussian 688/92, exponentially modified Gaussian 286/92, two Gaussians 202/90 (`vertical_model_compare.py`, scratch).

**Result (two Gaussians):** a core with σz = 3.10 mm at −0.5 mm, plus a broad component (σ 6.9 mm, centred +7.5 mm above the core) holding 16% of the atoms. The core matches σx = 3.05 mm, so **the core is round** (aspect 1.02). The tail is still 3% of the peak at +20 mm. It runs into the baseline window (the record starts 25 ms before the cloud), so 16% is a lower bound.

**The tail is atoms.** The dead-time-corrected PMT (τ 60 ns) shows the same early tail out to −16 ms (z ≈ +37 mm). A no-atom shot (10/17 0000 #7, lock lost) is flat in the PMT to ±0.2% of the peak (`pmt_vs_pda_profile.py`, `no_atom_baseline_check.py`, scratch). A probe-power hump would move the PMT background opposite to the PDA dip; here both rise.

**Consequences:** the pencil atom number uses one fixed-shape Gaussian per shot, so it misses the broad part and is low by at least about 16%. Use the core σz = 3.1 mm for the temperature. The axis label now says "transit time": what was called time of flight here is the transit through the sheet, z = v·(t\_c − t).

**The vertical profile has a leading tail.** The dip is 0.3–0.4% at z = +9 mm (7–9% of the peak) and zero at −9 mm, and the lower edge is steeper than a Gaussian. The tail is the same with quadratic or linear baselines from 6–12, 6–15 or 6+ ms (`vertical_tail_check.py`, scratch). The PMT envelope, which has no light baseline, shows the same early excess (0.29 of its amplitude at −4 ms against 0.09 at +4 ms). So these are atoms, not a baseline artifact. A slow detector would give a trailing tail, not a leading one. This is the "upper side high" residual noted on 10/16, so it is not new.

**What it changes:** the earlier σz of 3.27–3.29 mm and aspect ratio of 1.03–1.06 describe the core. The core is round, but extra atoms ride above it. The atom number is unchanged (it uses the time integral). A temperature from σz should use the core, or the lower side.

**Cause, open:** either faster atoms from the launch (a velocity tail) or a MOT with an upward tail (a position tail). The downward pass separates them. A velocity tail moves to the late side and broadens about 6× (556 against 87 ms of flight); a position tail stays early with the same width. The shots record only 0.2 s (PDA) and 0.15 s (PMT), so this needs a longer acquisition window over the 556 ms crossing.

The horizontal wings also sit above the Gaussian at ±8 mm (0.23 and 0.41% against 0.09 and 0.25%), which fits a broad component, but those points are close to the opening's edges.

**Input move (from the user):** the input was moved on 10/16, not the retro mirror. The opening stayed fixed (see the soft-edge entry), so it is not an aperture on the input optics.

### 2026-10-17 · Probe light profile, soft-edged fit: a fixed 16.4 mm opening; the beam centre is not measurable

**Data:** the light through the 2 mm iris from 10/15 0017–0033 (before realignment), 10/16 0023–0035 and 10/17 0003–0013 (after), each drift-corrected by its 0 mm visits. Script `noise_tests\beam_profile_softedge_fit.py`, figure `noise_tests\beam_profile_softedge.png`.

**Model:** a Gaussian beam times an opening whose edges are erf steps with their own blur, averaged over the iris. The hard-edged model (`iris_profile_fit.py`) could not fit 10/16 (χ²/dof 208/6): the iris saw 0.16 V at −10 mm, beyond the − edge.

|  | 10/15 | 10/16 + 10/17 |
| --- | --- | --- |
| − edge | −7.61 ± 0.45 mm, blur 1.1 ± 0.6 mm | −7.54 ± 0.09 mm, blur 0.96 ± 0.04 mm |
| + edge | +8.64 ± 0.28 mm, blur 0.25 ± 0.15 mm | +9.08 ± 0.16 mm, blur 0.2 ± 0.3 mm |
| χ²/dof | 14/7 | 29/12 |

**The opening did not move** while the beam crossed from the + side to the − side: −7.5 to +8.9 mm, 16.4 mm wide, centred +0.65 mm (the cloud is at +0.6 mm). Its − edge is soft (about 1 mm) and its + edge sharp (about 0.2 mm). Knife-edge diffraction blurs over roughly √(λL/2), which would put the + edge about 0.1 m before the iris (the LCP4S or cage parts in the detection arm) and the − edge 1–2 m away. That is suggestive only; divergence also blurs.

**The beam centre and width are not constrained.** Each date shows one flank, falling more slowly than a Gaussian of the expected size. The fits put the centre far behind an edge (+9 ± 4 mm on 10/15, −22 ± 5 mm after) with w 16–21 mm. One beam of a single width, moved between dates, fits worse (Δχ² = +33 for 5 more dof) and needs a 35 mm move, so those parameters are not physical.

**Centring without a beam model:** for a symmetric beam, ln(I(−6 mm)/I(+6 mm)) is proportional to the beam's offset. It was −1.7 on 10/15 and +2.2 after the realignment, so the balance point sits about 43% of the way from the old adjustment to the new one. Steer until the light at −6 and +6 mm matches.

**Supersedes** the hard-edged 10/15 numbers (opening −9.5 to +9.0 mm, centre +5.2 mm, w 11.0 mm). The −8 mm atom-scan points sit on the soft − edge, which fits their slight excess.

*Update 10/19:* the user moved the input on 10/16, not the retro mirror. The opening did not move with it, so it is not an aperture on the input optics.

*Update 10/19:* `beam_profile_softedge.png` is now plotted in volts (drift-corrected), each run with its own fitted curve, instead of normalised to the fit maximum. Size of the drift the 0 mm visits correct, over the \~25 min runs:

| Run | Range of the 0 mm level | rms | Error floor per point |
| --- | --- | --- | --- |
| 10/15 | 14% | 5.8% | 5.8% |
| 10/16 | 12% | 5.0% | 5.0% |
| 10/17 | 4% | 1.8% | 1.8% |

Shot-to-shot scatter within a folder is 1–7%. Drift between visits (about 10 min apart) is not followed, so each point carries the rms as an error floor. This is multiplicative probe-power drift. The additive dark offset (52.7 mV) is not corrected by the 0 mm division, but it only matters at the dark-side points (a few mV on 10/15) and is covered by the 5 mV error term.

**Reference level vs time** (`noise_tests\ref_level_vs_time.py`, figure `noise_tests\ref_level_vs_time.png`): the 0 mm mean versus minutes since each run's first shot (from labscript's `run time`), with error bars equal to the sd of the shot medians.

- **10/15:** 1.29 → 1.12 → 1.21 → 1.30 V.
- **10/16:** 2.68 → 3.03 → 2.89 V.
- **10/17:** 2.74 → 2.80 → 2.67 V.

The within-visit sd is 1–7%. The steps between visits (up to 13% in 10–15 min on 10/15 and 10/16) are several times larger than that, so the drift is real. It is also not monotonic: 10/15 fell and then recovered. Linear interpolation between visits 10–15 min apart only approximates it, which is why the visit rms is carried as an error floor. 10/17 was the most stable (±2%). The levels differ between days because of probe power (80 µW on 10/15, 150 µW after) and alignment.

**Plain Gaussian check** (`noise_tests\beam_profile_gauss_fit.py`, figure `noise_tests\beam_profile_gauss.png`). A Gaussian A·exp(−2(x − x₀)²/w²), with no edges, is fitted only to the points with |x| ≤ 6 mm and then extrapolated. The data are plotted divided by that Gaussian.

|  | Centre x₀ | w (1/e²) | χ²/dof |
| --- | --- | --- | --- |
| 10/15, \|x\| ≤ 6 | +8.5 ± 1.8 mm | 15.5 ± 1.6 mm | 2.8/4 |
| 10/15, \|x\| ≤ 4 | +7.7 ± 3.2 mm | 14.4 ± 3.0 mm | 2.0/2 |
| 10/16–17, \|x\| ≤ 6 | −16.0 ± 2.6 mm | 18.5 ± 1.5 mm | 30.9/10 |
| 10/16–17, \|x\| ≤ 4 | −19 ± 10 mm | 20 ± 5 mm | 26.4/6 |

Inside ±6 mm the data / Gaussian ratio stays within 0.95–1.10, so the light is Gaussian there. The ratio shows where the beam is cut, and it is the same before and after the input move:

- **− side:** 1.0 at −6 mm, 0.35–0.42 at −8 mm, 0.01–0.13 at −10 mm. 10/15's 0.09–0.13 beyond the edge is 6–10 mV, within the dark-level uncertainty.
- **+ side:** 0.83–1.03 at +8 mm, 0.01–0.09 at +10 mm.

This agrees with the soft-edge fit (−7.5 and +8.9 mm), without any edge model. The width is consistent between dates (w ≈ 15–19 mm). The centre swings from about +8 mm to about −16 mm, which suggests the input move shifted the beam by roughly 25 mm at the detector. Only one flank is seen on each date, so x₀ and w trade off; the 4 mm window gives the same answer with larger errors.

**How far off is the centre, if the beam is Gaussian?** (`gauss_centre_profile.py`, scratch; 10/16–17, |x| ≤ 6 mm.) I fixed the centre at a series of values, with the width and amplitudes free. The free optimum is χ² = 31 at x₀ = −16 mm.

| Fixed centre | χ² |
| --- | --- |
| −12 mm | 43 |
| −10 mm | 69 |
| −8 mm | 132 |
| −6 mm | 286 |
| −20 mm | 36 |
| −25 mm | 46 |

Even after inflating the errors by χ²/dof ≈ 3, a centre inside about −10 mm is excluded. Each day on its own gives the same answer (a log-quadratic gives −16 mm and −15 mm). The model-free fact behind it is that the light at −6 mm is 9–10× the light at +6 mm, and ln(light) is only weakly curved, which forces a wide beam (w ≈ 18 mm) centred far out.

The weak link is the Gaussian assumption. A 1/e² diameter of 30–37 mm is large, and the \~25 mm shift from the input move is large. A position-dependent loss (for example, polarization varying across the beam and changing the PBS reflection into the detection arm, or the PDA not moving with the iris) would make the same gradient without the beam being off-centre. Open questions for the user: the beam size out of the collimator, whether the PDA moves with the iris, and how large the input move was. A direct look with an IR card, iris open, would settle it.

**Correction (10/19, after discussion with the user):** the Gaussian beam parameters above are not physical and should not be quoted. That covers x₀ ≈ +8 / −16 mm and w ≈ 15–19 mm from the plain fit, x₀ ≈ +9 / −21 mm from the soft-edge fits, and the \~25 mm beam shift. A 30–37 mm beam centred 16 mm outside a 16 mm opening, moving 25 mm after an input nudge, means the model is wrong. The light inside the opening is not a clean clipped Gaussian; its slow 14 mm fall-off is the wrong shape. The earlier suggestion that centring could put a few times more light on the cloud depended on those parameters and is withdrawn.

**What stands (model-free):**

- The cut-offs are at ≈ −7.5 and +9 mm, the same on all three days.
- The light is lopsided inside the opening: ×5 brighter at +6 than at −6 on 10/15, and ×9–10 brighter at −6 than at +6 on 10/16–17.
- The input move flipped which side is bright.

For alignment, balance the light at −6 and +6 mm. Next: get the collimated beam size and refit with w fixed to it, and look at the spot with an IR card, iris open.

**The 5 mV term in the beam-profile errors (10/19).** Each position's error is √(sem² + (f·y)² + (5 mV)²). The 5 mV was carried over from iris\_profile\_fit.py, where it is labelled "dark-level drift". It was never measured: the only dark run (10/15 0010) spans 52.2–53.8 mV, and there is no dark run on 10/16 or 10/17. The term only matters where there is almost no light (±10–12 mm), where the sem is 0.7–4 mV. At the bright positions the sem and the drift term are tens to hundreds of mV.

To check the sensitivity, I refitted the soft-edge model with the term set to 0, 2, 5, 10 and 20 mV (scratch script abs\_floor\_sensitivity.py):

| Abs. term | 10/16–17 (B): lo / hi / s\_lo, χ² | 10/15 (A): lo / hi / s\_lo, χ² |
| --- | --- | --- |
| 0 mV | −7.54 / 9.07 / 0.96, 29.6/12 | −7.50 / 8.63 / 1.48, 151/7 |
| 2 mV | −7.54 / 9.07 / 0.96, 29.6/12 | −7.52 / 8.63 / 1.40, 37.5/7 |
| 5 mV | −7.54 / 9.08 / 0.96, 29.3/12 | −7.61 / 8.64 / 1.13, 14.4/7 |
| 10 mV | −7.54 / 9.08 / 0.96, 28.3/12 | −7.72 / 8.68 / 0.08, 6.8/7 |
| 20 mV | −7.54 / 9.07 / 0.95, 25.5/12 | −7.73 / 8.69 / 0.06, 3.6/7 |

The 10/16–17 edges and blur do not depend on the term. On 10/15 the edge positions move by at most 0.2 mm, but the low-edge blur jumps from 1.1–1.5 mm to <0.1 mm between 5 and 10 mV. The beam is dim on that side on 10/15, so the blur rests on a few mV-level points. **Don't quote the 10/15 s\_lo.** The 10/15 points beyond the edges also read +6 to +11 mV above the dark level. That could be the dark level having shifted since 0010, or stray light. The fix is to measure the dark level at the start and end of each scan, or to fit a constant background, rather than hiding it in the error.

**Drift-correction timing check (10/19).** Each shot is linearly interpolated, in time, between the 0 mm visit before it and the one after it. It is not scaled by the last visit alone. The relative floor f is a single number per run (10/15 5.8%, 10/16 5.0%, 10/17 1.8%) and is not a fit parameter. The scripts take shot times from file mtimes. On all three days the mtime is a constant 11.2–12.3 s after the labscript 'run time', so the mtimes are safe to use. Caveat: on 10/16 the −10 and +10 mm folders (0034, 0035) came after the last visit (0033), so they are held at 0033's level rather than interpolated. The −10 mm point (155 mV) is one of the few that pin the low-edge blur. Next time, end each scan with a 0 mm visit (scratch script visit\_bracketing.py).

**Deepest single-shot dip, 10/16–17 scans (10/19).** Ranked by a model-free measure: 1 − light/baseline, with the quadratic baseline fitted more than 10 ms from the PMT cloud centre, averaged in 1 mm (0.43 ms) bins, no PMT normalisation. The largest is **10/17 0003 shot 6 at iris 0 mm: 5.80%** (fixed-shape D 5.58%). It dips from 2.74 V by about 160 mV. The top 10 are all from 10/17 at 0 or +2 mm. The three deepest are from the first 0 mm visit, which has the highest PMT signal of the run. The 10/16 maximum is 2.73% (0033 shot 3). At 0 mm the median peak dip is 1.43% on 10/16 and 4.17% on 10/17, and the median PMT fluorescence is 0.67 vs 1.33 (relative to the median over both scans). So 10/17 has about 2× the PMT signal but 2.9× the dip: more atoms is most of the difference, but not all of it. Pitfall: the 10/16 +10 mm shots (10–27 mV of light, outside the lit range) give fake 6–29% "dips". The script now requires at least 100 mV of light. Script: noise\_tests/largest\_dip\_trace.py (figure largest\_dip\_trace.png; --show opens a window).

**Correction (10/19):** "more atoms is most of the difference" above is wrong. The 10/16 scan ran at `mod_offset` 0 (off resonance) and the 10/17 scan at 9.5 mV (line centre). From 0036 that detuning alone is worth 3.6 ± 0.4× in dip and about 2.3× in PMT signal (dead time). The observed 2.9× dip and 2.0× PMT signal therefore mean 10/17 had about the same atom number or slightly fewer (about 0.8–0.9×), not more.

**`mod_offset` and probe power (10/19).** `mod_offset` is the DC level of the 100 ms sine on `MAIN_JUMP_AMP_ao` (ao3). It sets the probe frequency only while `MAIN_REL_JUMP_do` is high: 70 ms, 43–113 ms into the PDA record, with the cloud near 70 ms. It does not change the probe power. In 0036 and 0037 the light inside the jump window divided by the light after it is 0.994–1.000 at every offset from 0 to 12 mV, with the same ±0.3–0.9% shot scatter (scratch light\_vs\_mod\_offset.py). There is no direct volts-to-MHz calibration. The 0036 ratio implies roughly 0.5 MHz per mV for a Lorentzian, and line centre is somewhere between 7 and 12 mV (9.5 mV is used).

**The dip's duration grows with horizontal offset (10/17, 10/20).** The fixed dip width (1.51 ms) is an unweighted median of pass-1 free widths over the 72 of 88 shots that pass the 8σ cut. Broken down by iris position (scratch shape\_width\_weighting.py):

| Iris (mm) | Clean / shots | Median D | Median free width (ms) | Typical error (ms) |
| --- | --- | --- | --- | --- |
| −8 | 2/8 | 0.25% | 2.77 | 0.21 |
| −6 | 6/8 | 0.42% | 2.24 | 0.07 |
| −4 | 8/8 | 1.22% | 1.61 | 0.04 |
| −2 | 8/8 | 2.94% | 1.45 | 0.02 |
| 0 | 24/24 | 4.02% | 1.41 | 0.01 |
| +2 | 8/8 | 3.83% | 1.43 | 0.02 |
| +4 | 8/8 | 2.23% | 1.67 | 0.04 |
| +6 | 6/8 | 0.94% | 2.43 | 0.14 |
| +8 | 2/8 | 0.76% | 3.05 | 0.32 |

The widening is far larger than the errors from ±4 mm outward, where every shot passes the cut, so it isn't noise. At ±8 only 2 of 8 shots pass, which could favour wide dips. The dip lasts about twice as long at ±6–8 mm as at the centre. That is consistent with a broader component (like the vertical profile's 6.9 mm halo) dominating the horizontal wings, but this is not established. Consequences: (1) one fixed width under-reads the time-integrated dip at the wings by up to about 2×, so the horizontal profile in D (σx 3.05 mm) is narrower than the profile in integrated dip, and the atom number from D × width is low. (2) The common width is pulled up by the wings: weighted by error it is 1.44 ms, and for |x| ≤ 2 mm it is 1.41 ms, vs 1.51 ms used. Next: per-position integrated dip (model-free area) to replace D for the profile and atom number.

### 2026-10-17 · Line-centre iris scan (0003–0013): σ\_x 3.1 mm, round cloud, N ≈ 1.1–1.5 × 10⁷ (isotropic)

**Shots:** 2026-10-17 0003–0013, 11 folders × 8 launched shots, order 0, +6, −2, +4, −8, 0, +2, −6, +8, −4, 0 mm. 2 mm iris, 70 dB, `mod_offset = 0.0095`, `mod_amp = 0`, after the lock fix (0002). Script `noise_tests\iris_atom_scan.py` (new optional dropout cut on the PMT, 0.6 × median: none cut), figures `noise_tests\iris_atom_scan_1017_0003-0013.png` and `noise_tests\cloud_size_1017_0003-0013.png`; drift cross-check in scratch `scan_drift_check.py`.

| Quantity | 10/17, line centre | 10/16, offset 0 (0023–0035) |
| --- | --- | --- |
| σ\_x | 3.06 ± 0.07 mm | 3.19 ± 0.13 mm |
| Centre (stage) | +0.58 ± 0.08 mm | −0.37 ± 0.13 mm |
| Peak dip | 4.01 ± 0.10% | 1.35 ± 0.05% |
| σ\_z (0003 traces) | 3.27 ± 0.07 mm | 3.29 ± 0.09 mm (0036) |
| Fit χ²/dof | 16.6/6 | 5.6/6 |

**Drift:** the 0 mm visits fell 4.51 → 3.97 → 3.34% over 22 min (MOT 42.5 → 39.8 mV), so the lock drifted again. Dividing by the PMT leaves 4.09 → 3.75% (it under-corrects, being dead-time compressed). Raw, PMT-normalised and reference-interpolated fits all give σ\_x 3.05–3.13 mm and centre +0.57 to +0.69 mm, so the shape is robust.

**Shape:** σ\_x matches the off-resonance scan, as expected if detuning only scales the dip. σ\_z/σ\_x = 1.06 ± 0.03: the cloud is round to about 5%. The ±8 mm points sit 0.1–0.2% above the Gaussian, either faint wings or the edge artefact; the point sum and the Gaussian area agree to 1.5%.

**Cloud moved about 1 mm:** the centre is +0.6 mm against −0.4 mm on 10/16. The light profile is unchanged (the −8 mm edge value is the same to 1%), so the detection geometry did not move; most likely the MOT position shifted with the relocked detuning.

**Atom number:** as measured, 1.14 × 10⁷ (isotropic, s = 0.2) and 4.6 × 10⁶ (cycling). Per MOT mV, the 0 mm dip ran at 0.88 → 0.70 of its 0036 value during the run (about 0.79 on average), so at line centre N ≈ 1.45 × 10⁷ (isotropic) or 0.6 × 10⁷ (cycling). Today's MOT (41.5 mV) holds about 1.25 × 10⁷ by MOT\_CAL. Transfer from MOT to probe is therefore 50% or more, with the isotropic case above 100%, the same tension as before.

### 2026-10-17 · Lock offset readjusted (0002): close to the old point

**Shots:** 2026-10-17 0002, 8 launched shots at iris 0 mm, `mod_offset` 0.0095, after the user readjusted the lock's DC offset (0001's shift came from that offset). Script `noise_tests\detuning_ratio_check.py`.

|  | 0036 / 0000 | 0001 | 0002 |
| --- | --- | --- | --- |
| MOT signal | 46–47 mV | 30 mV | 43.2 mV |
| Dip at 0 mm | 5.7–6.0% | 2.5% | 4.45 ± 0.21% |
| Dip per MOT mV | 0.12%/mV | 0.084%/mV | 0.103%/mV (0.85×) |
| Dip scatter | 5–11% | 26% | 13% |
| Correlation, dip vs PMT | +0.71 | +0.98 | +0.86 |

Mostly restored. The 0.85× per-atom dip and the remaining dip–PMT correlation suggest the probe still sits about 1 MHz from the peak, on the edge of the flat top. Good enough for the line-centre scan: the shape is unaffected, and the atom number would read up to about 15% low, inside the factor-2 accuracy. The 0 mm visits in the scan track it.

### 2026-10-17 · After relocking (0001): locked, but not at the old frequency

**Shots:** 2026-10-17 0001, 8 launched shots at iris 0 mm, `mod_offset` 0.0095, after the user relocked. Script `noise_tests\detuning_ratio_check.py`.

|  | Before (0036, 0000 good shots) | After relock (0001) |
| --- | --- | --- |
| MOT signal | 46–47 mV | 30.0 mV (4% rms) |
| Dip at 0 mm | 5.7–6.0% | 2.53 ± 0.23% |
| Dip per MOT mV | 0.12%/mV | 0.084%/mV (0.70×) |
| Dip scatter | 5–11% | 26% |
| Correlation, dip vs PMT / MOT | +0.71 / +0.5 | +0.98 / +0.77 |

No shot lost the lock entirely, but the MOT is 35% smaller and the probe absorbs 30% less per MOT atom with first-order scatter, as at offset 0. Both point to a lock point shifted by roughly 2 MHz (a Lorentzian at 0.70 of the peak sits about 0.65 half-widths off), with the shot-to-shot frequency jitter now acting on the slope. The dip following the MOT (+0.77, dip changing about 4× faster) fits both being driven by the frequency. Not ready for the scan: restore the lock point until the MOT reads about 46 mV and the 0 mm dip is back near 5.7%.

### 2026-10-17 · Lock check (0000): the dropouts lead into a lost lock

**Shots:** 2026-10-17 0000, 8 launched shots at iris 0 mm, `mod_offset` 0.0095, otherwise as 0036. Script `noise_tests\detuning_ratio_check.py`.

- 5 good shots: dip 5.25–5.95%, mean 5.68%, 4.7% scatter (better than 0036's 11%).
- 2 dropouts (2.65% and 1.56%) with a normal MOT signal and the PMT down in step, as in 0037.
- Last shot: no MOT (0 mV), no PMT counts, no dip. The laser had lost lock entirely; the user confirms.

The partial dropouts in 0037 and here are the lock slipping before it fails: MOT loading over seconds still works, while the probe frequency during the 70 ms window lands off resonance. This favours the frequency explanation over launch losses. Next: relock, repeat the 8-shot check, then the line-centre scan; log the lock error on ai3 if possible so slips can be flagged per shot.

### 2026-10-16 (late) · Offset scan 7–12 mV (0037): flat top, and some shots drop out

**Shots:** 0037, 40 launched shots at iris 0 mm, `mod_offset` 7, 8.5, 9.5, 10.5 and 12 mV × 8, shuffled, otherwise as 0036. Scripts `noise_tests\offset_scan_plot.py` (figure `noise_tests\offset_scan_0036-0037.png`) and `noise_tests\detuning_ratio_check.py`; scratch `dropout_check.py`.

| Offset | Median dip | Shot scatter | Dropouts |
| --- | --- | --- | --- |
| 7.0 mV | 5.00% | 15% | 0 of 8 |
| 8.5 mV | 4.72% | 32% | 1 of 8 |
| 9.5 mV | 5.09% | 30% | 1 of 8 |
| 10.5 mV | 5.52% | 43% | 3 of 8 |
| 12.0 mV | 5.51% | 45% | 2 of 8 |
| 9.5 mV (0036, 15 min earlier) | 5.97% | 11% | 0 of 12 |

**Line centre not resolved:** the normal shots sit at 5–7% across the whole range. With the 0 mV point the half-width is about 5–6 mV, so 7 and 12 mV should be only about 15% below the peak, inside the scatter. 9.5 mV stays; a 0–19 mV range would locate the peak.

**Dropouts** (below 60% of their offset's median; a few more sit at 55–65%):

- The PMT falls with the dip, and every shot, including the off-resonance 0036 ones, lies on one dip-vs-PMT curve (bent by PMT dead time). So fewer atoms scatter: either fewer atoms reach the probe, or the probe is off resonance on those shots.
- Normal on dropout shots: MOT signal, cloud arrival time (70.3 ± 0.2 ms), dip width (1.37 vs 1.45 ms), and the recorded jump drive.
- They began between 0036 and 0037. They are deeper at larger offsets and absent at 7 mV, which points weakly to frequency excursions toward higher offset; 8 shots per offset cannot rule out launch losses. The earlier dropouts (0070, 0080) were total losses from the stage.

**Next:** log the laser lock error (ai3) or a frequency monitor to see excursions directly; check the lock and the stage when dropouts appear; keep 9.5 mV.

### 2026-10-16 · Line centre vs offset 0 (0036): dip 3.6× deeper; the shot-to-shot scatter is probe frequency

**Shots:** 0036, 24 launched shots at iris 0 mm, `mod_offset` 0 and 0.0095 interleaved (12 each, shuffled), otherwise as the scan. Script `noise_tests\detuning_ratio_check.py`; scratch `deadtime_match.py`.

|  | Offset 0 | Line centre (9.5 mV) |
| --- | --- | --- |
| Peak dip at 0 mm | 1.58 ± 0.18% | 5.69 ± 0.18% |
| Dip scatter per shot | 40% | 11% |
| PMT integral scatter | 27% | 9% |
| Correlation, dip vs PMT | +0.98 | +0.71 |
| MOT signal | 46.3 mV | 46.1 mV |
| Dip width (free fit) | 1.42 ms | 1.47 ms |
| PMT peak, measured | 6.9 Mcps | 12.2 Mcps |

**Ratio:** 3.61 ± 0.43 from the raw dips. For a Lorentzian, offset 0 sits 1.6 half-widths from line centre, about 4–5 MHz, so roughly 0.5 MHz per mV of offset. The dip width does not change, so the scan's shape stands.

**The scatter is probe frequency:** off resonance the dip and the PMT move together (+0.98) and scatter 3–4× more than at line centre, while the MOT holds to 2.5%. Launch variation would scatter both settings equally, so it is at most the 9% left at line centre. The size implies a shot-to-shot frequency jitter of roughly 0.3–0.4 half-widths, about 1 MHz rms, from the lock or the frequency jump. At line centre the dip starts to follow the MOT (+0.51).

**PMT dead time:** the PMT ratio is only 2.3 because the brighter line-centre shots lose more counts (12 Mcps measured, close to 1/τ). A non-paralysable correction matches the dip ratio at τ ≈ 60 ns (52–64 ns over the ratio's error), close to the 67 ns estimate. Without an ND filter the PMT cannot count the upward pass at line centre.

**Atom number, revised:** the scan's shape (σ\_x 3.19 mm, σ\_t 1.48 ms) scaled to the line-centre dip at 0 mm gives ∫ I dx ≈ 1720 µs·mm. That is N ≈ 1.7 × 10⁷ (isotropic σ₀, s = 0.2) or 0.7 × 10⁷ (cycling); applying the interleaved ratio to the scan instead gives 1.4 × 10⁷ (isotropic), since the offset-0 dip drifts with the lock between sessions. The MOT holds 1.4 × 10⁷, so the cross-section is probably above the isotropic value or the MOT calibration is low; either way most of the MOT reaches the probe. Accuracy a factor of about 2, as before.

**Lock-in link:** the 1f lock-in runs at 6 mV sit on a steep part of the 1f curve (+20.5 at 6 mV, −9.3 at 10 mV in 0126), so the same jitter likely drives their 16–34% amplitude scatter. Fixes: tighten the lock or the jump, or log the lock error on ai3 and correct per shot.

**Cloud size:** the 12 line-centre traces, averaged and converted to height (z = v·(t\_c − t)), give σ\_z = 3.29 ± 0.09 mm (σ\_t 1.43 ± 0.04 ms, error from single-shot fits). Against σ\_x = 3.19 ± 0.13 mm from the scan, the aspect ratio is 1.03 ± 0.05, so the cloud is round at 87 ms. The averaged trace has a small baseline residual (about ±3% of the peak, the upper side high). Script `noise_tests\cloud_size_plot.py`, figure `noise_tests\cloud_size_0023-0036.png` (both profiles normalised on one axis, and a to-scale drawing of the cloud in the sheet).

*Update 10/19:* the "upper side high" residual is atoms, not baseline. The vertical profile has a leading tail (0.3% at z = +9 mm on 10/17, robust to the baseline choice and seen by the PMT). σ\_z here describes the round core. See the 10/19 dip-vs-offset entry.

**Next:** move the beam about +5 mm onto the cloud (scan entry); keep 0.0095 for DC-dip work.

### 2026-10-16 · Iris atom scan: cloud σ\_x 3.2 mm at −0.4 mm; probe was off resonance

**Shots:** 0023–0035, 13 folders × 8 launched shots, in the order 0, +4, −6, +8, −2, 0, +6, −8, +2, −4, 0, −10, +10 mm. 2 mm iris, 70 dB, 150 µW, `mod_amp = 0`, **`mod_offset = 0`** (the plan had 0.0095). 0018–0022 are earlier tests at 0 and +4 mm. Script `noise_tests\iris_atom_scan.py`, figure `noise_tests\iris_atom_scan_0023-0035.png`; checks in scratch `atom_scan_checks.py` and `pmt_scatter_compare.py`.

| Quantity | Value |
| --- | --- |
| Cloud centre (stage) | −0.4 ± 0.13 mm |
| Horizontal size σ\_x | 3.19 ± 0.13 mm |
| Peak dip | 1.35 ± 0.05% |
| Dip timing | −0.13 ms from the PMT centre, σ\_t 1.48 ms (PMT 1.93 ms), so σ\_z 3.4 mm |
| Fit | χ²/dof 5.6/6, positions within ±8 mm |
| ∫ I dx | 400 µs·mm (sum of points 394) |
| N as measured | 3.9 × 10⁶ (isotropic σ₀, s = 0.2); 1.6 × 10⁶ (cycling) |

**Off resonance:** the 0 mm dip is 1.32%, against 3.5–5.5% through the iris at line centre on 2026-10-15 (0007–0009, same stage position, MOT 41–43 mV vs 43 mV today). The ratio, about 0.3, is mostly the detuning (`mod_offset` 0 vs 9.5 mV, about 1.5 half-widths); the higher power adds about 10%. Corrected for it, N ≈ 1.2 × 10⁷ (isotropic), consistent with 0130 (0.9 × 10⁷) and the MOT (1.4 × 10⁷). Confirmed by 0036 (ratio 3.6 ± 0.4); revised N ≈ 1.4–1.7 × 10⁷ (isotropic), see that entry.

**Shape:** σ\_x 3.2 mm matches σ\_z 3.4 mm from the dip width, so the cloud is round. At 87 ms after release, 12 µK (vertical) needs an initial σ ≈ 2.1 mm; a 1 mm MOT would mean T\_x ≈ 20 µK.

**Shot-to-shot:**

- The MOT signal is steady (2.5% rms) and the dip does not follow it (correlation −0.2).
- The dip follows the PMT fluorescence (correlation +0.74 to +0.92 within ±2 mm). The PMT integral varies 21% rms, uncorrelated with the MOT, so the variation enters after the MOT: launch transfer or probe frequency. 0130 at line centre had 11% PMT scatter.
- Dividing each shot by its PMT integral (corrected for its light level) cuts the dip scatter from 29% to 19% at the bright positions.
- What remains is a floor of 0.18% of the light per shot at every position, 6× the white-noise fit error: probe intensity noise in the dip's \~100 Hz band.

**Beam after realignment:** the light peaks at −6 mm (7.1 V), then drops to 3.2 V at −8 and 0.16 V at −10 mm. The bright part now sits against the − edge, an overshoot from +5.2 mm. The cloud sees 2.9 V, 40% of the peak (about 65% before). ±10 mm are outside the lit window and not fitted; −10 reads −0.5 ± 0.1%, at a hard edge where a few µm of pointing mimics a dip.

**Overlap (2026-10-17):** in the detector frame the cloud's shadow is centred near stage 0 (the cage axis) and the probe's brightest light at −6 mm, or further out if the − edge clips it. The sheet still covers the whole cloud (±2σ = ±6.4 mm); the light falls about 5:1 across it (5.8 V at −4 mm, 1.2 V at +4 mm). Weighting by the measured light, a full-beam dip is about 0.48 × the dip at the cloud centre; a beam of the same width centred on the cloud would give about 0.5–0.55. Centring mainly raises the light on the atoms, about 1.2–1.4× for the same power (rough), and moves the bright part off the hard edge. Iris measurements are fractional per strip and unaffected.

**Next:** set `mod_offset` back to 0.0095 and check 0 mm; move the beam about +5 mm so its peak sits on the cloud, which is also the centre of the opening.

**Redo? (2026-10-16):** no full redo. Detuning scales every position by the same factor, so the centre and σ\_x stand. One interleaved folder before touching the beam: iris at 0 mm, `mod_offset` \[0, 0.0095\] × 12 reps, shuffled (24 shots, about 5 min). The raw dip ratio converts the whole scan to line centre (about ±10%, vs about ±20% from the 10/15 shots). The PMT ratio is an independent check, and the PMT scatter at each setting tests detuning against launch variation.

### 2026-10-16 · Power swing traced to the PM fiber; plan for the iris atom scan

**Power swing:** the user traced the probe-power swings to a polarization mismatch at the polarization-maintaining fiber. It is reduced, not gone. The PBS turns the fiber's polarization drift into a power change, so the atoms see a fixed polarization at varying power. This confirms the 0127 suspicion (fiber polarization or an etalon).

**Effect on the scan:** the fractional dip (dip ÷ light above dark, same shot) does not depend on probe power at s ≈ 0.1. The remaining swing cancels shot by shot, so the light level needs no drift references.

**Plan:** absorption depth versus iris position across the upward cloud.

| Setting | Value |
| --- | --- |
| Sequence | exp12, launched, upward pass only |
| Probe | `mod_amp = 0`, `mod_offset = 0.0095` (line centre, DC dip), as in 0130 and 0007–0009 |
| Detector | PDA100A2 at 70 dB; dark 52.7 mV (0010) |
| Iris | 2 mm, vertical stage as in 0017–0033; `iris_x_mm` set for each folder |
| Positions | −8 to +8 mm in 2 mm steps (9): the whole lit window (−9.5 to +9.0 mm), strips tile with no gaps |
| Order | 0, +4, −6, +8, −2, 0, +6, −8, +2, −4, 0 (11 folders) |
| Shots | 5 per folder, 55 in all |

**Expected:** the 1.66% full-beam dip (0130) and 3.7–4.2% through the iris (0007, 0008) fit a cloud of σ ≈ 3 mm near stage 0. So expect a peak of about 4%, with the profile inside the window. Detector noise (0.9 mV per shot) is 0.05% of the light at the bright positions and 0.8% at −8 mm; the 10–14% per-shot scatter dominates, about 5% per position.

**Analysis:** each shot gets the 0130 dip fit (baseline × (1 − Gaussian) at the PMT time), integrated over time and normalised by that shot's MOT signal. The atom number needs no assumed horizontal size (v = 2.30 m/s, Δx = 2 mm, double pass):

```latex
N = \frac{v}{2\sigma_{\mathrm{abs}}} \sum_i \left( \int \alpha_i \, dt \right) \Delta x
```

A Gaussian fit to the profile gives the centre, σ\_x and the area beyond the window. At 87 ms after release, 12 µK adds 2.3 mm of expansion, so σ\_x ≈ 2.5–3.1 mm for a 1–2 mm MOT; σ\_x against the 3.0 mm vertical dip width tests isotropy.

**What it answers:** the cloud's horizontal size and position, an atom number that replaces the isotropic assumption in the 2026-10-14 cloud table, and where the cloud sits in the beam (beam centre +5.2 mm). The repeat shots also show whether the per-shot scatter fell with the swing reduced.

### 2026-10-15 (night) · Beam profile fit: an 18.5 mm opening centred on the stage, beam 5 mm off-centre in it

**Shots:** 0017–0033, no launch, 70 dB, iris about 2 mm, 4 shots per stop. 13 positions from −12 to +12 mm in shuffled order, with returns to 0 mm at 0.3, 9.4, 18.3 and 26.2 min. Script: `noise_tests\iris_profile_fit.py`; figure `noise_tests\iris_profile_fit_0017-0033.png`.

**Superseded 2026-10-17:** a soft-edged fit of the same data puts the opening at −7.6 to +8.6 mm and leaves the beam centre and width unconstrained; see the 2026-10-17 profile entry.

**Drift:** the 0 mm reference read 1291, 1120, 1211 and 1295 mV (14% range, 5.8% rms). Each shot is divided by the reference interpolated linearly in time (at most ±6%). The references are about 9 min apart and the swing is faster, so 6% is kept as an error floor.

**Model:** a Gaussian I = A exp(−2(x−x₀)²/w²), cut by an opening lo < x < hi, averaged over a 2 mm top-hat for the iris. Five free parameters; χ²/dof = 50/9, and the errors below are scaled by √5.6. Residuals are within ±6% of the peak.

| Parameter | Value |
| --- | --- |
| Beam centre x₀ | +5.2 ± 0.7 mm |
| 1/e² radius w | 11.0 ± 0.7 mm (1/e² diameter 22 mm) |
| Lower edge | −9.5 ± 0.5 mm |
| Upper edge | +9.0 mm (about ±0.3 mm; the fit's ±0.02 mm ignores the uncertain iris size) |
| Opening | 18.5 mm wide, centred at −0.2 mm |

**Supersedes** the five-point estimate (centre about +8 mm, w about 15 mm), which was distorted by drift and by having no points near the edges.

**What the opening is:** 18.5 mm is narrower than the LCP4S's SM1 bore (about 23–25 mm), so the bare adapter is not the edge unless something narrower is threaded into it. It matches the 18 mm sheet width quoted earlier, so the slit's ends are now the leading candidate; that would put the cut before the atoms. The opening is centred on the stage's zero.

**Consequences:**

- 76% of the beam's horizontal profile passes; centred, 91% would (1.2× more light).
- The + edge cuts where the beam is at 78% of its peak, the − edge at 3%. The transmitted power then changes about 7% per mm of horizontal beam motion. Centred, both edges sit at 24% and the change vanishes to first order. Horizontal pointing jitter of about 0.7 µm/√Hz at 7%/mm would match the probe's 40–100 Hz noise (4.7 × 10⁻⁵/√Hz); unproven.
- The beam is wider than the opening, so even centred the sheet is a Gaussian cut at 24% on both sides, not a full Gaussian.

### 2026-10-15 (evening) · Coarse beam profile through a \~2 mm iris: off-centre, and clipped on the + side

**Shots:** 0012–0016, no launch, 70 dB, 4 shots each at stage offsets −12, −6, 0, +6, +12 mm (`iris_x_mm`; negative = stage pushed right when facing North). Iris about 2 mm. Script: scratch `iris_centre_fit.py`, dark 52.7 mV subtracted.

| Offset (mm) | Level above dark (mV) |
| --- | --- |
| −12 | 14 ± 2 |
| −6 | 479 ± 31 |
| 0 | 1546 ± 30 |
| +6 | 2662 ± 107 |
| +12 | 10 ± 2 |

**Light level:** up to 2.7 V at 70 dB, about 9× the 1 mm iris at offset 0. The scan dip should be measurable to an SNR of about 50 per shot near the bright side.

**Off-centre:** the light rises steadily toward +, so stage 0 is not the beam centre.

**Clipped:** a Gaussian through the −6, 0 and +6 points (superseded by the 0017–0033 fit: centre +5.2 mm, w 11.0 mm; here centre about +8 mm, 1/e² radius about 15 mm; solved exactly from three points, not fitted, so roughly 14–16 mm and +7 to +9.5 mm if the between-position drift matches the within-position spread; a least-squares Gaussian through all five points gave 6 mm and +3.8 mm but cannot follow the clipped ±12 mm points and misses −6 mm by 13%) predicts about 2.4 V at +12 mm; 10 mV is measured. Something cuts the beam sharply between +6 and +12 mm, and the −12 point is also below the Gaussian (14 vs about 80 mV). The light reaching the PDA spans roughly −10 to +10 mm. Candidates: the slit's ends, the viewport, or 1-inch optics. The user suspects the LCP4S 30-to-60 mm cage plate adapter: its central SM1 bore (about 25 mm, about 23 mm with an SM1 part threaded in) fits a lit region of about −10 to +10 mm, centred near stage 0, with the beam about 8 mm to the + side. It sits after the atoms (detection arm: PBS → LCP4S → horizontal and vertical stages → PDA, all on cage rods), so it limits only the PDA's view; the atoms see the full sheet. The cage keeps the parts coaxial with each other, not with the beam: the return lands about 8 mm off the cage axis, from a retro-mirror tilt (2θL, so about 8 mrad at L ≈ 0.5 m) or from an input beam that crosses the PBS off-centre. To tell which, check whether the return retraces the input before touching the mirror. A hard edge where the beam is bright turns pointing jitter and thermal pointing drift into power changes, which makes it a candidate for the probe's 40–100 Hz noise and part of the 5-minute power swing.

**Power swing during the run:** within each 4-shot group the level fell 6–11% in under a minute (for example 2.78 → 2.49 V at +6), so the five-point shape is uncertain at the 10–20% level. The clipping is far larger than that.

**Next:** find the clipping element with an IR card or viewer, and recentre if easy (that changes the sheet, so it comes before the atom scan). Then the atom scan over the lit region: −10 to +10 in 2 mm steps plus +7, +9 and +11 for the edge, 3 shots each, shuffled with returns to 0: 0, −6, +4, +6, −4, +9, 0, +8, −8, +2, −2, 0, +7, +11, −10, +10, 0.

### 2026-10-15 (evening) · Iris scan prep: dark level, double pass confirmed, 1 mm iris too dim

**Shots:** upward only, PDA100A2 at 70 dB, `mod_amp = 0`, `mod_offset = 0.0095`, 80 µW. By the data: **0010** (4 shots, no launch) is the dark run and **0011** (4 shots, no launch) is the mirror-blocked run. Confirmed by the user (first noted as 0009 and 0010). 0009 is a single launched shot with about 6 V of light. The iris may not be centred on the beam, so stage position 0 is not necessarily the beam centre; the scan fit finds both centres. 0007 and 0008 are single launched test shots, assumed to be through the 1 mm iris. Script: scratch `dark_and_retro_block.py`, `dark_dip_noise_70db.py`, `pencil_test_dips.py`.

**Dark level at 70 dB:** 52.7 mV (52.2–53.8 mV per shot).

**Mirror blocked:** 51.6 mV (46.8–55.6 mV per shot), the same as dark within about 2 mV. Light that skips the round trip is at most about 2% of the 1 mm-iris level, so the PDA sees only the double-passed beam, as the atom-number calculation assumed.

**1 mm iris light level:** 162–187 mV above dark (about 65–70 nW), about 1% of the sheet's light; I had estimated 7–9%. Candidate causes: the iris is not centred vertically on the 1 mm sheet, the iris clips the beam on both passes (diffraction widens the returning beam), or the iris is smaller than 1 mm.

**Atom dip through the iris (0007, 0008):** 6.8–6.9 mV, 3.7–4.2% of the light. 0009 also shows 4.1% (6.0 V of light, iris setting unknown).

**Noise at 70 dB:** a Gaussian dip (σ 1.3 ms) can be measured to 0.9 mV per shot in the dark and 1.2–1.3 mV with the iris light. The detector dominates, so the per-shot SNR at the centre is only about 5.5.

**Change to the plan:** open the iris to about 2 mm (about 4× the light, so about 4× the SNR; it adds only about 0.6 mm rms of blur to a 3 mm cloud), with 2 mm steps over ±12 mm, 3 shots each, shuffled: 0, −10, +4, +12, +8, 0, +2, −4, −8, −12, 0, +10, +6, −2, −6, 0.

### 2026-10-15 · MOT atom number: same formula as the Levine v2 note; aperture entry was 5× too large

**Reference:** J. Levine, "Calculate number of atoms from photodiode signal – v2 with collection optics" (PDA100A, f = 60 mm lens, 6.5 mm pinhole, 15 mm viewport at 180 mm, Ω = 5.4 × 10⁻³ sr).

**Formula:** the notebook's `mot_atom_number` is the note's Eq. 3, N = V·(4π/Ω)/(eG·hν·R), with R = (Γ/2)s/(1+s+(2Δ/Γ)²), s = 9/2.7 mW/cm² and Δ = −0.15Γ. `mot_collection_fraction` gives (1−cosθ)/2 = Ω/4π exactly.

**Fix:** `MOT_CAL["aperture_diameter_m"]` was 0.075; the note gives a 7.5 mm radius, so 0.015. The collection fraction falls from 0.0105 to 4.33 × 10⁻⁴ (Ω = 5.45 × 10⁻³ sr), and file 2026-10-15\_0002 (46.3 mV at coils off, 70 dB) goes from 5.8 × 10⁵ to **1.4 × 10⁷ atoms**.

**Still to confirm:**

- Load on ai2: the note's 2.35 × 10⁶ V/A is 70 dB into 50 Ω. `MOT_CAL` uses 4.75 × 10⁶ V/A (high impedance); with a 50 Ω terminator on ai2, N doubles to 2.8 × 10⁷.
- Optics: the note's Ω assumes the lens. A bare 10 × 10 mm diode against the viewport sees about 2.9 × 10⁻³ sr (equivalent aperture diameter 11.3 mm), which gives about 2.5 × 10⁷.

**Consistency:** 46 mV at high impedance is 23 mV into 50 Ω, close to the note's 21 mV, so the MOT is about the same size as in the note. With N\_upward ≈ 0.5 × 10⁷ from the absorption dip, about a third of the MOT reaches the probe on the way up (each number good to about 2×).

### 2026-10-15 · Probe loss explained: a 1 mm slit; the atoms are far from saturation

**Source:** the user. A 1 mm slit cuts the 1 mm × 18 mm sheet out of a round beam about 18 mm across.

**Light budget:** a centred 1 mm slit passes 9–11% of a Gaussian with 1/e² diameter 15–18 mm, so 7–10 µW of the 80–90 µW measured before it. The PDA sees 7.4 µW (6.1 V at 60 dB, 0.55 A/W). The slit accounts for essentially all of the "missing 90%", and the PDA level is the power in the sheet.

**Saturation:** about 8 µW in the sheet (horizontal 1/e² half-width 9 mm) gives a peak of about 0.07 mW/cm² per pass, 0.14 mW/cm² with the retroreflection: **s ≈ 0.05–0.15** (I\_sat 1.1–2.7 mW/cm²), not 0.5–2.

**Revises:**

- 2026-10-13/14 saturation estimate: "don't add probe power" was wrong. It assumed all 90 µW passed the slit.
- 0130: saturation explains at most about 7% of the 27% deeper dip at 5.5 V than at 9.1 V. Most of that trend is something else, most likely the laser frequency following the same thermal swing (the 2f/1f detuning indicator also tracks the PDA level). The detuning explanation of the scatter stands again, and a power servo alone would not fix it.
- Atom number (0130): with s ≈ 0.1 instead of 1 the cross-section is about 1.8× larger, so N\_upward ≈ 0.5 × 10⁷ (0.4–0.6 × 10⁷), about 1 × 10⁶ in the sheet on the return pass. Optical pumping into F=3 during the transit is only a few percent at this intensity.
- The 0130 recommendation to lower the probe to 30–40 µW is withdrawn.
- Expected return/upward ratio: the horizontal profile is the beam's own Gaussian, so the Gaussian end applies, about 0.44. The "missing half" of the returning atoms stands.

**Prediction, 3× the light in the sheet at 50 dB:** about 1.4–1.6× lock-in SNR. The noise at f\_m is mostly the detector's, so the signal grows faster than the noise; the absorption per atom falls only about 15% at s ≈ 0.3. At 50 dB, 3× keeps the PDA under 10 V across the thermal swing.

**Two ways to get it:** more power before the slit, or a cylindrical telescope of about 12:1 (for example +120 mm and −10 mm, 110 mm apart) that squeezes the beam to about 1.5 mm tall, so the slit passes about 80% instead of 9%. The sheet height drops out of the absorption signal; it only sets s.

**Slit checks:**

- Centre the beam vertically on the slit (maximize the PDA level). At 2 mm off-centre, vertical beam motion changes the transmitted power by about 10% per mm, so 0.5 µm/√Hz of motion would explain the probe's 40–100 Hz noise. Centred, the effect is second order.
- Profile the beam before the slit, in both directions: 10× more light, and the vertical size sets s.

### 2026-10-15 · Low-frequency noise: the probe beam carries most of it below 150 Hz

**Data:** probe-on no-launch shots 0059–0069 against blocked 0005, in 25 ms windows (40 Hz bins) along the record. **Script:** `noise_tests/low_freq_probe_noise.py`. Prompted by the user's reading of the blocked-vs-lit figure.

| Band (windows after 100 ms) | Probe on (µV/√Hz) | Blocked (µV/√Hz) | Probe beam's share of the noise power |
| --- | --- | --- | --- |
| 40–100 Hz | 295 | 71 | 94% |
| 120 Hz line | 103 | 20 | 96% |
| 140–400 Hz | 28 | 19 | 54% |
| 0.4–1 kHz | 8.8 | 8.1 | 16% |

- **40–100 Hz is steady probe noise:** about 300 µV/√Hz through the whole record (4.7 × 10⁻⁵/√Hz relative), only somewhat higher in the first 25 ms.
- **120 Hz and 140–400 Hz are mostly the shutter ringing down:** the 120 Hz component falls from 750 µV/√Hz at 12 ms to about 100 by 160 ms. At these frequencies the transient lasts well over 150 ms, not 40 ms, so the shutter should open at least 300 ms early.
- **Above 400 Hz** the detector electronics dominate.
- **For the lock-in** this noise is out of band: mixed onto the 243 µV carrier background it gives about 10⁻⁸ V/√Hz.
- **For DC absorption (modulation off, as in 0130) and slow normalization** it is the limit. A pick-off would help there if the noise arises before the pick-off point. Even the C12703-01 is good enough at 40–100 Hz, since its own floor is about 0.7–1.2 × 10⁻⁵/√Hz. Whether it arises in the laser and fiber or in the chamber and retro optics needs the coherence test.

### 2026-10-15 · Correction: the photodetector's own electronics dominate the noise at f\_m (blocked-detector run 0005)

**Shots:** 2026-10-15 **0005**, 20 files (10 upward, 10 return), PDA100A2 blocked at 60 dB. Settings were left from 0130 (`mod_amp` 0, `MRR_TRIG_BOOL` 1). Compared with the probe-on, no-launch shots 0059–0069. **Script:** `noise_tests/dark_vs_lit.py`; figure `dark_vs_lit_0005.png` beside it.

| Noise next to f\_m, enabled window | µV/√Hz | Share of noise power |
| --- | --- | --- |
| Probe on, total | 4.07 | 100% |
| Detector blocked (PDA100A2 electronics + NI ADC) | 3.51 | 74% |
| Shot noise at 6.1 V | 1.4–1.5 | 12–14% |
| Laser intensity noise (the remainder) | 1.4–1.5 | 12–13% |

- **Matches the spec:** the PDA100A2 manual gives 423 µV RMS over 9 kHz at 60 dB, about 3.6 µV/√Hz. Referred to the input that is 2.4 pA/√Hz, about twice the shot noise of the 4 µA photocurrent.
- **Laser RIN is about −130 dBc/Hz,** not −123. The earlier budget assigned the detector's noise to the laser by subtraction.
- **Low frequencies:** the large transient after the probe shutter opens is absent when blocked, so it is optical (beam motion or vibration). The steady 300 Hz–3 kHz noise and the 1.5–2 kHz bump are unchanged when blocked: detector electronics. Below 300 Hz about half the power needs light (150–300 Hz: 40 µV/√Hz lit, 24 blocked).
- **Lines:** the 38–121 kHz lines remain when blocked (electrical, about 3 µV, far from f\_m). The 3.17, 25.0 and 48 kHz lines need light. 120 Hz is mostly light-induced (119 µV/√Hz lit, 27 blocked).
- **Track:** blocked upward records (stage moving) and return records (stage still) are identical, so the moving stage adds no noise on this channel.
- **Not tested:** electrical crosstalk of the modulation drive into the PDA channel, since `mod_amp` was 0.

**What this changes:**

- The reference pick-off can only remove the laser part: at most about 13% of the noise power, under 7% in SNR. It stays useful for the slow drift only.
- Lowering the probe power would lower the SNR: the signal falls while the detector noise stays fixed.
- The levers now: a quieter detector (shot-noise-limited at the same light gives about 1.7× SNR), or more light on the detector without more at the atoms (3× the light at 50 dB gives about 1.8×). Only about 8 µW of the 80 µW probe reaches the PDA.

### 2026-10-15 · Noise sources from existing data: the "1/f" is mostly a probe-shutter transient, and f\_m is already optimal

**Data:** no-launch shots 0059–0069 (47 records) and launched 0049–0058. The probe is on for the whole 0.2 s record. **Scripts:** `noise_tests/noise_sources_existing.py` and `noise_tests/noise_time_resolved.py`; figure `noise_time_resolved_0059-0069.png` beside them.

- **Shutter transient:** for about 40 ms after the probe shutter opens (commanded 7 ms before the record starts), noise below 5 kHz is up to 50× higher (this holds above about 1 kHz; the 100–400 Hz part rings for more than 150 ms, see the low-frequency entry). A ring at 3.17 kHz falls from 311 µV/√Hz at 12 ms to about 7 by 50 ms, with a second peak at 1.3 kHz. The tail reaches into the enabled window: the 0.15–1 kHz band is 19 µV/√Hz at 45 ms versus 11 steady. Both passes show it, so it is the probe shutter, not the launch.
- **The lock-in band is not affected:** 6.5–9 kHz is steady at 3.3 µV/√Hz from 25 ms on.
- **Steady spectrum (after the sine ends):** rises below 500 Hz (about 60 µV/√Hz at 120 Hz, 13 at 500 Hz), a bump at 1.5–2 kHz (about 15), then about 4 µV/√Hz near f\_m. Shot noise is 1.7 µV/√Hz.
- **Not common electrical pickup:** the reference channel ai1 sits at 1.5–2 µV/√Hz with zero coherence with ai0 at every frequency. The narrow lines on ai0 (19.8, 25.0, 27.2 kHz and many at 40–120 kHz) are absent on ai1.
- **The noise at f\_m looks optical (wrong: see the blocked-detector entry above):** divided by the PDA response (single pole at 11 kHz), it is flat within ±1% from 6 to 14 kHz. That is a white input current noise 2.5× shot noise (RIN about 7 × 10⁻⁷/√Hz). Feedback-resistor Johnson noise would be about 0.2 µV/√Hz, and amplifier voltage noise would rise with frequency. Laser intensity noise is the likely source; a capped-PDA run would confirm it.
- **Light-level scaling inconclusive:** the DC level spanned only 5.6–6.7 V across these shots.

**Modulation frequency:** noise referred to the signal, relative to 6.1 kHz, is 1.09 at 4 kHz, 1.04 at 5 kHz, 0.99–1.01 at 8–14 kHz and 1.06 at 16 kHz. Changing f\_m is predicted to gain nothing. Keep it between 6 and 14 kHz and away from the 19.8, 25.0 and 27.2 kHz lines. This ignores the laser lock's FM response.

**Proposed runs** (exp12, return pass only, no launch, about 5 shots each):

1. PDA capped: the electronic floor, and whether the transient and lines are electrical.
2. Capped with the track controller off: which lines are the track's PWM.
3. ND 0.3 in front of the PDA, interleaved with normal: per-band scaling (optical noise halves, shot noise falls by √2, electronic noise stays).
4. Probe shutter opened about 200 ms earlier: should remove the transient (needs an exp12 change).
5. MOT PDA100A as a pick-off before the chamber, recorded on ai2 during the probe window: coherence with ai0 shows which noise is already in the beam before the chamber, and tests subtraction at f\_m.

The transient tail also sits in the first 20 ms of the enabled window that the DC-dip atom-number fit (0130) uses; opening the shutter earlier would clean that up too.

### 2026-10-15 · Reference pick-off: the C12703-01 can track the drift but not cancel the 6.1 kHz noise

**Question:** without a second PDA100A2, can a probe pick-off focused onto the Hamamatsu C12703-01 serve as the reference? The noise has several parts: Johnson noise, slow thermal drift, PWM-like pickup from the track, and large PDA baseline swings.

**C12703-01** ([datasheet](https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/c12703_series_kacc1215e.pdf)): Si APD module, φ3.0 mm, APD gain 30, −1.5 × 10⁸ V/W at 800 nm (negative output), DC to 100 kHz, NEP 20 fW/√Hz. Maximum input 50–60 nW. Gain stable to ±2.5% over 25 ± 10 °C.

|  | PDA100A2 (signal) | C12703-01 at 40 nW |
| --- | --- | --- |
| Light on detector | about 10 µW | 40 nW (limit 50–60 nW) |
| Photocurrent before gain | 5.3 µA | about 18 nA |
| Shot noise relative to DC (/√Hz) | 2.5 × 10⁻⁷ | 4–7 × 10⁻⁶ (APD excess-noise factor 1–3) |
| Versus laser noise at f\_m (about 6 × 10⁻⁷/√Hz) | below it | 7–12× above it |

- **It cannot cancel the 6.1 kHz laser noise.** Its own shot noise is 7–12× the laser noise it would have to remove. Even with the optimal weight it removes at most 2% of the laser-noise power. Sensitivity doesn't help: relative noise falls only with more photons, and the APD saturates at 60 nW.
- **Useful cancellation needs about 7 µW or more** on an ordinary photodiode with bandwidth well above 6 kHz and a few-volt output, so the NI's own floor (about 1.3 µV/√Hz) stays below the laser noise it carries.
- **It is a good slow monitor** of the probe power entering the chamber: 1–2 × 10⁻⁴ resolution in a 1 kHz bandwidth, far below the ±30% thermal swing. Its gain drifts about 0.25%/°C.
- **First test it allows:** is the 5-minute swing (PDA 5.8–10.2 V) already present before the chamber? If yes, the atoms see it: use it for the saturation correction or a power servo. If no, it comes from the detection path (retro mirror, focus onto the PDA, optics after the chamber) and the atoms see steady power.
- **Suggested setup:** a wedge pick-off just before the chamber (about 3 µW of 80 µW), ND 2.3 to bring it to about 15 nW (about −2 V), leaving room for the swing below the 50 nW limit. An 850 nm bandpass filter in a lens tube against room light. Record on ai3.

**Swap option for the fast reference:** the MOT-window PDA100A sees at most 24.5 nW above its dark level (17.9 nW of atom fluorescence at coils-off, 0130). A φ3 mm APD at the same spot would get about 1.3–1.7 nW, about 0.2 V: far from saturation. The C12703-01 could take over the MOT monitoring, freeing the PDA100A for the probe pick-off at about 10 µW. The MOT atom-number calibration would change (smaller aperture, new gain).

**Script:** scratch `mot_light_levels.py` (temporary).

### 2026-10-15 · Matched-filter noise calibration: the return-pass uncertainty is ±15%

- **Correction:** readings within one no-launch shot are correlated (correlation exp(−d²/4σ²) at spacing d). Counting them as independent overstated the precision of κ. The earlier notebook statement that readings one σ apart are close to independent was wrong.
- **Method:** bootstrap over whole no-launch shots (4000 resamples).
- **Return:** κ = 1.23, 68% range 1.00–1.36 (±15%, 10 shots). **Upward:** κ = 1.19, 68% range 1.14–1.23 (±4%, 37 shots).
- **Return matched-filter SNR:** 49, 68% range about 44–60 (previously quoted as 48.8 with no range).
- **Notebook:** `Matched_Filter_Shot_Averaging.ipynb` Step 4 now prints the bootstrap range for each pass.
- **To tighten:** about 30 more no-launch return shots would roughly halve the range (it scales as 1/√shots).
- **Script:** scratch `kappa_bootstrap.py` (temporary).

### 2026-10-15 · Correction: the upward-pass PMT is far outside its linear range

**Point raised by the user:** the H7421 photon-counting head is linear to about 1.5 Mcps (assumed 10% count loss there, implying a dead time of about 67 ns).

**Count rates:**

- **Upward:** the PMT reads 5.5–7.3 Mcps at the cloud peak (plus about 0.45 Mcps background). A measured 6.5 Mcps means about 11.5 Mcps true: **about 40% of counts lost at the peak**.
- **Return:** about 0.9 Mcps total, so about 6% loss; nearly linear.

**Model (dead time 67 ns, detector assumed to ignore photons during the dead time):**

- A true +10% atom-number change reads as **+5–7%** on the upward pass.
- Peak flattening widens a Gaussian fit: a true σ of 1.29 ms fits as 1.43–1.55 ms. That explains about half of the PMT (1.73 ms) vs absorption (1.29 ms) width difference; the rest is unexplained.

**Consequences:**

- **Withdrawn:** "the scatter left after the saturation correction (10%) matches the PMT's 10%". The upward PMT spread is compressed.
- **What stands:** absorption is linear in atom number (about 1.7% dip). Its 10% residual scatter is an **upper bound** on the atom-number variation; it also contains fit noise, imperfect saturation correction and cloud position changes.
- The upward matched-filter width should be about 1.3 ms (absorption width), not the PMT's 1.74 ms.
- The temperature (about 12 µK) is essentially unchanged, since it is set mostly by the return width, where the PMT is nearly linear.
- Correlations against the upward PMT amplitude keep their sign, but their sizes are compressed.

**Proposed calibration run:** about 15 upward shots at low probe power (10–20 µW, s ≪ 1, so absorption per atom barely depends on power), with an OD 1 filter in front of the PMT (peak about 0.7 Mcps, about 5% loss, correctable). That gives two linear, independent atom measures in the same shots. A short ND-filter step series would also measure the PMT dead time and how the detector behaves.

### 2026-10-15 · Absolute atom number from the direct absorption dip (0130)

**Shots:** exp12 sequence folder **0130**, 16 upward shots. `mod_amp = 0` (no modulation), `mod_offset` 9.5 mV (line centre), 80 µW, 60 dB. **Script:** `noise_tests\absorption_dip_atoms.py`.

**Method:** fit each shot's photodiode trace (10 kHz) inside the modulated window with a quadratic baseline times (1 − Gaussian dip). The dip is centred on the PMT cloud time. Then integrate the fractional dip and convert using the probe geometry, speed and cross-section.

**The dip:**

- Clearly resolved: peak 1.66% (free width) or 1.44% (width fixed to the PMT's); fit uncertainty about 0.02% per shot. Centre within 0.2 ms of the PMT.
- **The absorption profile is narrower than the PMT envelope:** σ = 1.29 ± 0.13 ms vs 1.73 ms. Unexplained; possibly the PMT's field of view weights the cloud differently. This also fits the upward matched-filter optimum near 1 ms.
- Integrated fractional dip: 53.9 ± 3.2 µs (free width).

**Probe saturation drives the dip depth:**

- On 13 clean shots (3 with noisy fits excluded), the dip anti-correlates with the PDA baseline (−0.72): deeper when less light reaches the PDA.
- At line centre, detuning jitter barely matters, so this is not detuning. The size (dip ratio 1.27 between 5.5 V and 9.1 V) matches saturation with s ≈ 1 at 9.1 V (withdrawn 2026-10-15: the 1 mm slit leaves s ≈ 0.05–0.15, which explains at most about 7%; see the slit entry). **The thermal swing changes the probe power at the atoms.**
- With that trend removed, the scatter falls from 14% to 10%, matching the PMT's 10% (plausibly the real atom-number variation).
- **Revises the previous entry:** the lock-in 1f–2f anticorrelation may come mostly from this saturation and power broadening rather than detuning jitter; the two can't yet be separated.

**Atom number, upward pass:**

- Inputs: speed 2.30 m/s; horizontal size 3.0 mm (assumed equal to the vertical size from the dip width; measured 3.2 ± 0.1 mm on 2026-10-16); sheet 1/e² half-width 9 mm; double pass; s ≈ 1.
- **N\_upward ≈ 0.9 × 10⁷ atoms** (0.7–1.1 × 10⁷ across polarization assumptions: σ₀ from 1.4 to 3.5 × 10⁻⁹ cm²).
- **Return pass:** about 0.2 × that in the sheet (about 2 × 10⁶). If the returning cloud is centred, the whole cloud is about 5 × 10⁶, roughly half the launched atoms. Temperature about 12 µK.
- **Accuracy: a factor of about 2.**
  - These push the true number higher: any detuning from line centre, and optical pumping into F=3 during the 0.4 ms in the sheet (repump shutter closed).
  - Also tens of percent: the sheet width, and the double-pass assumption (if the PDA sees only one pass, double N).

**Recommendation:** withdrawn 2026-10-15, since the atoms are at s ≈ 0.05–0.15 (see the slit entry). It read: lower the probe to about 30–40 µW (s ≈ 0.4–0.5), or add a probe power servo. That reduces the sensitivity to the thermal swing, but it would lower the lock-in SNR, because the noise at f\_m is detector-limited (corrected 2026-10-15).

### 2026-10-14 (night) · Cloud physics from the two passes (preliminary, hand calculation)

**Inputs:** the matched-filter means (0049–0058, 0127, 0129) and the PMT timing.

| Quantity | Value | How |
| --- | --- | --- |
| Speed through the probe | 2.30 m/s, both passes | crossings at 86.8 and 556.4 ms after drop; apex 321.6 ms, 27 cm above the probe |
| Vertical size at the probe | σ = 4.0 mm upward → 13.3 mm return | PMT transit σ 1.74 → 5.78 ms × speed |
| Cloud temperature | about 12 µK (velocity spread 27 mm/s) | growth of the size over the 470 ms between passes |
| Return / upward atoms in the sheet | about 0.19 | ratio of time-integrated 1f signals (peak × transit σ); same speed on both passes, so the per-atom signal is the same |
| Expected with no losses | about 0.43–0.51 | sheet ±9 mm wide vs horizontal σ 4 → 13 mm (isotropic assumption; upward σ\_x measured 3.2 mm on 2026-10-16), Gaussian vs flat-top sheet |

- **About half the atoms expected back through the sheet are missing.** Candidates: the returning cloud is off-centre from a tilted launch, background-gas collisions over about 0.5 s, or other losses. The PMT ratio is 0.27, the same range.
- **Assumptions:** the cloud expands equally in all directions; the probe detuning is the same on both passes.
- **Absolute atom number needs one calibration:** about 10 upward shots with `mod_amp = 0` and `mod_offset` at line centre (about 9.5 mV) to measure the DC absorption dip (about 1%). Combine it with beam area, cross-section (s ≈ 0.05–0.15, corrected 2026-10-15) and dwell time (0.43 ms in the 1 mm sheet). Then N\_return ≈ 0.19 × N\_upward. With `MOT_CAL` filled in, this also gives the launch efficiency.

### 2026-10-14 (night) · What causes the shot-to-shot scatter: detuning jitter, partly thermal

**User hypothesis:** the lab has no temperature control and the laser is temperature-sensitive.

**Data:** the same 25 return and 26 upward matched-filter shots (sessions 0049–0058, 0127, 0129).

**1. The probe power (PDA DC level as a thermal proxy) explains little of the scatter.** Removing a linear DC trend: return 22% → 18% pooled; upward 27% → 25%. Session 0049–0058 had a steady DC level (6.0–6.3 V) yet still 16% (return) and 34% (upward) scatter.

**2. The scatter is mostly laser-detuning jitter, not atom number.** 1f and 2f anti-correlate in every session:

| Session | Upward corr(1f, 2f) | Return corr(1f, 2f) | Upward corr(1f+2f, PMT) | Return corr(1f+2f, PMT) |
| --- | --- | --- | --- | --- |
| 0049–0058 | −0.72 | −0.21 | +0.94 | +0.76 |
| 0127 | −0.32 | −0.31 | +0.66 | +0.50 |
| 0129 | −0.69 | −0.22 | +0.86 | +0.92 |

Moving toward line center lowers 1f and raises 2f; atom-number changes would move them together. Combined (normalized 1f + 2f), the signal tracks the PMT closely, so atom number appears once the detuning effect is removed.

**3. Part of the detuning is thermal, part is fast jitter.**

- The detuning indicator (2f/1f, independent of atom number) correlates with the PDA DC swing in all six pass–session cases (+0.23 to +0.66), consistent with the user's temperature hypothesis acting on the laser frequency.
- Shot-to-next-shot correlation of 2f/1f is −0.09 to −0.45, so a large part of the jitter is fast rather than slow drift: lock noise, or incomplete settling after the launch jumps.

**Next:**

- Log the room or optics temperature per shot (thermistor on `ai3`; four channels at 250 kS/s exactly fill 1 MS/s) and the laser frequency or lock error signal.
- Analysis: a detuning-insensitive amplitude per shot from 1f and 2f together, which also gives a per-shot detuning estimate.
- Longer term: enclose the laser, fiber and probe optics.

### 2026-10-14 (night) · Cross-checks: shot-to-shot spread is not atom number; channel crosstalk found

**Script:** scratch cross-check of the matched-filter results (same 25 return and 26 upward shots).

- **The SNR stands.** The matched-filter SNR of about 49 (return) is a noise statement, measured directly on no-atom shots.
- **Correction: the 22–27% spread of the PDA amplitude is NOT atom number.** The PMT fluorescence amplitude does not follow it (correlation 0.27 return, −0.09 upward), and dividing by the PDA DC level does not reduce it.
- **It tracks the probe-power swing instead.** The amplitude falls when the PDA DC level rises: correlation −0.62 (0049–0058) and −0.68 (0127), but +0.08 in 0129, after the probe was lowered to 80 µW. Candidate cause: polarization drift changing both the light reaching the PDA and how the atoms absorb. **Fixing the swing should also reduce the shot-to-shot spread.**
- **NI crosstalk: `ai0` leaks into `ai1`.**
  - The DC level of the recorded reference (`RelJumpModulation`, `ai1`) is almost perfectly anti-correlated with the PDA level on `ai0`: −0.97 to −1.00 in 0126, 0127 and 0129, about −3 mV per V. That holds even before the sine starts, while the command is a fixed −150 mV.
  - This is crosstalk in the multiplexed ADC (channels sampled 1.3 µs apart; `ai1` likely driven from a high-impedance source).
  - **Harmless for the lock-in,** which uses only the sine's phase; the constant leak at f\_m is removed with the background.
  - **But the reference DC level cannot be used to monitor the laser offset.** It probably explains the old "3.2 mV measured vs 6–7 mV commanded" puzzle.
  - **Fix:** buffer `ai1` or check its source impedance; a grounded spare channel between them in the scan order also helps.
- The matched-filter notebook's statement that averaging is "limited by atom number" should read **limited by real shot-to-shot signal variation of unknown cause**, probe-power swing suspected.

### 2026-10-14 (night) · Matched filter: return-pass SNR ×4 from analysis alone

**What:** a new analysis, `Matched_Filter_Shot_Averaging.ipynb` with functions in `cloud_signal_estimation.py` (both fully commented).

- Each shot's lock-in signal is filtered with a Gaussian whose width equals the cloud transit width measured by the PMT, instead of the 0.5 ms triangular filter.
- The signal is read at the PMT-defined cloud time t\_c and projected onto the common signal phase.
- Noise is calibrated on no-launch shots.

**Data:** signal shots at `mod_offset` 6 mV, `mod_amp` 10 mV from 0049–0058, 0127 and 0129 (25 return, 26 upward; one clipped shot excluded). Noise shots: 0059–0068 and 0069.

**Transit widths (PMT):** upward σ = 1.74 ± 0.07 ms; return σ = 5.78 ± 0.25 ms. Arrival time varies by only 0.22 ms (upward) and 0.41 ms (return) from shot to shot.

| Per-shot SNR | Return | Upward |
| --- | --- | --- |
| Old metric (max \|Z\|, both quadratures) | 11.3 ± 1.9 | 135 ± 31 |
| Narrow filter at t\_c, projected | 12.0 ± 3.0 | 183 ± 53 |
| **Matched filter at t\_c, projected** | **48.8 ± 10.3** | **326 ± 82** |
| Gain (matched / narrow) | **4.1×** | 1.8× |

- **Noise calibration matters.** The white-noise prediction underestimates the wide-filter noise: κ = 1.23 (return) and 1.19 (upward), growing with width to 1.29–1.75 at 8 ms. There is extra slow noise near f\_m. All SNRs above are corrected for it.
- **The width scan confirms the matched-filter principle.** Return SNR is best at 4–6 ms (transit 5.8 ms). Upward is best at about 1 ms, slightly below the 1.74 ms transit because the slow noise grows faster at larger widths.
- **PMT alignment barely matters:** aligned and unaligned averages differ by less than 1%, because the arrival jitter is small compared with the widths.
- **Averaging 25 return shots:** noise-limited SNR 241, but the shot-to-shot scatter of the amplitude is 22% per shot, so the SNR including that scatter is 23.
  - The detection is no longer the limit: per-shot noise is about 2% of the signal.
  - The real atom-number fluctuation from shot to shot now dominates.
  - Part of the scatter here is from mixing three sessions (90 and 80 µW probe power).
- **Caveat:** the matched filter gives up time resolution, so use it for amplitude, not shape or timing. Recalibrate with new no-launch shots whenever probe power, gain or modulation change.

### 2026-10-14 (evening) · Modulation amplitude scan (0129): best at 10–14 mV

**Shots:** exp12 sequence folder **0129**, 80 files. `mod_amp` 4, 7, 10, 14, 20 mV × 8 reps × both passes, `mod_offset` 6 mV, same paired setup. **Probe power lowered from 90 to 80 µW:** no clipping, no MOT dropouts, nothing flagged. MOT 43–44 mV.

| `mod_amp` | Upward signed 1f | Upward 1f SNR | Upward 2f SNR | Return signed 1f | Return 1f SNR | Return 2f SNR |
| --- | --- | --- | --- | --- | --- | --- |
| 4 mV | 27.8 ± 2.6 mV | 89 ± 9 | 17 | 1.33 ± 0.39 mV | 6.3 ± 0.3 | 2.6 |
| 7 mV | 42.0 ± 3.9 | 106 ± 10 | 30 | 2.48 ± 0.16 | 9.7 ± 0.5 | 3.0 |
| 10 mV | 53.3 ± 4.2 | 127 ± 6 | 40 | 2.81 ± 0.24 | 10.6 ± 1.0 | 4.9 |
| 14 mV | 49.2 ± 9.1 | 121 ± 15 | 72 | 3.11 ± 0.24 | 11.2 ± 0.9 | 5.1 |
| 20 mV | 37.1 ± 4.6 | 94 ± 13 | 106 | 2.44 ± 0.32 | 8.8 ± 1.1 | 8.2 |

- **The best 1f amplitude is broad: 10–14 mV**, equal within errors on the return pass. 4 mV gives about half the SNR; 20 mV starts to fall off. **Keep `mod_amp = 0.010`** (12–14 mV might gain about 5%, within the noise).
- **2f grows with amplitude** and nearly matches 1f on the return pass at 20 mV (8.2 vs 8.8), even at an offset that favors 1f. Untested option: 2f at `mod_offset` ≈ 10 mV (line center) with `mod_amp` = 20 mV.
- **Conclusion:** the offset and amplitude are both near optimal and the maxima are broad, so the modulation settings are not a lever for large SNR gains. Remaining levers: laser intensity noise at f\_m, probe power stability (the 5-minute swing), and atom number and launch reliability.

### 2026-10-14 (evening) · Offset scan pass 2 repeated after relocking (0127): 6 mV is best

**Shots:** exp12 sequence folder **0127**, 80 files, same paired setup as 0126. The lasers were relocked first. MOT 44 mV, atoms in every shot.

**Clipping:** 10 of 80 files reached the PDA100A2's 10.2 V limit (runs 10, 26, 27, 29, 30, 57, 75, 76, 77, 79). Clipped upward shots had about 2.5× lower SNR than unclipped ones at the same offset. The scan notebook now flags and excludes any shot with a PDA sample at or above 9.95 V (`clipped`).

| `mod_offset` | Upward signed 1f | Upward SNR | Upward signed 2f | Return signed 1f | Return SNR (shots) |
| --- | --- | --- | --- | --- | --- |
| 2 mV | +51.8 ± 5.9 mV | 127 ± 14 | 27.9 mV | +2.53 ± 0.16 mV | 10.0 ± 0.7 (8) |
| **6 mV** | +48.5 ± 3.6 | 130 ± 10 | 32.5 | **+2.78 ± 0.19** | **11.7 ± 0.3 (7)** |
| 10 mV | −6.7 ± 10.8 | 115 ± 9 | 59.4 | +0.19 ± 0.83 | 6.9 ± 0.9 (6) |
| 14 mV | −30.0 ± 11.6 | 122 ± 16 | 36.2 | −2.24 ± 0.17 | 9.0 ± 0.6 (8) |
| 18 mV | −44.9 ± 3.3 | 119 ± 11 | 15.1 | −2.55 ± 0.19 | 9.6 ± 0.7 (5) |

- **Line center is at about 9.5–10 mV on both passes** (1f sign change; 2f largest at 10 mV).
- **Keep `mod_offset = 0.006`:** best return SNR, 11.7 ± 0.3, matching the afternoon baseline (11.6). The other offsets are 15–25% lower, and line center is worst for 1f.
- The upward SNR is nearly flat across offsets (115–130), so the return pass is what decides the setting.
- **The probe-power swing remains after relocking:** PDA 5.8–10.2 V, a smooth cycle of about 25 shots (about 5 min), likely thermal (fiber polarization or an etalon; 2026-10-16: the user traced it to a polarization mismatch at the PM fiber). Until it's fixed, lower the probe power about 20% or use 50 dB so the peaks stay under about 9 V.
- **Next:** `mod_amp` scan at 6 mV offset with the same paired setup.

### 2026-10-14 (evening) · Offset scan pass 2, clean run (0126): line center about 9 mV

**Shots:** exp12 sequence folder **0126**, 80 files. `mod_offset` 2, 6, 10, 14, 18 mV × 8 reps × both passes.

**Runmanager setup:** in the Axes tab, `t_PMT_delay` was the innermost axis and not shuffled, with `mod_offset` and `rep` shuffled. Every upward shot is followed by its return shot at the same offset and rep. Each block of 10 shots covers all five offsets in the same shuffled order, and the blocks run in shuffled rep order.

**No dropouts:** atoms in all upward shots, steady 11.2/12.4 s between shots, nothing flagged.

| `mod_offset` | Upward signed 1f | Upward SNR | Upward signed 2f | Return signed 1f | Return SNR |
| --- | --- | --- | --- | --- | --- |
| 2 mV | +24.3 ± 1.6 mV | 74 ± 2 | 12.7 mV | −0.93 ± 0.07 mV | 5.6 ± 0.2 |
| 6 mV | +20.5 ± 4.4 | 66 ± 10 | 29.1 | −0.84 ± 0.10 | 5.3 ± 0.5 |
| 10 mV | −9.3 ± 2.8 | 46 ± 6 | 34.7 | +0.59 ± 0.42 | 3.8 ± 0.4 |
| 14 mV | −20.5 ± 2.0 | 70 ± 8 | 12.6 | +0.75 ± 0.18 | 4.3 ± 0.3 |
| 18 mV | −21.2 ± 1.8 | 66 ± 8 | 8.3 | +1.26 ± 0.27 | 5.1 ± 0.5 |

- **Line center is now at about 8.8 mV** (1f zero crossing between 6 and 10 mV; 2f largest at 10 mV), about 2 mV lower than pass 1 earlier today. Both passes agree.
- **The 1f signal is flat on each side** (2 ≈ 6 mV, 14 ≈ 18 mV). The low side is best on the return pass. Keep 6 mV, or 2–4 mV for more margin from line center.
- **SNR roughly halved compared with the afternoon** (upward about 66 vs 146, return 5.3 vs 11.6 at 6 mV), because conditions changed, not because of the offsets:
  - the MOT signal is 33 mV instead of 42 mV;
  - **the probe power swings:** the PDA DC level cycles between 5.8 and 10.2 V over about 2 minutes, and 8 of 80 shots reach the PDA100A2's 10 V maximum, so those shots are clipped.
- **Next:** find and fix the probe power swing, and keep the PDA below about 9 V (lower the probe power slightly or use 50 dB). Then scan `mod_amp` at 6 mV with the same paired setup.

### 2026-10-14 · Offset scan pass 2 (0080): launch failed at run 20; failures coincide with a longer time between shots

**Shots:** exp12 sequence folder **0080**, 80 files. `mod_offset` 2, 6, 10, 14, 18 mV × 8 reps × both passes. **Shuffle was off**, so the shots ran in order: all upward by offset, then all return.

**Atoms vanish at run 20 (17:48:09, 10 mV) and never return** (60 shots). The MOT is normal in all 80 shots, and no shot is flagged.

**The clue: the time between shots grows by about 4 s at the moment of failure, in both failed scans.**

| Sequence | Last good shot | Time between shots before → after |
| --- | --- | --- |
| 0070 | run 05 | 9–10 s → 13–16 s |
| 0080 | run 19 | 11.2 s → 14–18 s |

- The sequence itself is unchanged, so the extra time is spent in BLACS between shots. Prime suspect: the BBD301 stage faulting or not moving, which would leave the MOT, shutters and probe normal but launch no atoms.
- The offsets at failure (+30 mV in 0070, 10 mV in 0080) are probably coincidental.
- The PDA DC also steps up from the same point (about 8 → 9.8 V); not yet explained.

**Usable data (upward only, runs 00–19):**

| Offset | Good shots | 1f SNR | Signed 1f |
| --- | --- | --- | --- |
| 2 mV | 8/8 | 142 ± 6 | +58 mV |
| 6 mV | 8/8 | 150 ± 10 | +56 mV |
| 10 mV | 4/8 | 34–116 | falling toward the zero crossing |

- 2 and 6 mV are equally good on the upward pass; 10 mV is close to line center (about 11 mV), as in pass 1.
- **No return-pass data.**
- **Next:**
  - check the BBD301 tab and BLACS log around 17:48:09, and whether the stage is at 190 mm and moves;
  - rerun pass 2 with **Shuffle on**, stopping if the time between shots jumps from about 11 s to about 15 s.

### 2026-10-14 · Repump unlocking: MOT lost on alternate shots (0071–0074)

- **0071 (before the relock):** both shots loaded a normal MOT (43 mV) but no atoms reached the probe.
- **0072–0074:** the first shot of each pair worked (upward 1f SNR 144 and 104; return SNR 12 in 0073). The second shot loaded **no MOT** (0.7–0.8 mV of atom signal), whichever pass it was.
- In the no-MOT shots the 61 → 75 mV scattered-light step is still there, so the cooling light and coils were fine and no atoms accumulated.
- **Cause (user):** the repump laser was unlocking. Fixed by relocking. An intermittent repump dropout may also explain the empty shots in 0070, rather than the +30 mV offset.
- **Scan notebook:** now flags shots with a MOT atom signal below 10 mV (`no MOT`) or `MRR_TRIG_BOOL = 0` (`no launch`), lists them, and by default excludes them from the summary and plots (`EXCLUDE_FLAGGED`). Tested on 0072–0074: all three no-MOT shots flagged.

### 2026-10-14 · Offset scan pass 1, second attempt: line center found, atoms lost after run 05

**Shots:** exp12 sequence folder **0070**, files `_00`–`_26`. Upward only, `MRR_TRIG_BOOL = 1`, `mod_offset` −18 to +30 mV × 3 reps, shuffled.

**Only runs 00–05 have atoms.**

- From run 06 on, all 21 shots show only noise (1f SNR about 2.3, PMT SNR about 2), at every offset. That includes 6 mV, which worked at run 04.
- Unchanged in the failed shots: the MOT signal (about 43 mV), the MRR shutter wait (no timeouts), the PDA DC level, the readback following each offset, and the commanded stage settings.
- So the loss happens between the MOT and the probe. The first failed shot (run 06) was also the first at +30 mV, and every later shot failed regardless of offset. The user's diagnosis: the lock was lost at +30 mV and did not recover. The unchanged MOT signal suggests only the jumped part of the lock failed (PGC, launch ramp and probe all use the relative jump; the MOT does not). Until the limit is known, keep \`mod\_offset\` within −12 to +18 mV, which worked before the failure.

**Good runs 00–05 (upward 1f, signed relative to the phase at 6 mV):**

| `mod_offset` | Signed 1f | 1f SNR | Signed 2f |
| --- | --- | --- | --- |
| −12 mV (2 shots) | +6.0 mV | 23 | −1.0 mV |
| 0 mV | +44 mV | 139 | −5 mV |
| 6 mV | +58 mV | 170 | +8 mV |
| 12 mV | −15 mV | 54 | +23 mV |
| 18 mV | −55 mV | 156 | −6 mV |

- **Line center is at about `mod_offset` = 11 mV.** The 1f signal changes sign between 6 and 12 mV (interpolated zero at 10.7 mV), and the 2f signal is largest at 12 mV.
- **The 1f maxima sit about ±6 mV from center.** The current 6 mV setting is already near the low-side maximum; 18 mV gives about the same magnitude with the opposite sign. 1f beats 2f by about 2.5×.
- **Next:**
  - check the stage and lock, and take one test shot at 6 mV;
  - then run pass 2 with both passes at about 2, 6, 10, 14 and 18 mV, 8 shots each.

### 2026-10-14 · Offset scan pass 1, first attempt: no launch (unusable as a scan)

- **Shots:** exp12 sequence folder **0069**, files `_00`–`_26`. Upward only (`t_PMT_delay` = 0.016), `mod_offset` = −18 to +30 mV in 6 mV steps × 3 reps, shuffled.
- **`MRR_TRIG_BOOL` was still 0** from the no-launch set, so no atoms reached the probe. The scan is being rerun with the launch on.
- **Usable as a no-atom reference at each offset:** the noise-only SNR is 2.2–3.0 and the PMT SNR is 2.1–2.7 at every offset. The offset does not change the probe background.
- **Scan notebook fixes:** a runmanager sequence writes every shot into one folder (`_00`, `_01`, …), so the notebook now:
  - takes the pass from each shot's `t_PMT_delay` instead of the filename;
  - ignores `t_PMT_delay` when auto-detecting the scan variable;
  - warns when any shot has `MRR_TRIG_BOOL = 0`.

  For one sequence folder, set `FIRST_SHOT = LAST_SHOT` to its number.

### 2026-10-14 (afternoon) · Repeat baseline with MOT monitor: 10 launched + 10 no-launch

**Shots:** exp12 2026-10-14 **0049–0058** launched and **0059–0068** no-launch (`MRR_TRIG_BOOL = 0`).

**Settings:** 250 kS/s (ai2 in use), 60 dB, 90 µW, `mod_amp` 0.01 V, `mod_freq` 6125 Hz, `mod_offset` 0.006 V.

**Scripts:** `noise_tests\baseline_exp12.py` (SNR now uses background >20 ms from cloud, plus MOT columns) and `noise_tests\no_atom_reference.py`.

| Returning `_1` | Launched 0049–0058 | No launch 0059–0068 | Morning launched 0005–0014 |
| --- | --- | --- | --- |
| 1f SNR, background >20 ms | **11.6 ± 1.7** | 2.4 ± 0.2 (noise only) | 10.3 ± 1.5 |
| 1f SNR, old ±8 ms window | 9.7 ± 1.4 | 2.4 | 8.8 ± 0.9 |
| 1f peak | 3.14 ± 0.43 mV | 0.62 mV (noise max) | 2.65 ± 0.38 mV |
| Lock-in noise, background >20 ms | 0.271 mV | 0.258 mV | 0.265 mV |
| Peak / PDA DC | 512 ppm | — | 400 ppm |
| PMT SNR | 6.4 | 2.6 | 5.2 |
| Noise at f\_m | 4.5 µV/√Hz | 4.1 µV/√Hz | 4.2 µV/√Hz |
| MOT atom signal | 42.4 ± 1.8 mV | 42.6 ± 0.8 mV | — |

- **About 13% better than the morning baseline, all from a larger signal**: more atoms reach the probe (PMT up too) and absorption per probe volt is higher. The recorded settings are the same.
- **250 kS/s costs nothing at f\_m.** The floor above 60 kHz rose from 0.87 to 1.3 µV/√Hz, as expected.
- **The far-background method is validated:** the launched noise measured >20 ms from the cloud (0.271 mV) is within 5% of the true no-atom noise (0.258 mV). With the old window it read 0.315 mV, about 20% too high.
- **No repeatable background, confirmed again:** 10 no-launch shots average down as √N (0.083 vs 0.082 mV), and subtracting them leaves the averaged launched SNR unchanged (12.9 → 12.9).
- **The scatter comes from after the MOT.**
  - The MOT atom signal varies 4% shot to shot; the return lock-in peak varies 14% and the upward one 28%.
  - The lock-in peak does not track the MOT (correlation 0.12 return, 0.40 upward), so dividing by the MOT signal doesn't reduce the spread.
  - On the upward pass, the PDA and PMT signals track each other (0.87), so atoms are lost or gained between MOT and probe: PGC, launch or transfer. The stage trigger jitter (0.1–6 ms) is a candidate.
- **For scans:** neither the MOT nor the PMT normalizes the return signal, so rely on shuffling and about 8 shots per setting.
- **MOT loading:** τ = 1.40 ± 0.13 s, scattered light 13 mV, so the 4 s load reaches about 94% of the asymptote.

### 2026-10-14 · First MOT loading curves; MOT shutter closes early

**Shots:** exp12 2026-10-14 0025–0026 (`_0` and `_1`), 250 kS/s, MOT PDA100A at 70 dB. **Script:** `noise_tests\mot_shutter_close.py`.

| MOT load (4 files) | Value |
| --- | --- |
| Dark level | 60.5–61.0 mV, the same before opening and after closing |
| Scattered light (B − dark) | 13.0–13.7 mV |
| Asymptote A | 44–45 mV |
| Loading time constant τ | 1.23–1.33 s (4 s load ≈ 3τ) |
| Atom signal at coils off | 42–44 mV |
| Fit residual | 1.1–1.3 mV |

**End of the trace, relative to coils off:**

- Coils off at 0 ms: no visible change.
- **PGC jump at +10 ms:** 116 → 81 mV. The atoms scatter less further red.
- **The shutter actually closes at about +17.7 ms:** 81 → 61 mV over about 1 ms, **2.3 ms before** the commanded close at +20 ms.

**Shutter timing:**

- The MOT shutter's opening matches the command to within about 0.1 ms, so the 6 ms open delay is right.
- The close delay is about **3.7 ms**, not the 6 ms in `connection_table_eedmW52.py` (`MOT_SHUTTER_do`, `delay=(6e-3, 6e-3)`, marked as a guess).
- **Effect:** PGC on the MOT arm runs about 7.7 ms instead of `t_redshift` = 10 ms.
- **Fix (not yet applied):** `delay=(6e-3, 3.7e-3)`, then recompile the connection table in BLACS.

**Analysis fixes:**

- "Dark after close" was NaN because it waited 10 ms after the commanded close, and the trace ends there; it now uses 3 ms after.
- The ipympl toolbar moved to the top (`toolbar_position = "top"`) in both notebooks, so it no longer covers the y-axis.

**Unexplained:** narrow downward spikes about every 0.73 s during loading (about 10 mV), seen in the 0026 `_0` curve.

### 2026-10-14 · MOT loading channel added (no data yet)

- **Hardware:** PDA100A pressed against a MOT window, cabled to NI6363 `ai2`.
- **Connection table:** added `AnalogIn ai2`. With three AI channels the 1 MS/s total allows at most 333 kS/s per channel, so `acquisition_rate` went 500e3 → **250e3** (divides the 100 MHz timebase exactly). This applies to every script using the table. BLACS must recompile the connection table because a channel was added.
- **exp12:** `ai2.acquire(label='MOT_fluorescence', start_time = t-0.01, end_time = t+t_load+0.04)`. The window covers 10 ms of dark before the MOT shutter opens, the full load, coils off, PGC and about 15 ms of dark after the shutter closes. That is about 1 M samples, roughly 13 MB per file.
- **Analysis:** `load_mot_trace` block-averages to 10 kHz and times the trace from the commanded shutter opening. `fit_mot_loading` fits V = B + A(1 − e^(−(t−t\_open)/τ)) from the measured opening to coils off. `mot_atom_number` converts volts to atoms.
- **Notebooks:** both have a new "MOT loading curve and atom number" section with a `MOT_CAL` dictionary. `None` entries skip the atom-number conversion.
- **Tested on a synthetic curve:** recovered dark, scatter, A, τ and the coils-off signal exactly. Old shots without ai2 skip the section cleanly.
- **Why it matters for the scans:** the MOT atom number is measured before the probe, so it normalizes atom-number swings independently of probe power and modulation depth. The PMT contrast cannot do that.
- **Needed for atom number:** PDA100A gain, MOT–aperture distance and clear aperture, window transmission, total MOT intensity, MOT detuning.

### 2026-10-14 · No-atom reference shots (no launch)

**Shots:** exp12 2026-10-14 0015–0019 with `MRR_TRIG_BOOL = False` (0015–0016 at 1:45 pm, 0017–0019 at 1:53 pm), compared with the launched baseline 0005–0014. Same settings otherwise. **Script:** `noise_tests\no_atom_reference.py`.

| Returning `_1` | Launched (10) | No launch (5) |
| --- | --- | --- |
| PMT SNR (noise-only max ≈ 2.5) | 5.2 | 2.5 |
| PMT smoothed peak / median | 1.9–2.8 | 1.19–1.36 |
| 1f "SNR" (largest noise excursion when no atoms) | 8.8 ± 0.9 | 2.4 ± 0.3 |
| Lock-in noise, background >10 ms from cloud | 0.299 mV | 0.233 mV |
| Lock-in noise, background >20 ms from cloud | 0.265 mV | 0.235 mV |
| 1f SNR, background >20 ms from cloud | **10.3 ± 1.5** | 2.4 ± 0.3 |

- **Clean reference:** no fluorescence bump on either pass (upward PMT SNR 2.4 vs 31 launched). Dropped atoms never reach the probe.
- **No repeatable instrument background.** Five no-launch shots average down exactly as √N (0.109 mV vs 0.104 expected on `_1`; 0.114 vs 0.111 on `_0`), and subtracting them from the launched average changes nothing.
- **The excess seen when averaging launched shots is the cloud's own wings.** The launched average is √N-limited more than 20 ms from the cloud (0.063 vs 0.059 mV), but 2.5× too noisy 10–20 ms from it (0.146 vs 0.059 mV).
- **The ±8 ms window lets cloud signal into the "background"**, so the per-shot noise was overestimated. With the background taken >20 ms from the cloud, the baseline SNR is 10.3 ± 1.5 rather than 8.8.
- A largest noise excursion of about 2.4 is what pure noise gives in a 16 ms window, so that's the false-detection level for a single shot.
- Probe level drifted: PDA DC 9.2 V for 0015–0016, 6.1 V for 0017–0019.
- **Conclusions:**
  - No-atom subtraction isn't needed at this level.
  - Widen the background guard (cloud ±8 ms plus 12 ms guard) for future SNR numbers.
  - Align shots on the PMT peak when averaging.

### 2026-10-14 · exp12 baseline, 10 shots

**Shots:** exp12 2026-10-14 0005–0014 (`_0` and `_1`). **Settings:** 500 kS/s, 60 dB, 90 µW, `mod_amp` 0.01 V, `mod_freq` 6125 Hz, `mod_offset` 0.006 V.

| Metric | Returning `_1` | Upward `_0` |
| --- | --- | --- |
| 1f SNR | 8.8 ± 0.9 | 151 ± 31 |
| 1f peak | 2.65 ± 0.38 mV | 54 ± 11 mV |
| Lock-in noise (complex RMS) | 0.30 ± 0.02 mV | 0.36 ± 0.05 mV |
| Peak / PDA DC | 400 ± 40 ppm | 7900 ± 900 ppm |
| 2f SNR | 3.2 ± 0.7 | 31 ± 16 |
| PMT SNR | 5.2 ± 0.6 | 31 ± 2 |

- Noise at f\_m 4.2 µV/√Hz (RIN −124 dBc/Hz); floor at 60–200 kHz 0.87 µV/√Hz.
- 1f signal phase 44 ± 5° on every shot: the laser sits consistently on one side of the line.
- **Analysis-method test:**
  - Projecting onto the fixed 44° phase did not raise the SNR (8.8 → 8.6), because the SNR already divides by the two-quadrature RMS. It still keeps the sign of the signal.
  - Averaging the 10 shots coherently gave SNR 15.2, not the 27 that √N predicts. The averaged peak drops to 1.7 mV because the cloud arrival time varies (PDA peaks 552.6–558.3 ms after drop), and the averaged noise stays at 0.11 mV. Corrected after the no-launch shots: that excess sits 10–20 ms from the cloud and comes from the cloud's own wings, not from a repeating background.
- **Conclusion:** use 5–8 shots per scan setting; align shots on the PMT peak before averaging (no-atom subtraction turned out unnecessary).

### 2026-10-14 · exp12 first shots: offset typo

- Shots exp12 0000–0001 had `mod_offset = 0.06` instead of 0.006. The readback mean during modulation was +40 to +48 mV, against −21 mV for exp11 0005 taken the same session, and the signal was poor. Those shots had 50 dB gain; the gain was set back to 60 dB before the baseline (confirmed by the user).
- Shots 0002–0004 were test shots after the fix (0004 had the PDA at 9.0 V).

### 2026-10-14 · exp12 created

- Copied `exp11_switchArmsAndLaunch.py` to `exp12_switchArmsAndLaunch_Absorption.py`; exp11 is unchanged.
- In exp12 only, the sine line reads `MAIN_JUMP_AMP_ao.sine(t, 0.1, mod_amp, 2*np.pi*mod_freq, 0, mod_offset, 1e5)`.
- With `mod_amp = 0.01`, `mod_freq = 6123`, `mod_offset = 0.006`, exp12 is identical to exp11.

### 2026-10-13/14 · Probe saturation estimate

- **Beam:** about 1 mm × 18 mm sheet, retroreflected, focused onto the PDA. Probe power measured: 90 µW.
- **On the detector:** 2.5 V at 50 dB (0.475 MV/A) ≈ 5.3 µA ≈ 10 µW at 0.55 A/W. About 90% of the probe does not reach the PDA (open question).
- **Saturation:** with 90 µW in the sheet, peak about 1.3 mW/cm² per pass, 2.5 mW/cm² with the retroreflection, so **s ≈ 0.5–2**.
- **Conclusions:**
  - Don't add probe power; it could reduce absorption per atom. (Wrong, 2026-10-15: a 1 mm slit passes only about 9% of the beam, so s ≈ 0.05–0.15.)
  - At most try 0.5×.
  - Choose the gain to put the PDA near 7–9 V at the chosen power, rather than as low as possible.
  - Radiation-pressure push is cancelled by the retroreflection; transit time through 1 mm is about 0.4 ms.

### 2026-10-13 · Sample rate 100 → 500 kS/s

**Change:** `acquisition_rate=500e3` on NI6363 in `connection_table_eedmW52.py`. Two channels (ai0, ai1) cap the per-channel rate at 500 kS/s.

**Problems hit on the way:**

- Editing `max_AO_sample_rate` on the NI6738 by mistake gave the BLACS error "connection table … not a subset". It's a connection-table property; reverted to 100 kHz.
- Setting the sine's sample rate argument to 500e3 gave a "clock rate … cannot support" error. That argument is the AO output rate, capped at 100 kHz by the clock line; reverted to 1e5.
- A runmanager `KeyError: 'labscriptlib'` after editing the connection table came from a stale compile subprocess. Fix: Restart subprocess.

**Shots:** exp11 2026-10-13 0000–0012 (100 kS/s) vs 0016–0029 (500 kS/s). Shots 0019–0023 have no data (never ran); 0024 has zero PMT counts (PMT, gating or counter fault, not lock loss).

| Metric | 100 kS/s (5 returning shots) | 500 kS/s (8 returning shots) |
| --- | --- | --- |
| Lock-in noise | 0.286–0.308 mV (mean 0.30) | 0.248–0.288 mV (mean 0.27) |
| Floor 20–45 kHz | 2.4–2.5 µV/√Hz | 1.6–1.8 µV/√Hz |
| PDA SNR ÷ PMT SNR | 1.38 | 1.33 |

- **Conclusion:** about 10% lower lock-in noise, as predicted (aliased ADC noise removed). The 500k shots looked worse because the clouds were weaker (lock trouble); normalized by the PMT they match.
- Analysis code was updated for 500 kS/s (see Files).

### 2026-10-12/13 · Initial noise budget (exp11)

**Shots:** exp11 2026-10-12 0013–0015 and 2026-10-13 0000–0012. Shots 2026-10-12 0010–0012 had no modulation (reference amplitude about 0) and are excluded.

- PDA at about 5.4 V (60 dB, about 7 µW on the detector). Noise near f\_m 4.7 µV/√Hz; white floor 2.4 µV/√Hz; 10–100 µV/√Hz below 1 kHz.
- The predicted lock-in noise 2 × 4.7 µV/√Hz × √1360 Hz = 0.35 mV matched the measured 0.355 mV, so the noise is stationary broadband noise at f\_m.
- The excess noise is the same with modulation on and off, so the FM doesn't create it.
- The PMT is shot-noise limited on about 380 kcounts/s of scattered light: predicted 64 kcounts/s lock-in noise, measured 62.
- Peak/background SNR 2.5–7 per shot; lock-in peak about 1–2 mV.
- **Conclusion:** chopping the probe would not help; balanced detection would (but only one detector is available); faster sampling helps a little.

## Files

All paths are under `C:\Users\EEDM\labscript-suite\userlib\labscriptlib\eedmLabscriptOps\` unless they start with `C:\Experiments`.

**Code changes**

| File | Change | Date |
| --- | --- | --- |
| `connection_table_eedmW52.py` | NI6363 `acquisition_rate` 100e3 → 500e3 (10-13) → 250e3 with \`ai2\` added (10-14) (NI6738 `max_AO_sample_rate` left at 100 kHz) | 2026-10-13 |
| `controlsAndMotDiagnosticAndLaunchTest_3_10_26\exp12_switchArmsAndLaunch_Absorption.py` | New; copy of exp11 with `mod_amp`, `mod_freq`, `mod_offset` globals in the sine line | 2026-10-14 |
| `Sandbox\Data_Labscript_Share\Scripts\Florescence_Absorption_10-6-26\fluorescence_absorption_lockin.py` | `load_shot`: AI rate read from file; separate `expected_pmt_fs_hz = 100 kHz`; accepts the exp12 sine line; returns `mod_amp`, `mod_freq`, `mod_offset`, `pda_gain_db`, `probe_power_uW` (None for exp11) | 2026-10-13/14 |
| `…\Returning_Fluorescence_Absorption_Lockin_10-6-26.ipynb` | `EXPECTED_FS_HZ = None`; spectrum plotted to Nyquist; `FM_HZ` taken from `mod_freq` when present; metadata table shows the test globals | 2026-10-13/14 |

Upward notebook: `Fluorescence_Absorption_Lockin_10-6-26.ipynb` (upward-pass notebook) set `EXPECTED_FS_HZ = 100_000.0` until 2026-10-14; it now uses None, takes FM\_HZ from mod\_freq, and has the MOT loading section. Both notebooks import the MOT functions (load\_mot\_trace, fit\_mot\_loading, mot\_collection\_fraction, mot\_atom\_number) from the helper.

**Scan notebook (2026-10-14):** `Absorption_Scan_Analysis.ipynb`, in the same folder. Set `DAY_DIR`, `FIRST_SHOT` and `LAST_SHOT`; the scan variable is the global that changes across the shots (`rep` is ignored). For each pass it reports, per scan value, the mean ± standard error of:

- signed 1f and 2f signals (projected on one reference phase, so they change sign across line center);
- 1f SNR (background >20 ms from the cloud), PMT SNR, MOT signal and PDA DC.

It writes `per_shot.csv`, `summary_by_scan_value.csv` and two figures to `DAY_DIR\scan_analysis_outputs\<first>-<last>\`. Tested on 0049–0068 (`MRR_TRIG_BOOL` 1 vs 0); the numbers match `baseline_exp12.py`.

**Analysis scripts** in `Sandbox\Data_Labscript_Share\Scripts\Florescence_Absorption_10-6-26\noise_tests\`. Each takes the helper's folder as its first argument and shot files after it, and runs with `C:\Users\EEDM\miniconda3\python.exe`.

| Script | What it does |
| --- | --- |
| `noise_budget.py` | PDA DC level and ASD in bands (0.3–1 kHz, near f\_m, at f\_m and 2f\_m, 20–45 kHz) for pre-modulation, enabled and post-enable windows; takes shot files only |
| `compare_rates.py` | Per-shot lock-in noise, peak, SNR, PMT SNR and band ASDs; used for the 100k vs 500k comparison |
| `baseline_exp12.py` | Per-shot table plus mean, std and spread per acquisition (`_0`, `_1`), including RIN, 1f phase and fractional absorption |
| `projection_average.py` | Fixed-phase projection and coherent multi-shot average of the 1f lock-in |
| `inspect_exp12.py` | Prints the globals, programmed AO range, PDA levels and reference mean per time window; takes shot files only |
| no\_atom\_reference.py | Launched vs no-launch shots: PMT check, per-shot noise, √N test of each group's average, and subtraction of the no-launch average. Arguments: helper dir, shots dir (one day), launched range, no-launch range, acquisition (0 or 1) |

**Data folders**

- exp11: `C:\Experiments\eedmLabscriptOps\exp11_switchArmsAndLaunch\2026\10\12` and `…\13`
- exp12: `C:\Experiments\eedmLabscriptOps\exp12_switchArmsAndLaunch_Absorption\2026\10\14`

## Noise budget reference

Correction (2026-10-15): the noise at f\_m is dominated by the PDA100A2's own electronics, not by laser intensity noise. With the detector blocked it reads 3.5 of the 4.1 µV/√Hz. The table below assigned that term to the laser by subtraction: read its excess row as mostly detector electronics, with a laser part of about 1.4 µV/√Hz.

Lock-in noise from a detector noise density S (V/√Hz), with the notebook's factor-2 normalization and the filter's noise bandwidth (about 1360 Hz for the 0.25 ms triangular filter):

```latex
\sigma_{|Z|} = 2\,S\,\sqrt{\mathrm{ENBW}}
```

Shot noise at the PDA output for DC level V and transimpedance gain G (60 dB ≈ 1.5 MV/A, 50 dB ≈ 0.475 MV/A, high-impedance load):

```latex
S_{\mathrm{shot}} = \sqrt{2\,e\,V\,G}
```

ADC noise aliased into the band (no anti-aliasing filter on the NI6363; σ ≈ 230–315 µV rms per sample on ±10 V) at sample rate f\_s:

```latex
S_{\mathrm{ADC}} = \sigma \,/\, \sqrt{f_s/2}
```

PMT background Poisson noise for count rate r, in (counts/s)/√Hz:

```latex
S_{\mathrm{PMT}} = \sqrt{2\,r}
```

| Term at f\_m (returning, 60 dB) | Oct 12–13 (5.4 V) | Oct 14 baseline (6.7 V) | Scales with probe power |
| --- | --- | --- | --- |
| Total measured | 4.7 µV/√Hz | 4.2 µV/√Hz | — |
| Excess (laser) | ≈ 3.7 µV/√Hz | ≈ 3.7 µV/√Hz | ∝ P |
| Shot noise | ≈ 1.6 µV/√Hz | ≈ 1.8 µV/√Hz | ∝ √P |
| Electronics + ADC | ≈ 1.5–2.4 µV/√Hz (100 kS/s) | 0.87 µV/√Hz floor (500 kS/s) | constant |

- **SNR vs probe power:**
  - laser-noise-limited: ∝ 1/(1+s);
  - shot-noise-limited: ∝ √s/(1+s), best at s = 1;
  - PMT against scatter background: also ∝ √s/(1+s).
- **Lock-in harmonics:** 1f follows the slope of the absorption line (zero at line center, sign flips across it); 2f follows the curvature (largest at center).
- Cs D2 numbers used: Γ/2π = 5.2 MHz, I\_sat = 1.1 (cycling) to 2.7 (isotropic) mW/cm², recoil 3.5 mm/s = 4.1 kHz Doppler per photon.

## Open questions and planned tests

**Planned tests, in order** (shuffle settings in runmanager; 5–8 shots per setting):

- [x] No-atom reference: about 10 shots with `MRR_TRIG_BOOL = False` at baseline settings; confirm no PMT bump
- [ ] Noise vs frequency, no atoms, modulation off: dark and probe-on at 60 dB, probe-on at about 0.5× power, 50 dB with matching power; 3 shots each
- [ ] `mod_offset` scan, e.g. 0, 0.003, 0.006, 0.009, 0.012 V; watch where the 1f X quadrature changes sign
- [ ] `mod_amp` scan at the best offset, e.g. 0.003, 0.006, 0.01, 0.02 V
- [ ] `mod_freq` test, only if the noise spectrum is lower elsewhere; above about 12 kHz needs a faster NI6738 output clock (connection-table change and BLACS recompile)
- [ ] Probe power 0.5× vs 1×, adjusting gain to keep the PDA at 7–9 V
- [ ] Locate the excess laser noise with no atoms: PDA after the fiber, before the chamber, at the return; tap the table and tilt the λ/4 at the return

**Analysis to build:** a batch script grouping shots by globals, with PMT-aligned averaging, no-atom background subtraction and error bars per setting.

**Open questions:**

- Where is the 90 µW probe power measured, and why does only about 10 µW reach the PDA?
- Laser lock drops out intermittently (suspected in 2026-10-13 shots 0016–0029); fix before the scans.
- Why did 2026-10-13 0024 record zero PMT counts?
- Is the probe path free of repump light (affects optical pumping into F=3)?
