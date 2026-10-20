"""Fluorescence/absorption lock-in helpers for the exp11 HDF5 acquisition setup.

These functions read one file, handle cumulative PMT counts, fit the recorded laser
reference, and demodulate either acquisition without any upward/returning labels.
All times are acquisition-relative seconds internally. Drop-relative times use the
stored commanded start and drop_time. Detector outputs retain their native units.
The loader verifies the 100 ms sine / 70 ms enable sequence used by these files.
No files are written, no figures are created, and no embedded shot code is executed.
"""

# Numerical and HDF5 tools shared by both notebooks; plotting stays in the notebooks.
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import signal, optimize

def decode(value):
    return value.decode() if isinstance(value, (bytes, np.bytes_)) else str(value)


def sample_rate(times, channel, expected_fs_hz=None):
    """Verify a uniform, increasing sampling grid and return its rate."""
    intervals = np.diff(times)
    dt = float(np.median(intervals))
    if dt <= 0 or not np.allclose(intervals, dt, rtol=1e-5, atol=1e-12):
        raise ValueError(f"{channel}: timestamps are not uniformly increasing")
    fs = 1 / dt
    if expected_fs_hz is not None and not np.isclose(fs, expected_fs_hz, rtol=1e-5):
        raise ValueError(f"{channel}: expected {expected_fs_hz:g} Hz, found {fs:g} Hz")
    return fs


def load_shot(path, expected_fs_hz=None, sine_duration_s=0.100, enable_duration_s=0.070,
              expected_pmt_fs_hz=100_000.0):
    """Read the three stored traces and align them to the commanded acquisition start.

    expected_fs_hz applies to the NI analog inputs; None accepts the rate stored in the file.
    The photon counter has its own fixed sample rate, checked against expected_pmt_fs_hz.
    """
    path = Path(path)
    with h5py.File(path, "r") as h5:
        # Copy datasets so the arrays remain available after closing the read-only file.
        pda = h5["data/traces/TOF_florescence"][()]
        reference = h5["data/traces/RelJumpModulation"][()]
        pmt = h5["data/traces/photon_counts"][()]
        ai_rows = h5["devices/NI6363/AI"][()]
        ai = {decode(row["label"]): row for row in ai_rows}
        drop_s = float(h5["shot_properties"].attrs["drop_time"])
        acquisition_start_s = float(ai["TOF_florescence"]["start"])
        counter_start_s = float(h5["devices/photon_counter"].attrs["t_start"])
        # Read the PMT delay explicitly so the plot origin can be checked against the counter trigger.
        pmt_delay_s = float(h5["globals"].attrs["t_PMT_delay"])
        metadata_fs = float(h5["devices/NI6363"].attrs["acquisition_rate"])
        markers = {decode(row["label"]): float(row["time"]) for row in h5["time_markers"][()]}
        embedded_script = decode(h5["script"][()])
        # exp12 takes these from the absorption_tests globals; exp11 shots lack them (None).
        test_globals = {name: (float(h5["globals"].attrs[name]) if name in h5["globals"].attrs else None)
                        for name in ("mod_amp", "mod_freq", "mod_offset", "pda_gain_db", "probe_power_uW")}

    # Preserve analog sample spacing but remove the different stored shot-time origins.
    t_pda = np.asarray(pda["t"], dtype=float) - float(pda["t"][0])
    t_ref = np.asarray(reference["t"], dtype=float) - float(reference["t"][0])
    t_counter = np.asarray(pmt["t"], dtype=float)
    fs_pda = sample_rate(t_pda, "TOF_florescence", expected_fs_hz)
    fs_pmt = sample_rate(t_counter, "photon_counts", expected_pmt_fs_hz)
    sample_rate(t_ref, "RelJumpModulation", expected_fs_hz)
    assert np.isclose(metadata_fs, fs_pda, rtol=1e-5)
    assert np.isclose(float(ai["RelJumpModulation"]["start"]), acquisition_start_s)
    assert np.allclose(t_pda, t_ref, atol=1e-12, rtol=0)
    assert decode(ai["TOF_florescence"]["units"]) == "Volts"
    if not np.isclose(counter_start_s - drop_s, pmt_delay_s, atol=1e-9, rtol=0):
        raise ValueError("The stored PMT delay disagrees with its acquisition start")

    # Difference cumulative counts in signed integers to avoid unsigned underflow.
    cumulative = np.asarray(pmt["values"], dtype=np.int64)
    increments = np.diff(cumulative)
    if np.any(increments < 0):
        raise ValueError("Counter reset or rollover: inspect before interpreting count increments")

    # An increment belongs to the interval between two counter samples, so use its midpoint.
    t_pmt = (t_counter[:-1] + t_counter[1:]) / 2 + counter_start_s - acquisition_start_s
    pmt_rate = increments / np.diff(t_counter)

    # This sequence has a 100 ms sine and a 70 ms enable interval; verify timing overrides.
    if not np.isclose(sine_duration_s, 0.100) or not np.isclose(enable_duration_s, 0.070):
        raise ValueError("Recheck the embedded script before changing this sequence timing")
    # Locate the sine in the saved shot schedule, rather than assuming 40–140 ms.
    sine_start_s = markers["Jump to blue"] - acquisition_start_s
    # exp11 hard-codes the amplitude; exp12 passes mod_amp. Both use a 0.1 s sine.
    if ("MAIN_JUMP_AMP_ao.sine(t, 0.1, 0.01" not in embedded_script
            and "MAIN_JUMP_AMP_ao.sine(t, 0.1, mod_amp" not in embedded_script):
        raise ValueError("The embedded sine settings differ; recheck the timing configuration")
    if "MAIN_REL_JUMP_do.go_low(t+0.07)" not in embedded_script:
        raise ValueError("The embedded enable interval differs; recheck the timing configuration")

    # Validate finite signals before passing them to FFTs or nonlinear fitting.
    for values in (pda["values"], reference["values"], pmt_rate):
        if not np.isfinite(values).all():
            raise ValueError(f"Non-finite samples in {path.name}")
    return {"path": path, "t_pda": t_pda, "pda_v": pda["values"].astype(float),
            "reference_v": reference["values"].astype(float), "t_pmt": t_pmt,
            "pmt_rate": pmt_rate, "fs_pda": fs_pda, "fs_pmt": fs_pmt,
            "cumulative_counts": cumulative, "drop_s": drop_s,
            "acquisition_start_s": acquisition_start_s,
            "acquisition_after_drop_s": acquisition_start_s - drop_s,
            "pmt_delay_s": pmt_delay_s,
            "stored_first_analog_s": float(pda["t"][0]),
            "sine_start_s": sine_start_s, "sine_stop_s": sine_start_s + sine_duration_s,
            "enable_start_s": sine_start_s,
            "enable_stop_s": sine_start_s + enable_duration_s,
            **test_globals}


def triangular_kernel(fs_hz, frequency_hz, n_periods):
    """Two boxcars combined into one symmetric, odd-length FIR averaging kernel."""
    # A boxcar spans an integer number of modulation periods, suppressing carrier harmonics.
    n = int(round(n_periods * fs_hz / frequency_hz))
    if n < 2:
        raise ValueError("Choose an averaging interval with at least two samples")
    box = np.ones(n) / n
    kernel = np.convolve(box, box)
    # The combined length is 2N-1, giving an exact central sample even when N is even.
    # Convolving the combined kernel once avoids the one-sample shift from two 'same'
    # convolutions of even-length boxcars.
    return kernel, n


def fit_modulation_reference(shot, frequency_guess_hz=2_000.0):
    """Fit frequency nonlinearly and offset/drift/sine/cosine coefficients linearly."""
    t = shot["t_pda"]
    # Exclude the start/stop steps of the 100 ms analog-output sine.
    mask = ((t >= shot["sine_start_s"] + 0.002)
            & (t <= shot["sine_stop_s"] - 0.002))
    u = t[mask] - shot["sine_start_s"]
    v = shot["reference_v"][mask]
    center = float(np.mean(u))

    def linear_fit(f_hz):
        # Including a linear drift keeps slow DC motion from biasing the sine fit.
        theta = 2 * np.pi * f_hz * u
        matrix = np.column_stack([np.ones_like(u), u - center,
                                  np.sin(theta), np.cos(theta)])
        coefficients = np.linalg.lstsq(matrix, v, rcond=None)[0]
        residual = v - matrix @ coefficients
        return coefficients, float(np.mean(residual**2))

    # A narrow search around the programmed frequency avoids distant periodic fit minima.
    best = optimize.minimize_scalar(lambda f: linear_fit(f)[1],
                                   bounds=(frequency_guess_hz - 5, frequency_guess_hz + 5), method="bounded",
                                   options={"xatol": 1e-7})
    if not best.success:
        raise RuntimeError("Modulation reference fit did not converge")
    coefficients, mse = linear_fit(best.x)
    offset, slope, sine, cosine = coefficients
    # a*sin(theta)+b*cos(theta) = hypot(a,b)*sin(theta+atan2(b,a)).
    phase = float(np.arctan2(cosine, sine))
    result = {"f_hz": float(best.x), "amp_v": float(np.hypot(sine, cosine)),
              "phase_rad": phase, "offset_v": float(offset), "slope_v_per_s": float(slope),
              "residual_rms_v": float(np.sqrt(mse)), "center_s": center}
    # Store the fitted full trace for inspection; the fit is used only within the sine interval.
    u_all = t - shot["sine_start_s"]
    result["fit_v"] = (offset + slope * (u_all - center)
                       + result["amp_v"] * np.sin(2 * np.pi * best.x * u_all + phase))
    return result


def complex_lockin(values, times, fs_hz, shot, reference, harmonic, n_periods, switch_guard_s=0.001):
    """Return the peak-amplitude-normalized complex harmonic envelope and edge-valid mask."""
    kernel, n = triangular_kernel(fs_hz, reference["f_hz"], n_periods)
    theta = (2 * np.pi * reference["f_hz"] * (times - shot["sine_start_s"])
             + reference["phase_rad"])
    # Center before mixing to reduce leakage of the large photodiode DC voltage.
    mixed = (values - np.median(values)) * np.exp(-1j * harmonic * theta)
    z = 2 * signal.fftconvolve(mixed, kernel, mode="same")
    half_support_s = (len(kernel) - 1) / (2 * fs_hz)
    # This mask excludes both zero-padding artifacts and physically disabled modulation.
    valid = ((times >= max(times[0] + half_support_s,
                           shot["enable_start_s"] + half_support_s + switch_guard_s))
             & (times <= min(times[-1] - half_support_s,
                            shot["enable_stop_s"] - half_support_s - switch_guard_s)))
    return z, valid, n, half_support_s


def analyze_shot(shot, reference, n_periods=2, cloud_center_s=None,
                 cloud_half_width_s=0.008, background_guard_s=0.002,
                 switch_guard_s=0.001, harmonics=(1, 2)):
    """Locate the PMT transit, then find harmonic peaks in that independent search window."""
    # Define a cloud center from PMT fluorescence, independently of the PDA absorption waveform.
    kernel, _ = triangular_kernel(shot["fs_pmt"], reference["f_hz"], n_periods)
    smooth = signal.fftconvolve(shot["pmt_rate"], kernel, mode="same")
    half_support = (len(kernel) - 1) / (2 * shot["fs_pmt"])
    pmt_valid = ((shot["t_pmt"] >= shot["enable_start_s"] + half_support + switch_guard_s)
                 & (shot["t_pmt"] <= shot["enable_stop_s"] - half_support - switch_guard_s))
    pmt_peak_index = np.flatnonzero(pmt_valid)[np.argmax(smooth[pmt_valid])]
    if cloud_center_s is None:
        cloud_center_s = float(shot["t_pmt"][pmt_peak_index])
    results, rows = {}, []
    for channel, time_key, value_key, fs_key, unit in [
        ("PMT", "t_pmt", "pmt_rate", "fs_pmt", "counts/s"),
        ("PDA100A2", "t_pda", "pda_v", "fs_pda", "V"),
    ]:
        times = shot[time_key]
        for harmonic in harmonics:
            z, valid, n, support = complex_lockin(shot[value_key], times, shot[fs_key],
                                                shot, reference, harmonic, n_periods, switch_guard_s=switch_guard_s)
            cloud = valid & (np.abs(times - cloud_center_s) <= cloud_half_width_s)
            background_mask = valid & (np.abs(times - cloud_center_s)
                                      >= cloud_half_width_s + background_guard_s)
            if not np.any(cloud) or background_mask.sum() < len(kernel):
                raise ValueError("Insufficient settled cloud/background data")
            # Subtract the coherent off-cloud pickup before taking a positive magnitude.
            background = np.mean(z[background_mask])
            clean = z - background
            amplitude = np.abs(clean)
            noise_rms = float(np.sqrt(np.mean(np.abs(clean[background_mask])**2)))
            peak_index = np.flatnonzero(cloud)[np.argmax(amplitude[cloud])]
            full_peak_index = np.flatnonzero(valid)[np.argmax(amplitude[valid])]
            # Keep the largest off-cloud fluctuation as a useful empirical comparison.
            background_max = float(np.max(amplitude[background_mask]))
            peak_time_s = float(times[peak_index])
            results[(channel, harmonic)] = {"t": times, "z": clean, "r": amplitude,
                "valid": valid, "cloud": cloud, "background_mask": background_mask,
                "background": background, "noise_rms": noise_rms, "peak_index": peak_index,
                "peak_time_s": peak_time_s, "peak_amp": float(amplitude[peak_index]),
                "boxcar_ms": 1e3 * n / shot[fs_key], "support_ms": 2e3 * support}
            rows.append({"Channel": channel, "Harmonic": harmonic, "Units": unit,
                "Peak after acquisition (ms)": 1e3 * peak_time_s,
                "Peak after drop (ms)": 1e3 * (peak_time_s + shot["acquisition_after_drop_s"]),
                "Peak amplitude": float(amplitude[peak_index]),
                "Background complex RMS": noise_rms, "Background max amplitude": background_max,
                "Peak/background RMS": float(amplitude[peak_index]) / noise_rms,
                "Largest peak in full settled gate after acquisition (ms)": 1e3 * float(times[full_peak_index]),
                "Largest peak in full settled gate after drop (ms)":
                    1e3 * (float(times[full_peak_index]) + shot["acquisition_after_drop_s"]),
                "Phase at peak (deg)": np.rad2deg(np.angle(clean[peak_index])),
                "Boxcar length (ms)": 1e3 * n / shot[fs_key]})
    return {"results": results, "table": pd.DataFrame(rows),
            "cloud_center_s": cloud_center_s, "pmt_smooth": smooth,
            "pmt_peak_s": float(shot["t_pmt"][pmt_peak_index])}


# ---------------------------------------------------------------------------------------
# MOT loading curve (PDA100A against a window, NI6363 ai2, exp12 shots from 2026-10-14 on)
# ---------------------------------------------------------------------------------------

MOT_TRACE = "MOT_fluorescence"
# Cs D2 constants for the fluorescence-to-atom-number conversion.
CS_D2_WAVELENGTH_M = 852.347e-9
CS_D2_GAMMA_HZ = 5.234e6          # natural linewidth Γ/2π
PLANCK_J_S = 6.62607015e-34
LIGHT_SPEED_M_S = 299792458.0


def load_mot_trace(path, decimated_fs_hz=10_000.0):
    """Read the MOT photodiode trace, block-average it to decimated_fs_hz, and time it from
    the commanded MOT shutter opening ("Loading MOT" marker). Returns None if the file has
    no MOT trace (shots taken before ai2 was wired)."""
    with h5py.File(Path(path), "r") as h5:
        if f"data/traces/{MOT_TRACE}" not in h5:
            return None
        trace = h5[f"data/traces/{MOT_TRACE}"][()]
        markers = {decode(row["label"]): float(row["time"]) for row in h5["time_markers"][()]}
    t = np.asarray(trace["t"], dtype=float)
    v = np.asarray(trace["values"], dtype=float)
    fs = sample_rate(t, MOT_TRACE)
    # Averaging whole blocks lowers the noise as well as the rate; a partial last block is dropped.
    n = max(1, int(round(fs / decimated_fs_hz)))
    m = len(v) // n * n
    t0 = markers["Loading MOT"]
    return {"t": t[:m].reshape(-1, n).mean(axis=1) - t0,
            "v": v[:m].reshape(-1, n).mean(axis=1),
            "fs_raw_hz": fs, "fs_hz": fs / n, "block": n,
            "coils_off_s": markers["MOT coils off"] - t0,
            "shutter_close_s": markers["Close MOT shutter -- dark time"] - t0}


def fit_mot_loading(mot, settle_s=0.003, steady_window_s=0.05):
    """Fit V(t) = B + A·(1 − exp(−(t − t_open)/τ)) between the measured shutter opening and
    coils off. B − dark is scattered light with (almost) no atoms; A·(1 − exp(−T/τ)) is the
    atom fluorescence when the coils switch off."""
    t, v = mot["t"], mot["v"]
    pre = t < -0.003                       # before the commanded opening, allowing for shutter jitter
    if pre.sum() < 5:
        raise ValueError("Too little dark data before the MOT shutter opens")
    dark_v, dark_std_v = float(np.median(v[pre])), float(np.std(v[pre]))
    # The MOT shutter measured closing ~2.3 ms before its commanded time (1 ms edge), so 3 ms
    # after the command is safely dark; the trace ends 10 ms after the command.
    closed = t > mot["shutter_close_s"] + 0.003
    dark_after_v = float(np.median(v[closed])) if closed.sum() >= 5 else np.nan
    load = (t >= -0.003) & (t < mot["coils_off_s"])
    end = (t >= mot["coils_off_s"] - steady_window_s) & (t < mot["coils_off_s"])
    end_v = float(np.mean(v[end]))
    # Measured opening: first sample clearly above the dark level.
    threshold = dark_v + max(6 * dark_std_v, 0.02 * (end_v - dark_v))
    above = np.flatnonzero(load & (v > threshold))
    if len(above) == 0:
        raise ValueError("No light step found after the MOT shutter opens")
    t_open = float(t[above[0]])
    fit = (t >= t_open + settle_s) & (t < mot["coils_off_s"])

    def model(tt, b, a, tau):
        return b + a * (1 - np.exp(-(tt - t_open) / tau))

    b0 = float(v[fit][0])
    popt, pcov = optimize.curve_fit(model, t[fit], v[fit], p0=[b0, max(end_v - b0, 1e-6), 0.5],
                                    bounds=([-np.inf, 0.0, 1e-3], [np.inf, np.inf, 100.0]))
    b, a, tau = popt
    err = np.sqrt(np.diag(pcov))
    load_time = mot["coils_off_s"] - t_open
    return {"dark_v": dark_v, "dark_std_v": dark_std_v, "dark_after_close_v": dark_after_v,
            "t_open_s": t_open, "scatter_v": b - dark_v, "baseline_v": b, "baseline_err_v": err[0],
            "asymptote_v": a, "asymptote_err_v": err[1], "tau_s": tau, "tau_err_s": err[2],
            "load_rate_v_per_s": a / tau,
            "atoms_v_at_coils_off": a * (1 - np.exp(-load_time / tau)),
            "end_minus_baseline_v": end_v - b,
            "residual_rms_v": float(np.std(v[fit] - model(t[fit], *popt))),
            "fit_mask": fit, "model": lambda tt: model(tt, *popt)}


def mot_collection_fraction(aperture_diameter_m, distance_m):
    """Fraction of isotropic fluorescence entering a circular aperture on axis: (1 − cos θ)/2."""
    theta = np.arctan(aperture_diameter_m / 2 / distance_m)
    return (1 - np.cos(theta)) / 2


def mot_atom_number(atoms_v, gain_v_per_a, responsivity_a_per_w, collection_fraction,
                    total_intensity_mw_cm2, detuning_hz, isat_mw_cm2=2.7, transmission=1.0):
    """Atom number from MOT fluorescence voltage. Scattering rate per atom:
    R = (Γ/2)·s/(1 + s + (2δ/Γ)²), s = I_total/I_sat. isat 2.7 mW/cm² is the isotropic-
    polarization value usually used for a six-beam MOT (1.1 for the cycling transition)."""
    power_w = atoms_v / gain_v_per_a / responsivity_a_per_w / transmission
    photon_j = PLANCK_J_S * LIGHT_SPEED_M_S / CS_D2_WAVELENGTH_M
    gamma = 2 * np.pi * CS_D2_GAMMA_HZ
    s = total_intensity_mw_cm2 / isat_mw_cm2
    rate = gamma / 2 * s / (1 + s + (2 * detuning_hz / CS_D2_GAMMA_HZ) ** 2)
    return power_w / (photon_j * rate * collection_fraction)
