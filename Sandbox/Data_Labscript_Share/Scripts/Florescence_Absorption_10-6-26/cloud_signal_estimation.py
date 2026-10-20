"""Matched-filter estimation of the cloud's absorption signal, and averaging of many shots.

Background: what the lock-in already does
-----------------------------------------
The photodiode voltage v(t) contains the absorption signal as a small oscillation at the laser
modulation frequency f_m (6.125 kHz). The lock-in in fluorescence_absorption_lockin.py

    1. "mixes" v(t) with exp(-i*theta(t)), where theta(t) is the modulation phase fitted to the
       recorded reference. This shifts the f_m component down to zero frequency (baseband), so
       the absorption becomes a slowly varying complex number Z(t) = X(t) + i*Y(t);
    2. low-pass filters the mixed signal to remove everything that is not near baseband.

Step 2 uses a short triangular filter (two boxcars of 2 modulation periods, about 0.5 ms long).
That keeps good time resolution, but the cloud takes much longer than 0.5 ms to cross the probe:
the PMT shows a Gaussian transit of sigma ~1.75 ms on the upward pass and ~5.8 ms on the return.
A filter much shorter than the signal lets through far more noise than necessary.

What this module adds
---------------------
1. A **matched filter**. For a pulse of known shape in white noise, the linear filter with the best
   signal-to-noise ratio is one whose impulse response has the same shape as the pulse. The cloud
   transit is close to Gaussian, so the matched filter is a Gaussian whose width equals the transit
   width measured independently by the PMT. Theory for a Gaussian pulse (width s) and a Gaussian
   filter (width w): peak signal ~ s / sqrt(s^2 + w^2) and noise ~ 1 / sqrt(w), so the SNR is
   largest at w = s. At w = s the peak is reduced to 71% of the raw amplitude, but the noise is
   reduced much more.

2. **Evaluation at a known time.** Instead of searching for the largest |Z| (which always finds a
   positive number, even in pure noise), each shot's signal is read at the cloud center time t_c
   found from the PMT. Combined with projection onto a fixed signal phase (X along the signal
   direction), the result is an unbiased amplitude with Gaussian noise: pure noise gives values
   centred on zero with standard deviation 1 in SNR units.

3. A **noise model**. Wide filters leave too little cloud-free data in the 70 ms window to measure
   noise directly. Instead, the noise is measured with the narrow triangular filter far from the
   cloud, and scaled to the wide filter assuming the noise is white near f_m (verified earlier: the
   lock-in noise matched 2*S*sqrt(ENBW)). For a normalized FIR kernel k, the output noise variance
   is proportional to sum(k^2), so sigma_wide = sigma_narrow * sqrt(sum(k_wide^2)/sum(k_narrow^2)).
   The notebook checks this on no-launch shots.

4. **Shot averaging aligned on the PMT.** Every shot is evaluated at its own t_c, which is the same
   as shifting each shot's trace so the cloud centres line up before averaging.

All functions are pure: they read the arrays in a `shot` dictionary from load_shot() and return
new arrays or dictionaries. Nothing is written to disk and nothing is plotted here.
"""

import numpy as np
from scipy import optimize, signal

from fluorescence_absorption_lockin import triangular_kernel


# ---------------------------------------------------------------------------------------------
# Demodulation and filtering
# ---------------------------------------------------------------------------------------------

def demodulate(shot, reference, harmonic=1):
    """Mix the photodiode signal down to baseband, *without* any low-pass filtering yet.

    This is exactly the mixing step of complex_lockin() in fluorescence_absorption_lockin.py, so
    the phase convention and the factor of 2 (which makes a sinusoid of peak amplitude A give
    |Z| = A after filtering) are identical. Returns (times, raw complex baseband signal).
    """
    t = shot["t_pda"]
    # Modulation phase at every sample, from the fit to the recorded reference channel.
    theta = 2 * np.pi * reference["f_hz"] * (t - shot["sine_start_s"]) + reference["phase_rad"]
    # Removing the median first keeps the large DC level (~6 V) from leaking into baseband.
    z_raw = 2 * (shot["pda_v"] - np.median(shot["pda_v"])) * np.exp(-1j * harmonic * theta)
    return t, z_raw


def narrow_kernel(fs_hz, f_mod_hz, n_periods=2):
    """The lock-in's existing short triangular filter (two boxcars of n_periods each)."""
    kernel, _ = triangular_kernel(fs_hz, f_mod_hz, n_periods)
    return kernel


def gaussian_kernel(fs_hz, sigma_s, truncate=3.0):
    """A Gaussian filter of width sigma_s, truncated at +-truncate*sigma and normalized to sum 1.

    Normalizing to sum 1 means a constant input passes unchanged, so filtered values stay in the
    same units (mV of lock-in amplitude). A Gaussian wider than ~0.1 ms also suppresses the 2*f_m
    term left over from mixing by many orders of magnitude, so no extra filter is needed.
    """
    n_half = int(np.ceil(truncate * sigma_s * fs_hz))
    x = np.arange(-n_half, n_half + 1) / fs_hz
    kernel = np.exp(-0.5 * (x / sigma_s) ** 2)
    return kernel / kernel.sum()


def apply_kernel(z_raw, kernel):
    """Filter a complex baseband signal with a symmetric kernel (centred, no time shift)."""
    return signal.fftconvolve(z_raw, kernel, mode="same")


def valid_samples(t, shot, kernel, fs_hz, switch_guard_s=0.001):
    """Samples whose whole filter window lies inside the modulation-enabled interval.

    Outside this interval the filter would mix in the switch-on/off transients or zero padding.
    """
    half = (len(kernel) - 1) / (2 * fs_hz)
    return ((t >= shot["enable_start_s"] + half + switch_guard_s)
            & (t <= shot["enable_stop_s"] - half - switch_guard_s))


def noise_scale(kernel_new, kernel_ref):
    """Ratio of output noise for two normalized kernels, assuming white input noise.

    For white noise of variance s^2 per sample, a FIR filter k gives output variance
    s^2 * sum(k^2). The ratio of standard deviations is therefore sqrt(sum(k_new^2)/sum(k_ref^2)).
    """
    return float(np.sqrt(np.sum(kernel_new ** 2) / np.sum(kernel_ref ** 2)))


# ---------------------------------------------------------------------------------------------
# Cloud timing from the PMT (independent of the photodiode)
# ---------------------------------------------------------------------------------------------

def pmt_cloud_fit(shot, f_mod_hz, n_periods=2):
    """Fit a Gaussian plus constant to the smoothed PMT count rate inside the enabled window.

    The PMT sees fluorescence from atoms in the probe sheet, so its envelope gives the time the
    cloud is centred in the probe (t_c) and how long it takes to cross (sigma). Smoothing over a
    whole number of modulation periods removes the f_m flicker of the fluorescence first.
    Returns a dict with center_s, sigma_s, amp_cps, offset_cps, ok (False if the fit failed).
    """
    kernel, _ = triangular_kernel(shot["fs_pmt"], f_mod_hz, n_periods)
    rate = signal.fftconvolve(shot["pmt_rate"], kernel, mode="same")
    t = shot["t_pmt"]
    inside = (t > shot["enable_start_s"] + 0.002) & (t < shot["enable_stop_s"] - 0.002)
    tt, rr = t[inside], rate[inside]
    i_max = int(np.argmax(rr))

    def model(x, amp, center, sigma, offset):
        return offset + amp * np.exp(-0.5 * ((x - center) / sigma) ** 2)

    p0 = [rr[i_max] - np.median(rr), tt[i_max], 0.003, np.median(rr)]
    try:
        p, _ = optimize.curve_fit(model, tt, rr, p0=p0, maxfev=5000)
        amp, center, sigma, offset = p
        ok = (amp > 0) and (tt[0] < center < tt[-1]) and (0.0002 < abs(sigma) < 0.03)
        return {"center_s": float(center), "sigma_s": float(abs(sigma)), "amp_cps": float(amp),
                "offset_cps": float(offset), "ok": bool(ok)}
    except RuntimeError:
        return {"center_s": np.nan, "sigma_s": np.nan, "amp_cps": np.nan,
                "offset_cps": np.nan, "ok": False}


# ---------------------------------------------------------------------------------------------
# Per-shot estimate
# ---------------------------------------------------------------------------------------------

def shot_traces(shot, reference, t_c, sigma_match_s, harmonic=1, n_periods=2,
                background_guard_s=0.020):
    """Filter one shot both ways and return everything needed for estimates and plots.

    t_c                 cloud centre (analysis-local seconds), normally from pmt_cloud_fit
    sigma_match_s       matched-filter width (the transit sigma for this pass)
    background_guard_s  background = narrow-filter samples more than this far from t_c

    Steps:
      1. demodulate the raw photodiode signal;
      2. narrow filter -> measure the constant background Z_bg (residual amplitude modulation
         and pickup at f_m without atoms) and the noise of each quadrature far from the cloud;
      3. subtract Z_bg from the raw baseband signal, so both filters see only the cloud;
      4. matched filter;
      5. predict the matched-filter noise from the narrow-filter noise (white-noise scaling).
    The complex traces are returned unprojected; projecting onto the signal phase happens once the
    common phase for the pass is known (see project()).
    """
    fs = shot["fs_pda"]
    f_mod = reference["f_hz"]
    t, z_raw = demodulate(shot, reference, harmonic)

    k_narrow = narrow_kernel(fs, f_mod, n_periods)
    z_narrow = apply_kernel(z_raw, k_narrow)
    valid_narrow = valid_samples(t, shot, k_narrow, fs)
    background = valid_narrow & (np.abs(t - t_c) > background_guard_s)
    if background.sum() < 100:
        raise ValueError("Too little cloud-free data for the background estimate")
    z_bg = np.mean(z_narrow[background])
    # Noise per quadrature (X and Y have equal variance for noise with no preferred phase);
    # averaging the two variances halves the scatter of the estimate.
    noise_narrow = float(np.sqrt(0.5 * (np.var(z_narrow[background].real)
                                        + np.var(z_narrow[background].imag))))

    k_match = gaussian_kernel(fs, sigma_match_s)
    z_match = apply_kernel(z_raw - z_bg, k_match)
    valid_match = valid_samples(t, shot, k_match, fs)
    return {"t": t, "z_narrow": z_narrow - z_bg, "valid_narrow": valid_narrow,
            "z_match": z_match, "valid_match": valid_match, "background": background,
            "z_bg": z_bg, "noise_narrow": noise_narrow,
            "noise_match": noise_narrow * noise_scale(k_match, k_narrow),
            "k_narrow": k_narrow, "k_match": k_match, "t_c": t_c}


def value_at(t, z, t_c):
    """Complex value of a filtered trace at time t_c (linear interpolation between samples)."""
    return complex(np.interp(t_c, t, z.real) + 1j * np.interp(t_c, t, z.imag))


def common_phase(complex_values):
    """Signal phase shared by a set of shots: the angle of their vector sum.

    Summing the complex values weights each shot by its own amplitude, so strong shots dominate
    and noisy, nearly-zero shots barely affect the result.
    """
    return float(np.angle(np.sum(complex_values)))


def project(z, phase_rad):
    """Component of z along the signal phase: the signed amplitude X_signal."""
    return np.real(np.asarray(z) * np.exp(-1j * phase_rad))
