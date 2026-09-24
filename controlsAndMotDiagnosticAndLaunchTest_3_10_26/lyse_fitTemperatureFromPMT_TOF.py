# lyse_fitTemperatureFromPMT_TOF.py
# Plot PMT TOF fluorescence vs time since drop, fit a positive Gaussian.
# The PMT trace is CUMULATIVE photon counts, so differentiate to get a rate.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import lyse
import h5py

plt.rc('font', family='serif')
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', titlesize=16, labelsize=16)
plt.rc('legend', fontsize=30)
plt.rc('figure', titlesize=18)


def pos_gaussian(t, offset, amp, t0, sigma):
    return offset + amp * np.exp(-(t - t0) ** 2 / (2 * sigma ** 2))

def fmt(value, error, sig=2):
    """Format as value(error) with error in the last sig figures, e.g. 5.869(4)."""
    if error == 0 or not np.isfinite(error):
        return f"{value:.4g}"
    # number of decimal places so the error has `sig` significant figures
    exp = int(np.floor(np.log10(abs(error))))
    decimals = max(0, sig - 1 - exp)
    err_int = int(round(error / 10 ** (exp - (sig - 1))))
    return f"{value:.{decimals}f}({err_int})"


decimate = 10   # average this many adjacent bins together (1 = no decimation)

run = lyse.Run(lyse.path)

try:
    t_raw, cum_counts = run.get_trace('photon_counts')
except Exception:
    print("No photon_counts trace in this shot. Skipping.")
    raise SystemExit

# drop_time is an HDF5 attribute on the 'shot_properties' group -> read with h5py
with h5py.File(lyse.path, 'r') as f:
    try:
        t_PMT_delay = f['globals'].attrs['t_PMT_delay']
    except (KeyError, AttributeError):
        print("No t_PMT_delay global in this shot. Skipping.")
        raise SystemExit
# PMT delay is an estimate of the time between drop and when the PMT sees the fluorescence

# --- Cumulative counts -> count rate ---------------------------------------
# np.diff shortens the array by one, so use bin midpoints for the time axis.
dt = np.diff(t_raw)
rate = np.diff(cum_counts) / dt          # counts per second in each interval
t_mid = 0.5 * (t_raw[1:] + t_raw[:-1])   # midpoint of each interval
t = t_mid + t_PMT_delay                  # time since drop
values = rate


if decimate > 1:
    n = (len(values) // decimate) * decimate
    t = t[:n].reshape(-1, decimate).mean(axis=1)
    values = values[:n].reshape(-1, decimate).mean(axis=1)
 
# --- Initial guesses (positive peak) ---------------------------------------
offset0 = np.median(values)
amp0 = values.max() - offset0
t00 = t[np.argmax(values)]
sigma0 = (t.max() - t.min()) / 10.0
p0 = [offset0, amp0, t00, sigma0]

fit_ok = False
try:
    popt, pcov = curve_fit(pos_gaussian, t, values, p0=p0)
    perr = np.sqrt(np.diag(pcov))
    offset, amp, t0, sigma = popt
    sigma = abs(sigma)
    fit_ok = True

    residuals = values - pos_gaussian(t, *popt)
    resid_std = np.std(residuals, ddof=len(popt))   # ddof = number of fit params
    snr = amp / resid_std
    print(f"  residual std = {resid_std:.4g}")
    print(f"  SNR (B/std)  = {snr:.4g}")
    run.save_result('pmt_resid_std', resid_std)
    run.save_result('pmt_snr', snr)

    print("Positive Gaussian fit (PMT):")
    print(f"  offset = {fmt(offset, perr[0])}")
    print(f"  amp    = {fmt(amp, perr[1])}")
    print(f"  t0     = {fmt(t0, perr[2])}")
    print(f"  sigma  = {fmt(sigma, perr[3])}")
    print(f"  FWHM   = {2.3548 * sigma:.4g}")
    integrated_counts = amp * sigma * np.sqrt(2 * np.pi)
    print(f"  integrated counts (bkg subtracted) = {integrated_counts:.6g}")
    print(f"  estimate of atoms (1000photons/atom, 0.08% efficiency) = {integrated_counts/1000/0.008/0.1:.6g}")

    integrated_counts = amp * sigma * np.sqrt(2 * np.pi)   # area under Gaussian (bkg subtracted)
    print(f"  integrated counts (bkg subtracted) = {integrated_counts:.6g}")
    run.save_result('pmt_integrated_counts', integrated_counts)
    run.save_result('pmt_offset', offset)
    run.save_result('pmt_amp', amp)
    run.save_result('pmt_t0', t0)
    run.save_result('pmt_sigma', sigma)
    run.save_result('pmt_fwhm', 2.3548 * sigma)
except RuntimeError:
    print("Fit failed to converge.")

# --- Plot ------------------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(t, values, '.')

if fit_ok:
    t_fit = np.linspace(t.min(), t.max(), 1000)
    fit_label = (
        r'$f(t) = A + B\,e^{-(t-t_0)^2 / 2\sigma^2}$' + '\n'
        f"$A$ = {offset:.4g} cts/s\n"
        f"$B$ = {amp:.4g} cts/s\n"
        f"$t_0$ = {t0 * 1000:.4g} ms\n"
        f"$\\sigma$ = {sigma * 1000:.4g} ms\n"
        f"SNR = {snr:.3g}\n"
        f"Cts. - bkg = {integrated_counts:.4g} cts\n"
    )
    plt.plot(t_fit, pos_gaussian(t_fit, *popt), '-', label=fit_label)

plt.xlabel('Time since drop (s)')
plt.ylabel('PMT count rate (1/s)')
plt.legend(loc='upper right')
plt.grid(True)
plt.tight_layout()