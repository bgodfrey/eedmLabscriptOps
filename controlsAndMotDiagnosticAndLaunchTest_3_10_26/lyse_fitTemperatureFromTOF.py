# singleshot_tof_fluorescence.py
# Plot TOF fluorescence vs time since drop, fit a negative Gaussian.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import lyse
import h5py

path_to_data = 'data/traces/TOF_florescence'
#path_to_data = 'data/traces/photon_counts'

def neg_gaussian(t, offset, amp, t0, sigma):
    return offset - amp * np.exp(-(t - t0) ** 2 / (2 * sigma ** 2))


with h5py.File(lyse.path, 'r') as f:
    print('DUMB IS DUMB: ' + str(lyse.path))
    try:
        data = np.array(f[path_to_data])
        modData = np.array(f['data/traces/RelJumpModulation'])
    except KeyError:
        print("No TOF_florescence trace in this shot. Skipping.")
        raise SystemExit


    # shot_properties are stored as attrs on the 'shot_properties' group
    try:
        drop_time = f['shot_properties'].attrs['drop_time']
    except (KeyError, AttributeError):
        print("No drop_time saved in this shot. Skipping.")
        raise SystemExit
    # drop_time is an HDF5 attribute on the 'shot_properties' group -> read with h5py
    with h5py.File(lyse.path, 'r') as f:
        try:
            t_PMT_delay = f['globals'].attrs['t_PMT_delay']
        except (KeyError, AttributeError):
            print("No t_PMT_delay global in this shot. Skipping.")
            raise SystemExit

t = data['t'] - drop_time  # time since drop
values = data['values']
modValues = modData['values']
print('values = ', values)

print("TOF dataset:", data)
print("TOF shape:", data.shape)
print("TOF dtype:", data.dtype)

print("Mod dataset:", modData)
print("Mod shape:", modData.shape)
print("Mod dtype:", modData.dtype)

# --- Initial guesses ---
offset0 = np.median(values)
amp0 = offset0 - values.min()
t00 = t[np.argmin(values)]
sigma0 = (t.max() - t.min()) / 10.0
p0 = [offset0, amp0, t00, sigma0]

try:
    popt, pcov = curve_fit(neg_gaussian, t, values, p0=p0)
    perr = np.sqrt(np.diag(pcov))
    offset, amp, t0, sigma = popt
    sigma = abs(sigma)
    fit_ok = True
    print("Negative Gaussian fit:")
    print(f"  offset = {offset:.4g} +/- {perr[0]:.2g}")
    print(f"  amp    = {amp:.4g} +/- {perr[1]:.2g}")
    print(f"  t0     = {t0:.4g} +/- {perr[2]:.2g}")
    print(f"  sigma  = {sigma:.4g} +/- {perr[3]:.2g}")
    print(f"  FWHM   = {2.3548 * sigma:.4g}")
except RuntimeError:
    fit_ok = False
    print("Fit failed to converge.")

plt.figure(figsize=(8, 5))
# plt.plot(t, values, '.', label='data')
plt.plot(t, values)
#plt.plot(t, modValues+9.662, zorder = -10000)


if fit_ok:
    t_fit = np.linspace(t.min(), t.max(), 1000)
    fit_label = (
        r'$f(t) = A - B\,e^{-(t-t_0)^2 / 2\sigma^2}$' + '\n'
        f"$A$ = {offset:.4g} V\n"
        f"$B$ = {amp*1000:.4g} mV\n"
        f"$t_0$ = {t0*1000:.4g} ms\n"
        f"$\\sigma$ = {sigma*1000:.4g} ms\n" +
        f"baseline $\\sigma$: {np.std(values[:300]*1e6):.4g} $\\mu$V\n"
    )
    plt.plot(t_fit, neg_gaussian(t_fit, *popt), '-', label=fit_label)

plt.xlabel('Time since drop (s)')
plt.ylabel('Photodiode Voltage')
plt.legend(loc='lower right', fontsize=21)
plt.tight_layout()