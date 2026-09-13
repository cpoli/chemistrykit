r"""
Gaussian, Lorentzian, and Voigt lineshapes: two broadening mechanisms and their convolution
================================================================================================

Compares the three standard spectral lineshapes at matched FWHM: a
Lorentzian (the shape a classically radiating damped oscillator, or a
finite-lifetime quantum transition, produces -- "homogeneous"
broadening, identical for every molecule in the sample), a Gaussian (the
shape a Maxwell-Boltzmann spread of line-of-sight velocities produces via
the Doppler effect -- "inhomogeneous" broadening, different molecules
contributing different, momentarily Doppler-shifted, sub-lines), and
their convolution, the Voigt profile, which is what a real spectrum
typically shows once both mechanisms contribute together. The Lorentzian
has visibly heavier tails than the Gaussian at the same FWHM -- the
qualitative signature used to diagnose which mechanism dominates a real
measured line.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.spectro.utils.lineshapes import gaussian, lorentzian, voigt

fwhm = 4.0  # shared FWHM, arbitrary wavenumber-like units
x = np.linspace(-20.0, 20.0, 2001)

g = gaussian(x, x0=0.0, fwhm=fwhm)
lorentz = lorentzian(x, x0=0.0, fwhm=fwhm)
v = voigt(x, x0=0.0, fwhm_gaussian=fwhm, fwhm_lorentzian=fwhm)

print(f"Peak height -- Gaussian: {np.max(g):.4f}, Lorentzian: {np.max(lorentz):.4f}, Voigt: {np.max(v):.4f}")

# All three are normalized to unit area (over a wide enough grid to
# capture the Lorentzian's slowly-decaying tails):
x_wide = np.linspace(-2000.0, 2000.0, 400_001)
area_g = np.trapezoid(gaussian(x_wide, 0.0, fwhm), x_wide)
area_l = np.trapezoid(lorentzian(x_wide, 0.0, fwhm), x_wide)
area_v = np.trapezoid(voigt(x_wide, 0.0, fwhm, fwhm), x_wide)
print(f"Areas (Gaussian, Lorentzian, Voigt): {area_g:.4f}, {area_l:.4f}, {area_v:.4f}")

# %%
# The Lorentzian's heavier tails, at exactly matched FWHM, are the
# textbook way to tell homogeneous (lifetime/collisional) broadening
# apart from inhomogeneous (Doppler) broadening in a measured spectrum --
# far from line center the Lorentzian sits orders of magnitude above the
# Gaussian, even though both have fallen to the same half-maximum at
# +/-fwhm/2:

far_from_center = np.abs(x) > 3.0 * fwhm
ratio = np.mean(lorentz[far_from_center] / g[far_from_center])
print(f"\nAt |x| > 3*FWHM, the Lorentzian exceeds the Gaussian by a factor of order 1e{np.log10(ratio):.0f}")

# %%
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(x, g, label="Gaussian (Doppler / inhomogeneous)", color="steelblue")
ax.plot(x, lorentz, label="Lorentzian (lifetime / homogeneous)", color="crimson")
ax.plot(x, v, "--", label="Voigt (convolution of both)", color="black")
ax.set_yscale("log")
ax.set_ylim(1e-5, 1.0)
ax.set_xlabel("x - x0 (arbitrary units)")
ax.set_ylabel("intensity (log scale)")
ax.set_title(f"Lineshapes at matched FWHM = {fwhm:.0f}")
ax.legend()
fig.tight_layout()
plt.show()
