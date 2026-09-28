r"""
Scherrer equation: crystallite size from peak broadening
==========================================================

Scherrer (1918) showed that tiny crystallites give broadened powder
lines, with width :math:`\beta=K\lambda/(\tau\cos\theta)` for
crystallite size :math:`\tau`.
:func:`~chemistrykit.crystal.systems.xrd.scherrer_crystallite_size`
inverts this: here, copper powder patterns are broadened for three
crystallite sizes, and each size is recovered from the measured width of
the (111) line.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.crystal.systems.xrd import powder_xrd_peaks, scherrer_crystallite_size

wavelength_nm = 0.15418  # Cu-Kalpha
K = 0.9
peaks = powder_xrd_peaks("FCC", a=0.3615, wavelength=wavelength_nm, hkl_max=2)
two_theta = np.linspace(35.0, 80.0, 9001)


def broadened_pattern(size_nm):
    """Sum of Gaussian lines, each with the Scherrer width for crystallites of `size_nm`."""
    y = np.zeros_like(two_theta)
    for p in peaks:
        fwhm = np.degrees(K * wavelength_nm / (size_nm * np.cos(np.radians(p.two_theta / 2.0))))
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        y += p.multiplicity * p.relative_intensity * np.exp(-0.5 * ((two_theta - p.two_theta) / sigma) ** 2)
    return y


def measured_fwhm(y, center):
    """Full width at half maximum of the line nearest `center`, from half-maximum crossings."""
    window = np.abs(two_theta - center) < 3.0
    x, yy = two_theta[window], y[window]
    above = x[yy >= yy.max() / 2.0]
    return above[-1] - above[0]


fig, ax = plt.subplots(figsize=(8, 4))
first = peaks[0]
print(f"(111) line at 2theta = {first.two_theta:.2f} deg")
for size in (5.0, 20.0, 100.0):
    y = broadened_pattern(size)
    ax.plot(two_theta, y / y.max(), label=f"{size:.0f} nm crystallites")
    fwhm = measured_fwhm(y, first.two_theta)
    tau = scherrer_crystallite_size(fwhm, first.two_theta, wavelength_nm, shape_factor=K)
    print(f"true size {size:6.1f} nm: FWHM = {fwhm:.3f} deg -> Scherrer size = {tau:6.1f} nm")
    assert abs(tau - size) / size < 0.05
ax.set_xlabel(r"$2\theta$ (degrees)")
ax.set_ylabel("normalized intensity")
ax.set_title("Cu powder pattern: smaller crystallites, broader lines")
ax.legend()
plt.tight_layout()
plt.show()
