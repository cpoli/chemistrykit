r"""
Calibration curves, and IUPAC LOD/LOQ
==========================================

:func:`~chemistrykit.analytical.systems.calibration.fit_calibration` fits
a straight-line calibration by ordinary least squares, and
:class:`~chemistrykit.analytical.systems.calibration.LinearCalibration`'s
:meth:`~chemistrykit.analytical.systems.calibration.LinearCalibration.lod`/:meth:`~chemistrykit.analytical.systems.calibration.LinearCalibration.loq`
implement the IUPAC :math:`3.3\sigma/m`/:math:`10\sigma/m` convention.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.analytical.systems.calibration import fit_calibration
from chemistrykit.analytical.visualizers.analytical_plots import plot_calibration_curve

rng = np.random.default_rng(42)
concentration = np.array([0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0])
true_slope, true_intercept, noise_sigma = 12.4, 0.8, 1.5
signal = true_slope * concentration + true_intercept + rng.normal(scale=noise_sigma, size=concentration.shape)

cal = fit_calibration(concentration, signal)
print(f"Fitted slope:     {cal.slope:.4f} (true: {true_slope})")
print(f"Fitted intercept: {cal.intercept:.4f} (true: {true_intercept})")
print(f"R^2:              {cal.r_squared:.6f}")
print(f"Residual std error s_y/x: {cal.residual_std_error:.4f}")
print(f"\nLOD = 3.3*s_y/x/m = {cal.lod():.4f}")
print(f"LOQ = 10*s_y/x/m  = {cal.loq():.4f}")
print(f"LOQ/LOD ratio (always exactly 10/3.3): {cal.loq() / cal.lod():.6f}")

# %%
# A noisier calibration set has a larger LOD:
noisy_signal = true_slope * concentration + true_intercept + rng.normal(scale=8.0, size=concentration.shape)
cal_noisy = fit_calibration(concentration, noisy_signal)
print(f"\nLow-noise LOD:  {cal.lod():.4f}")
print(f"High-noise LOD: {cal_noisy.lod():.4f}")

# %%
# Predicting an unknown concentration from a measured signal:
measured_signal = 55.0
predicted_conc = cal.predict_concentration(measured_signal)
print(f"\nMeasured signal {measured_signal} -> predicted concentration {predicted_conc:.4f}")

# %%
ax = plot_calibration_curve(concentration, signal, calibration=cal, lod=cal.lod())
plt.tight_layout()
plt.show()
