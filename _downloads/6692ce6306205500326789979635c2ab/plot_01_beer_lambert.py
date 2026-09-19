r"""
Beer-Lambert absorbance, and its deviation from linearity at high concentration
==================================================================================

Absorbance is exactly linear in concentration by the ideal Beer-Lambert
law. A small amount of instrumental stray light (light that reaches the
detector without having passed through the fully absorbing sample path)
causes the apparent absorbance to fall increasingly below the true value
as concentration -- and hence true absorbance -- increases, the classic
"rolling over" deviation seen in real UV-Vis instruments.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.spectro.systems.beer_lambert import absorbance, apparent_absorbance_with_stray_light
from chemistrykit.spectro.visualizers.spectro_plots import plot_beer_lambert_deviation

epsilon = 8000.0  # L mol^-1 cm^-1
path_length = 1.0  # cm
concentrations = np.linspace(1e-6, 5e-4, 200)

true_absorbance = absorbance(epsilon, concentrations, path_length)
apparent = np.array([apparent_absorbance_with_stray_light(epsilon, c, path_length, stray_light_fraction=0.002) for c in concentrations])

print(f"At the lowest concentration, true A = {true_absorbance[0]:.4f}, apparent A = {apparent[0]:.4f}")
print(f"At the highest concentration, true A = {true_absorbance[-1]:.4f}, apparent A = {apparent[-1]:.4f}")

# %%
# The deviation (true minus apparent) grows sharply once the true
# absorbance gets large enough that the fixed stray-light fraction is no
# longer negligible next to the (exponentially shrinking) transmitted
# intensity:

deviation = true_absorbance - apparent
print(f"Deviation at A~0.1: {deviation[np.argmin(np.abs(true_absorbance - 0.1))]:.5f}")
print(f"Deviation at A~3.5: {deviation[np.argmin(np.abs(true_absorbance - 3.5))]:.5f}")

# %%
ax = plot_beer_lambert_deviation(molar_absorptivity=epsilon, path_length=path_length, concentrations=concentrations, stray_light_fraction=0.002)
plt.show()
