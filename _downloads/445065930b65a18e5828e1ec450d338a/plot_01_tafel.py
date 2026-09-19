r"""
Butler-Volmer kinetics and the Tafel-plot linearization
==========================================================

:func:`~chemistrykit.electrochem.systems.butler_volmer.butler_volmer_current_density`
gives the full nonlinear current-overpotential relationship for an
electrode reaction. At high overpotential its Tafel-linearized form
(:func:`~chemistrykit.electrochem.systems.butler_volmer.tafel_overpotential`)
becomes an excellent approximation; this example checks that
convergence numerically and then fits synthetic Tafel-regime data back
to recover the exchange current density and Tafel slope.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.electrochem.systems.butler_volmer import (
    butler_volmer_current_density,
    fit_tafel_plot,
    tafel_overpotential,
    tafel_slope,
)
from chemistrykit.electrochem.visualizers.electrochem_plots import plot_tafel

i0, alpha, n = 1.0e-6, 0.5, 1

# %%
# Convergence of the Tafel approximation to the full Butler-Volmer
# equation as overpotential grows.
etas = np.array([0.05, 0.10, 0.15, 0.20, 0.30])
i_full = butler_volmer_current_density(i0, etas, alpha=alpha, n=n)
eta_tafel = tafel_overpotential(i_full, i0, alpha=alpha, n=n, branch="anodic")
rel_error = np.abs(eta_tafel - etas) / etas
for eta, err in zip(etas, rel_error, strict=True):
    print(f"eta = {eta:.2f} V: Tafel relative error = {err:.2%}")

# %%
# Fit synthetic high-overpotential data to recover (i0, Tafel slope).
eta_fit_data = np.linspace(0.20, 0.40, 12)
i_fit_data = butler_volmer_current_density(i0, eta_fit_data, alpha=alpha, n=n)
fit = fit_tafel_plot(eta_fit_data, i_fit_data)
print(f"\ntrue i0   = {i0:.3e}, fitted i0   = {fit.exchange_current_density:.3e}")
print(f"true b    = {tafel_slope(alpha=alpha, n=n):.4f} V/decade, fitted b = {fit.tafel_slope:.4f} V/decade")

ax = plot_tafel(i0=i0, eta_range=(-0.4, 0.4), alpha=alpha, n=n, fit=fit)
plt.tight_layout()
plt.show()
