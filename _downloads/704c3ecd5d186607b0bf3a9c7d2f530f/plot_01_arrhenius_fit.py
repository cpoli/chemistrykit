r"""
Recovering an activation energy from an Arrhenius plot
=========================================================

Synthetic rate-constant-vs-temperature data is generated from a known
activation energy :math:`E_a` and pre-exponential factor :math:`A`, then
:func:`~chemistrykit.kinetics.systems.arrhenius.fit_arrhenius` recovers
both parameters back out via the classic Arrhenius-plot linearization
(:math:`\ln k` vs. :math:`1/T`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.systems.arrhenius import arrhenius_rate_constant, fit_arrhenius
from chemistrykit.kinetics.visualizers.kinetics_plots import plot_arrhenius

A_true, Ea_true = 4.2e12, 65_000.0  # J/mol
T = np.linspace(280.0, 360.0, 9)

rng = np.random.default_rng(0)
k_exact = arrhenius_rate_constant(A=A_true, Ea=Ea_true, T=T)
# A small amount of multiplicative noise, as a real kinetics measurement would have.
k_noisy = k_exact * (1.0 + rng.normal(0.0, 0.02, size=T.shape))

fit = fit_arrhenius(T, k_noisy)
print(f"true:   Ea={Ea_true:.0f} J/mol, A={A_true:.3e}")
print(f"fitted: Ea={fit.Ea:.0f} J/mol, A={fit.A:.3e}, R^2={fit.r_squared:.5f}")

# %%
# The Arrhenius plot: a straight line of slope :math:`-E_a/R`.

ax = plot_arrhenius(T, k_noisy, fit=fit)
ax.set_title(f"Arrhenius plot (fitted Ea = {fit.Ea / 1000:.1f} kJ/mol)")
plt.tight_layout()
plt.show()
