r"""
Onsager's exact solution of the two-dimensional Ising model
==============================================================

Onsager (1944) computed the free energy of the square-lattice Ising model
exactly: the first exact solution of a model with a genuine phase
transition. :class:`~chemistrykit.statmech.Ising2DOnsager` evaluates his
free energy, internal energy, and heat capacity, together with Yang's
(1952) spontaneous magnetization :math:`m=[1-\sinh^{-4}(2J/k_BT)]^{1/8}`.
The heat capacity diverges logarithmically at :math:`T_c`, and the
magnetization vanishes with exponent :math:`\beta=1/8`. Neither the
mean-field (Bragg-Williams) theory, shown for comparison, nor Ising's own
one-dimensional chain reproduces these features.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from chemistrykit.constants import K_B
from chemistrykit.statmech import Ising2DOnsager

J = 100.0 * K_B
model = Ising2DOnsager(coupling=J)
Tc = model.critical_temperature
print(f"T_c = {Tc:.2f} K")

T = np.linspace(20.0, 2.0 * Tc, 500)
T = T[np.abs(T - Tc) > 1e-6 * Tc]


# Mean-field comparison: m = tanh(4 J m / k_B T), T_c(MF) = 4 J / k_B
def mean_field_m(t):
    if t >= 4 * J / K_B:
        return 0.0
    return brentq(lambda m: m - np.tanh(4 * J * m / (K_B * t)), 1e-9, 1.0)


fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
axes[0].plot(T / Tc, model.spontaneous_magnetization(T), color="steelblue", label="Onsager-Yang (exact)")
axes[0].plot(T / Tc, [mean_field_m(t) for t in T], color="gray", linestyle="--", label="mean field")
axes[0].set_xlabel(r"$T / T_c$")
axes[0].set_ylabel("spontaneous magnetization")
axes[0].set_title(r"Order parameter, $\beta = 1/8$")
axes[0].legend()

axes[1].plot(T / Tc, model.internal_energy_per_spin(T) / J, color="steelblue")
axes[1].plot(T / Tc, model.free_energy_per_spin(T) / J, color="crimson")
axes[1].annotate("u / J", (1.5, -0.8), color="steelblue")
axes[1].annotate("f / J", (1.5, -2.9), color="crimson")
axes[1].axvline(1.0, color="gray", linestyle=":", linewidth=0.8)
axes[1].set_xlabel(r"$T / T_c$")
axes[1].set_ylabel("energy per spin / J")
axes[1].set_title("Energy and free energy")

axes[2].plot(T / Tc, model.heat_capacity_per_spin(T) / K_B, color="steelblue")
axes[2].axvline(1.0, color="gray", linestyle=":", linewidth=0.8)
axes[2].set_ylim(0.0, 4.0)
axes[2].set_xlabel(r"$T / T_c$")
axes[2].set_ylabel(r"$c / k_B$ per spin")
axes[2].set_title("Logarithmically divergent heat capacity")
fig.tight_layout()

# %%
# Close to T_c the heat capacity grows like -(8/pi) K_c^2 ln|1 - T/T_c|:
# each factor-of-ten step closer to T_c adds a constant amount.

for eps in (1e-2, 1e-3, 1e-4, 1e-5):
    print(f"1 - T/T_c = {eps:.0e}: c/k_B = {model.heat_capacity_per_spin(Tc * (1 - eps)) / K_B:.4f}")

plt.show()
