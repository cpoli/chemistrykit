r"""
Nernst's heat theorem: reaction entropy vanishes at absolute zero
=================================================================

For a reaction between perfect crystals, Nernst found that
:math:`\Delta G` and :math:`\Delta H` approach each other with zero slope
as :math:`T \to 0`, so :math:`\Delta S = -\partial\Delta G/\partial T
\to 0`. With Debye :math:`T^3` heat capacities,
:math:`\Delta C_p = \alpha T^3`, the curves are
:math:`\Delta H = \Delta H_0 + \alpha T^4/4`,
:math:`\Delta S = \alpha T^3/3`, and
:math:`\Delta G = \Delta H_0 - \alpha T^4/12`. Planck's stronger form,
:math:`S \to 0` for every perfect crystal, gives absolute entropies from
calorimetry, :math:`S(T) = \int_0^T C_p/T'\,dT'`. Those absolute
entropies underlie the standard Gibbs energies of formation that
:func:`~chemistrykit.thermo.systems.equilibrium.gibbs_energy_of_mixture`
takes as input.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

from chemistrykit.constants import R
from chemistrykit.thermo.systems.equilibrium import gibbs_energy_of_mixture

dH0 = -2000.0  # J/mol at T = 0
alpha = 2.0e-4  # J/(mol K^4), Delta Cp = alpha T^3
T = np.linspace(0.0, 40.0, 400)
dH = dH0 + alpha * T**4 / 4.0
dG = dH0 - alpha * T**4 / 12.0
dS = alpha * T**3 / 3.0

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
axes[0].plot(T, dH, label=r"$\Delta H$")
axes[0].plot(T, dG, label=r"$\Delta G$")
axes[0].set_xlabel("T (K)")
axes[0].set_ylabel("J/mol")
axes[0].set_title(r"$\Delta G$ and $\Delta H$ meet with zero slope at T = 0")
axes[0].legend()
axes[1].plot(T, dS, color="crimson", label=r"$\Delta S = -\partial \Delta G / \partial T$")
axes[1].plot(T, -np.gradient(dG, T), "k:", label=r"numerical slope of $\Delta G$")
axes[1].set_xlabel("T (K)")
axes[1].set_ylabel(r"$\Delta S$ (J/(mol K))")
axes[1].set_title(r"Reaction entropy $\to 0$ as $T \to 0$")
axes[1].legend()
fig.tight_layout()

# %%
# Planck's third law: the absolute entropy of a Debye crystal
# (:math:`C_p = \frac{12\pi^4}{5} R (T/\theta_D)^3` at low T) from
# numerical integration of :math:`C_p/T`, compared with the closed form
# :math:`S = C_p/3`.

theta_D = 150.0
T_cal = np.linspace(1e-6, 12.0, 2000)
Cp = 12.0 * np.pi**4 / 5.0 * R * (T_cal / theta_D) ** 3
S_num = cumulative_trapezoid(Cp / T_cal, T_cal, initial=0.0)
print(f"S(12 K): integrated {S_num[-1]:.5f} J/(mol K), closed form Cp/3 = {Cp[-1] / 3:.5f} J/(mol K)")

# %%
# Third-law absolute entropies turn calorimetry into equilibrium data.
# For :math:`N_2 + 2O_2 \to N_2O_4(g)` at 298.15 K, the standard entropies
# :math:`S^\circ(N_2) = 191.61`, :math:`S^\circ(O_2) = 205.14` and
# :math:`S^\circ(N_2O_4) = 304.29` J/(mol K), with
# :math:`\Delta_fH^\circ = 9.16` kJ/mol, give the Gibbs energy of formation
# without measuring any equilibrium constant; it is the kind of input
# :func:`~chemistrykit.thermo.systems.equilibrium.gibbs_energy_of_mixture`
# expects.

dfS = 304.29 - 191.61 - 2 * 205.14
dfG = 9160.0 - 298.15 * dfS
print(f"Delta_f S = {dfS:.1f} J/(mol K); Delta_f G = {dfG / 1000:.2f} kJ/mol (tabulated: 97.89 kJ/mol)")
print(f"G of 1 mol N2O4 at 1 bar relative to its elements: {gibbs_energy_of_mixture([1.0], [dfG], 298.15) / 1000:.2f} kJ")

plt.show()
