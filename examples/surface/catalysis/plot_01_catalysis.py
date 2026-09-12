r"""
Catalysis: turnover frequency and Arrhenius rate enhancement
================================================================

A catalyst lowers a reaction's activation energy without changing its
thermodynamics, which -- via the Arrhenius equation
(:mod:`chemistrykit.kinetics.systems.arrhenius`) -- can speed the
reaction up by many orders of magnitude even for a modest reduction in
:math:`E_a`. :func:`~chemistrykit.surface.systems.catalysis.compare_catalyzed_rate`
quantifies this, and :func:`~chemistrykit.surface.systems.catalysis.turnover_frequency`/
:func:`~chemistrykit.surface.systems.catalysis.turnover_number` are the
standard measures of catalytic activity and durability.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.surface.systems.catalysis import compare_catalyzed_rate, turnover_frequency, turnover_number

Ea_uncatalyzed = 100.0e3  # J/mol
T = 298.15  # K

# %%
# Scan a range of activation-energy reductions and see how sharply the
# rate enhancement grows.
delta_Ea_values = np.array([10e3, 20e3, 30e3, 40e3, 50e3])
enhancements = []
for delta_Ea in delta_Ea_values:
    comparison = compare_catalyzed_rate(Ea_uncatalyzed=Ea_uncatalyzed, Ea_catalyzed=Ea_uncatalyzed - delta_Ea, T=T, A_uncatalyzed=1e13)
    enhancements.append(comparison.rate_enhancement)
    print(f"Delta Ea = {delta_Ea / 1e3:5.1f} kJ/mol -> rate enhancement = {comparison.rate_enhancement:.3e}")

# %%
# Turnover frequency and turnover number for a specific catalytic system.
rate = 5.0e-3  # mol/L/s
active_sites = 2.0e-6  # mol/L
tof = turnover_frequency(rate, active_sites)
print(f"\nTurnover frequency: {tof:.1f} 1/s")

moles_converted, moles_catalyst = 0.5, 1.0e-4
ton = turnover_number(moles_converted, moles_catalyst)
print(f"Turnover number after full conversion: {ton:.0f} cycles")

# %%
fig, ax = plt.subplots()
ax.semilogy(delta_Ea_values / 1e3, enhancements, "o-")
ax.set_xlabel(r"$\Delta E_a$ (kJ/mol)")
ax.set_ylabel("rate enhancement $k_{cat}/k_{uncat}$")
ax.set_title("Catalytic rate enhancement vs. activation-energy reduction")
plt.tight_layout()
plt.show()
