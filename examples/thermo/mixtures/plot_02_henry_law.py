r"""
Henry's law: the dilute-solute limit of a real solution
=======================================================

Henry's law, :math:`P_B = x_B K_H`
(:func:`~chemistrykit.thermo.systems.mixtures.henry_law_pressure`),
says the vapor pressure of a dilute solute (equivalently, the amount of
a gas that dissolves at a given pressure) is proportional to its mole
fraction, with an empirical constant :math:`K_H` that differs from the
pure solute's vapor pressure. Below, an illustrative real solute's
partial pressure follows the Henry's-law line at high dilution and
leaves it as the solute becomes concentrated, where Raoult's law
(:func:`~chemistrykit.thermo.systems.mixtures.raoult_vapor_pressure`)
takes over.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.thermo.systems.mixtures import henry_law_pressure, raoult_vapor_pressure

P_B_star = 40.0  # kPa, pure-B vapor pressure
K_H = 120.0  # kPa, Henry's law constant (slope at infinite dilution)

x_B = np.linspace(0.0, 1.0, 200)
# The quadratic that starts with slope K_H at x_B = 0 and reaches P_B* at x_B = 1.
P_B_real = x_B * (K_H + (P_B_star - K_H) * x_B)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
axes[0].plot(x_B, P_B_real, color="black", linewidth=2, label="real (illustrative) partial pressure")
axes[0].plot(x_B, henry_law_pressure(x_B, K_H), "--", color="darkorange", label=r"Henry's law ($x_B K_H$)")
axes[0].plot(x_B, raoult_vapor_pressure(x_B, P_B_star), ":", color="steelblue", label=r"Raoult's law ($x_B P_B^*$)")
axes[0].set_ylim(0, 60)
axes[0].set_xlabel("mole fraction of solute B")
axes[0].set_ylabel("partial pressure of B (kPa)")
axes[0].set_title("Henry's law holds for the dilute solute")
axes[0].legend()

x_dilute = np.linspace(1e-4, 0.1, 100)
error = (henry_law_pressure(x_dilute, K_H) - x_dilute * (K_H + (P_B_star - K_H) * x_dilute)) / (x_dilute * (K_H + (P_B_star - K_H) * x_dilute))
axes[1].plot(x_dilute, 100 * error, color="darkorange")
axes[1].set_xlabel("mole fraction of solute B")
axes[1].set_ylabel("Henry's-law error (%)")
axes[1].set_title("The error vanishes as the solution becomes dilute")
fig.tight_layout()

# %%
# Read the other way, Henry's law gives the solubility of a gas: the
# dissolved mole fraction is :math:`x_B = P_B/K_H`, doubling when the gas
# pressure above the liquid doubles.

for P in (1.0, 2.0, 4.0):
    print(f"P_B = {P:4.1f} kPa -> dissolved x_B = {P / K_H:.4f}")

plt.show()
