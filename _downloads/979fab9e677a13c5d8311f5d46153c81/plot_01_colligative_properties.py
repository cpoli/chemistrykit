r"""
Raoult's law P-x-y diagrams and colligative properties
========================================================

A P-x-y diagram for an ideal binary solution
(:class:`~chemistrykit.thermo.systems.mixtures.BinaryIdealSolution`),
showing how the vapor phase is enriched in the more volatile component,
plus the three classic colligative properties -- freezing-point
depression, boiling-point elevation, and osmotic pressure -- as a
function of solute concentration.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.thermo.systems.mixtures import (
    BinaryIdealSolution,
    boiling_point_elevation,
    freezing_point_depression,
    osmotic_pressure,
)

# Benzene (A) / toluene (B) at ~25 degC, roughly (kPa).
solution = BinaryIdealSolution(P_A_star=12.7, P_B_star=3.8)
x_A = np.linspace(0.0, 1.0, 200)
P_total = solution.total_pressure(x_A)
y_A = solution.vapor_composition(x_A)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(x_A, P_total, label="liquid (P-x)", color="steelblue")
ax.plot(y_A, P_total, label="vapor (P-y)", color="darkorange")
ax.set_xlabel("mole fraction benzene")
ax.set_ylabel("P (kPa)")
ax.set_title("P-x-y diagram for an ideal benzene/toluene mixture")
ax.legend()
fig.tight_layout()

# %%
# Colligative properties of an aqueous NaCl solution (van't Hoff factor
# i=2 for complete dissociation), as a function of molality/molarity:

b = np.linspace(0.0, 1.0, 100)
Kf_water, Kb_water = 1.86, 0.51
dTf = freezing_point_depression(Kf_water, b, i=2.0)
dTb = boiling_point_elevation(Kb_water, b, i=2.0)
Pi = osmotic_pressure(M=b * 1000.0, T=298.15, i=2.0)  # M in mol/m^3

fig2, axes = plt.subplots(1, 3, figsize=(14, 4))
axes[0].plot(b, dTf, color="steelblue")
axes[0].set_xlabel("molality (mol/kg)")
axes[0].set_ylabel("freezing-point depression (K)")
axes[1].plot(b, dTb, color="darkorange")
axes[1].set_xlabel("molality (mol/kg)")
axes[1].set_ylabel("boiling-point elevation (K)")
axes[2].plot(b, Pi / 1000.0, color="seagreen")
axes[2].set_xlabel("molarity (mol/L)")
axes[2].set_ylabel("osmotic pressure (kPa)")
fig2.suptitle("Colligative properties of aqueous NaCl (i=2)")
fig2.tight_layout()

plt.show()
