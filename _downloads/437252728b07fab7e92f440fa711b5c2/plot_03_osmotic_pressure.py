r"""
van't Hoff's osmotic pressure law, Pi = iMRT
============================================

van't Hoff found that a dilute solution's osmotic pressure
(:func:`~chemistrykit.thermo.systems.mixtures.osmotic_pressure`) equals
the pressure the solute particles would exert as an ideal gas at the
same concentration and temperature. The factor :math:`i` counts
particles per formula unit: 1 for sucrose, about 2 for NaCl. The same
factor multiplies the other colligative effects, freezing-point
depression and boiling-point elevation
(:func:`~chemistrykit.thermo.systems.mixtures.freezing_point_depression`,
:func:`~chemistrykit.thermo.systems.mixtures.boiling_point_elevation`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.thermo.systems.equations_of_state import IdealGas
from chemistrykit.thermo.systems.mixtures import (
    CRYOSCOPIC_CONSTANTS,
    boiling_point_elevation,
    freezing_point_depression,
    osmotic_pressure,
)

T = 298.15
c = np.linspace(0.0, 0.5, 100)  # mol/L
M = c * 1000.0  # mol/m^3

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
axes[0].plot(c, osmotic_pressure(M, T, i=1.0) / 1e5, label="sucrose (i = 1)")
axes[0].plot(c, osmotic_pressure(M, T, i=2.0) / 1e5, label="NaCl (i = 2)")
Vm_gas = 1.0 / np.maximum(M, 1e-12)
axes[0].plot(c[::10], IdealGas().pressure(Vm_gas[::10], T) / 1e5, "k.", label="ideal gas at same n/V")
axes[0].set_xlabel("solute concentration (mol/L)")
axes[0].set_ylabel(r"osmotic pressure $\Pi$ (bar)")
axes[0].set_title("Solute particles behave like an ideal gas")
axes[0].legend()

water = CRYOSCOPIC_CONSTANTS["water"]
b = c  # dilute aqueous: molality close to molarity
axes[1].plot(b, freezing_point_depression(water["Kf"], b, i=2.0), label="freezing-point depression")
axes[1].plot(b, boiling_point_elevation(water["Kb"], b, i=2.0), label="boiling-point elevation")
axes[1].set_xlabel("NaCl molality (mol/kg)")
axes[1].set_ylabel(r"$\Delta T$ (K)")
axes[1].set_title("Other colligative effects share the factor i = 2")
axes[1].legend()
fig.tight_layout()

# %%
print(f"0.1 M sucrose at 25 degC: Pi = {osmotic_pressure(100.0, T) / 1e5:.3f} bar")
print(f"Isotonic saline (0.154 M NaCl, i = 2) at 37 degC: Pi = {osmotic_pressure(154.0, 310.15, i=2.0) / 1e5:.2f} bar")

plt.show()
