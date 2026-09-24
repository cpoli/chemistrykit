r"""
The Clausius-Clapeyron vapor-pressure curve of water
====================================================

Water's liquid-vapor phase boundary from the integrated
Clausius-Clapeyron equation
(:class:`~chemistrykit.thermo.systems.phase_equilibria.ClausiusClapeyron`),
anchored at the normal boiling point. Plotting :math:`\ln P` against
:math:`1/T` gives a straight line of slope :math:`-\Delta H_{vap}/R`;
:meth:`~chemistrykit.thermo.systems.phase_equilibria.ClausiusClapeyron.from_two_points`
recovers :math:`\Delta H_{vap}` from two measured points, and
:meth:`~chemistrykit.thermo.systems.phase_equilibria.ClausiusClapeyron.boiling_point`
inverts the curve to find the boiling point at altitude.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.thermo.systems.phase_equilibria import ClausiusClapeyron

water = ClausiusClapeyron(delta_h_vap=40700.0, T_ref=373.15, P_ref=101325.0)

T = np.linspace(280.0, 373.15, 200)
P = water.pressure(T)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
axes[0].plot(T, P / 1000.0, color="steelblue")
axes[0].scatter([373.15], [101.325], color="crimson", zorder=3, label="normal boiling point")
axes[0].set_xlabel("T (K)")
axes[0].set_ylabel("vapor pressure (kPa)")
axes[0].set_title("Water's liquid-vapor phase boundary")
axes[0].legend()

axes[1].plot(1000.0 / T, np.log(P), color="steelblue")
axes[1].set_xlabel("1000 / T (1/K)")
axes[1].set_ylabel("ln(P / Pa)")
axes[1].set_title(r"Straight line of slope $-\Delta H_{vap}/R$")
fig.tight_layout()

# %%
# Two points on the curve fix :math:`\Delta H_{vap}`, and inverting the
# model gives the boiling point at a mountain-altitude pressure of about
# 70 kPa:

fit = ClausiusClapeyron.from_two_points(353.15, float(water.pressure(353.15)), 373.15, 101325.0)
print(f"Enthalpy of vaporization from two points: {fit.delta_h_vap / 1000.0:.2f} kJ/mol")
T_boil = water.boiling_point(70_000.0)
print(f"At 70 kPa, water boils at {T_boil:.2f} K ({T_boil - 273.15:.2f} degC)")

plt.show()
