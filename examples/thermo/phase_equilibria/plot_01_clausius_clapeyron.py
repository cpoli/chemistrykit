r"""
The Clausius-Clapeyron phase boundary and the Gibbs phase rule
================================================================

Water's liquid-vapor phase boundary from the (integrated)
Clausius-Clapeyron equation
(:class:`~chemistrykit.thermo.systems.phase_equilibria.ClausiusClapeyron`),
anchored at the normal boiling point, plus a check of the Gibbs phase
rule (:func:`~chemistrykit.thermo.systems.phase_equilibria.gibbs_phase_rule`)
at a few notable points on a pure substance's phase diagram.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.thermo.systems.phase_equilibria import ClausiusClapeyron, gibbs_phase_rule

water = ClausiusClapeyron(delta_h_vap=40700.0, T_ref=373.15, P_ref=101325.0)

T = np.linspace(280.0, 373.15, 200)
P = water.pressure(T)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(T, P / 1000.0, color="steelblue")
ax.scatter([373.15], [101.325], color="crimson", zorder=3, label="normal boiling point")
ax.axhline(101.325, color="gray", linestyle=":", linewidth=0.8)
ax.set_xlabel("T (K)")
ax.set_ylabel("vapor pressure (kPa)")
ax.set_title("Water's liquid-vapor phase boundary (Clausius-Clapeyron)")
ax.legend()
fig.tight_layout()

# %%
# Inverting the model gives the boiling point at any pressure -- e.g. at
# a mountain-altitude pressure of about 70 kPa, water boils several
# degrees below 100 degC:

P_altitude = 70_000.0
T_boil = water.boiling_point(P_altitude)
print(f"At {P_altitude / 1000.0:.0f} kPa, water boils at {T_boil:.2f} K ({T_boil - 273.15:.2f} degC)")

# %%
# The Gibbs phase rule, F = C - P + 2, gives the number of intensive
# degrees of freedom at a few notable points on a pure substance's phase
# diagram:

print("Pure substance, single phase (F):", gibbs_phase_rule(n_components=1, n_phases=1))
print("Pure substance, two phases coexisting (F):", gibbs_phase_rule(n_components=1, n_phases=2))
print("Pure substance, triple point (F):", gibbs_phase_rule(n_components=1, n_phases=3))

plt.show()
