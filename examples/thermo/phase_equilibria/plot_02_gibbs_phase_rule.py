r"""
Gibbs's phase rule on the phase diagram of water
================================================

Gibbs's phase rule, :math:`F = C - P + 2`
(:func:`~chemistrykit.thermo.systems.phase_equilibria.gibbs_phase_rule`),
says how many intensive variables can change while a set of phases
coexists. On a one-component phase diagram, single-phase regions are
areas (:math:`F=2`), two-phase boundaries are curves (:math:`F=1`), and
three phases meet only at a point (:math:`F=0`). The sketch below draws
water's diagram around its triple point (vaporization and sublimation
curves from the Clausius-Clapeyron equation, melting line from the
Clapeyron slope) and labels each feature with its degrees of freedom.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.thermo.systems.phase_equilibria import ClausiusClapeyron, gibbs_phase_rule

T_tp, P_tp = 273.16, 611.657  # water triple point
vaporization = ClausiusClapeyron(delta_h_vap=45050.0, T_ref=T_tp, P_ref=P_tp)
sublimation = ClausiusClapeyron(delta_h_vap=51060.0, T_ref=T_tp, P_ref=P_tp)
# Clapeyron slope of the melting line: dP/dT = dH_fus / (T dV_fus), negative for water.
dH_fus, dV_fus = 6010.0, -1.63e-6
slope = dH_fus / (T_tp * dV_fus)

T_vap = np.linspace(T_tp, 300.0, 100)
T_sub = np.linspace(250.0, T_tp, 100)
P_melt = np.logspace(np.log10(P_tp), 7, 100)
T_melt = T_tp + (P_melt - P_tp) / slope

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(T_vap, vaporization.pressure(T_vap), label="liquid-vapor (F = 1)")
ax.plot(T_sub, sublimation.pressure(T_sub), label="solid-vapor (F = 1)")
ax.plot(T_melt, P_melt, label="solid-liquid (F = 1)")
ax.plot([T_tp], [P_tp], "ko", label="triple point (F = 0)")
ax.set_yscale("log")
ax.set_xlim(250.0, 300.0)
ax.set_ylim(50.0, 1e7)
for x, y, text in [(258.0, 1e5, "ice"), (288.0, 1e5, "liquid"), (288.0, 200.0, "vapor")]:
    ax.text(x, y, f"{text}\nF = {gibbs_phase_rule(1, 1)}", ha="center")
ax.set_xlabel("T (K)")
ax.set_ylabel("P (Pa)")
ax.set_title("Water near its triple point: F = C - P + 2")
ax.legend(loc="lower right")
fig.tight_layout()

# %%
# The rule for one and two components, and with a reaction constraint
# (e.g. :math:`N_2O_4 \rightleftharpoons 2NO_2` in the gas phase, where
# the two species count as only one independent component):

for C in (1, 2):
    print(f"C = {C}: " + ", ".join(f"P = {P}: F = {gibbs_phase_rule(C, P)}" for P in range(1, C + 3)))
print("N2O4/NO2 gas with the reaction at equilibrium: F =", gibbs_phase_rule(n_components=2, n_phases=1, reactions=1))

plt.show()
