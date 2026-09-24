r"""
Margules activity coefficients: deviations from Raoult's law
============================================================

Margules wrote the activity coefficients of a binary liquid as power
series in mole fraction. The one-parameter form used by
:class:`~chemistrykit.thermo.systems.mixtures.MargulesSolution`,
:math:`\ln\gamma_1 = A x_2^2` and :math:`\ln\gamma_2 = A x_1^2`, turns
Raoult's straight lines into curves: :math:`A > 0` (unlike molecules
repel) bows the total pressure upward and can produce a maximum-pressure
azeotrope, while :math:`A < 0` bows it downward. The model still obeys
the Gibbs-Duhem relation, and each component tends to Henry's law at
infinite dilution.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.thermo.systems.mixtures import MargulesSolution

P1_star, P2_star = 30.0, 20.0  # kPa
x1 = np.linspace(0.0, 1.0, 300)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
for A, color in [(-1.5, "steelblue"), (0.0, "gray"), (1.5, "crimson")]:
    sol = MargulesSolution(A=A, P1_star=P1_star, P2_star=P2_star)
    g1, g2 = sol.activity_coefficients(x1)
    axes[0].plot(x1, np.log(g1), color=color, label=f"A = {A}: ln gamma_1")
    axes[0].plot(x1, np.log(g2), "--", color=color)
    axes[1].plot(x1, sol.total_pressure(x1), color=color, label=f"A = {A}")
    if A != 0.0:
        axes[2].plot(x1, sol.total_pressure(x1), color=color, label=f"P-x, A = {A}")
        axes[2].plot(sol.vapor_composition(x1), sol.total_pressure(x1), ":", color=color, label=f"P-y, A = {A}")
axes[0].set_xlabel(r"$x_1$")
axes[0].set_ylabel(r"ln $\gamma$ (solid: 1, dashed: 2)")
axes[0].set_title("Activity coefficients")
axes[0].legend(fontsize=8)
axes[1].set_xlabel(r"$x_1$")
axes[1].set_ylabel("total pressure (kPa)")
axes[1].set_title("Positive and negative deviations from Raoult")
axes[1].legend()
axes[2].set_xlabel(r"$x_1$, $y_1$")
axes[2].set_ylabel("total pressure (kPa)")
axes[2].set_title("Azeotropes where P-x and P-y touch")
axes[2].legend(fontsize=8)
fig.tight_layout()

# %%
# The azeotrope is where :math:`x_1 = y_1`, and the Henry's-law constant
# of component 1 at infinite dilution is :math:`P_1^* e^{A}`:

for A in (-1.5, 1.5):
    sol = MargulesSolution(A=A, P1_star=P1_star, P2_star=P2_star)
    diff = sol.vapor_composition(x1[1:-1]) - x1[1:-1]
    i = int(np.argmin(np.abs(diff)))
    print(f"A = {A:+.1f}: azeotrope near x1 = {x1[1:-1][i]:.3f}, P = {float(sol.total_pressure(x1[1:-1][i])):.2f} kPa; K_H(1) = {sol.henry_constant_1:.1f} kPa")

plt.show()
