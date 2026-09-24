r"""
Guldberg and Waage's law of mass action: Q = K at equilibrium
=============================================================

For :math:`N_2O_4 \rightleftharpoons 2NO_2` the reaction quotient
:math:`Q = x_{NO_2}^2/x_{N_2O_4}` (at :math:`P = P^\circ`), computed by
:func:`~chemistrykit.thermo.systems.equilibrium.reaction_quotient`,
rises monotonically with the extent of reaction. The law of mass action
says the mixture stops changing where :math:`Q` reaches the
equilibrium constant :math:`K`. Starting from different mixtures of
reactant and product, the equilibrium compositions found by
:func:`~chemistrykit.thermo.systems.equilibrium.solve_equilibrium_composition`
all give the same :math:`Q = K`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.constants import R
from chemistrykit.thermo.systems.equilibrium import reaction_quotient, solve_equilibrium_composition

T = 298.15
K = 4.0
gibbs_formation = [0.0, -R * T * np.log(K) / 2.0]  # [N2O4, NO2]
nu = np.array([-1.0, 2.0])

xi = np.linspace(0.01, 0.99, 300)
Q = [reaction_quotient([(1 - x) / (1 + x), 2 * x / (1 + x)], nu) for x in xi]

fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogy(xi, Q, color="steelblue", label=r"$Q(\xi)$ starting from pure $N_2O_4$")
ax.axhline(K, color="crimson", linestyle="--", label=f"K = {K}")
ax.set_xlabel(r"extent of reaction $\xi$ (mol)")
ax.set_ylabel("reaction quotient Q")
ax.set_title(r"$N_2O_4 \rightleftharpoons 2NO_2$: the reaction stops where Q = K")
ax.legend()
fig.tight_layout()

# %%
# Different starting mixtures end with different amounts of each gas, but
# always with the same value of the mass-action ratio (at fixed total
# pressure, that also pins the mole fractions):

for n0 in ([1.0, 0.0], [0.05, 1.9], [1.0, 1.0], [0.2, 3.0]):
    res = solve_equilibrium_composition(("N2O4", "NO2"), [[-1.0], [2.0]], n0, gibbs_formation, T)
    Q_eq = reaction_quotient(res.x, nu)
    print(f"start n = {n0}: equilibrium n = {np.round(res.n, 4)} mol, Q = {Q_eq:.4f}")

plt.show()
