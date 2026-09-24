r"""
Coulson's molecular-orbital bond order in butadiene, hexatriene, and benzene
============================================================================

Charles Coulson (1939) defined the pi bond order between neighbouring
atoms from the Huckel molecular orbitals themselves,

.. math::

    p_{ij} = \sum_k n_k\, c_{ik} c_{jk},

summing over occupied orbitals :math:`k` with :math:`n_k` electrons. It
gives fractional orders that no single Lewis structure shows. In
butadiene the terminal bonds come out at 0.894 and the central bond at
0.447, so the "single" bond has real double-bond character; in benzene
every bond has exactly 2/3. This example computes these values with
:func:`~chemistrykit.structure.systems.bonding.coulson_pi_bond_order`
from :class:`~chemistrykit.quantum.systems.huckel.HuckelSystem` orbitals.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.quantum.systems.huckel import HuckelSystem
from chemistrykit.structure.systems.bonding import coulson_pi_bond_order


def pi_bond_orders(system, n_atoms, bonds):
    """Coulson bond order of each bond, with the n_atoms pi electrons filling the lowest orbitals in pairs."""
    result = system.solve()
    order = np.argsort(result.energies)
    occupations = np.zeros(n_atoms)
    occupations[order[: n_atoms // 2]] = 2.0
    return [coulson_pi_bond_order(result.coefficients, occupations, i, j) for i, j in bonds]


molecules = {
    "ethene": (HuckelSystem.linear_polyene(2), 2, [(0, 1)]),
    "butadiene": (HuckelSystem.linear_polyene(4), 4, [(i, i + 1) for i in range(3)]),
    "hexatriene": (HuckelSystem.linear_polyene(6), 6, [(i, i + 1) for i in range(5)]),
    "benzene": (HuckelSystem.cyclic_polyene(6), 6, [(i, (i + 1) % 6) for i in range(6)]),
}
orders = {}
for name, (system, n, bonds) in molecules.items():
    orders[name] = pi_bond_orders(system, n, bonds)
    print(f"{name:10s}: pi bond orders " + ", ".join(f"{p:.3f}" for p in orders[name]))

assert np.allclose(orders["butadiene"], [2 / np.sqrt(5), 1 / np.sqrt(5), 2 / np.sqrt(5)], atol=1e-3)
assert np.allclose(orders["benzene"], 2 / 3)

# %%
# The pi bond order adds to the one sigma bond, so benzene's total C-C
# bond order is :math:`1 + 2/3 = 5/3` for every ring bond. Bars show the
# pi order bond by bond:

fig, axes = plt.subplots(1, 4, figsize=(13, 3.5), sharey=True)
for ax, (name, p) in zip(axes, orders.items()):
    ax.bar([f"{i + 1}-{(i + 1) % len(p) + 1 if name == 'benzene' else i + 2}" for i in range(len(p))], p, color="C0")
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
    ax.set_title(name)
    ax.tick_params(axis="x", labelsize=8)
axes[0].set_ylabel(r"Coulson $\pi$ bond order $p_{ij}$")
fig.suptitle("Coulson bond orders from Huckel orbitals (dashed: a full double bond)")
fig.tight_layout()
plt.show()
