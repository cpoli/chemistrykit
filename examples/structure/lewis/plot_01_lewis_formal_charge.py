r"""
Lewis's shared electron pair: choosing between Lewis structures by formal charge
================================================================================

G. N. Lewis (1916) described a covalent bond as a pair of electrons shared
by two atoms, and each atom as seeking an octet of eight valence electrons.
Often several octet-satisfying structures can be drawn for one molecule.
Formal charge ranks them: split every bonding pair evenly, count the
electrons each atom then owns, and compare with the free atom,

.. math::

    FC = V - N - \frac{B}{2}.

The preferred structure keeps formal charges small and puts any negative
charge on the more electronegative atom. This example builds competing
Lewis structures of :math:`\mathrm{CO_2}`, :math:`\mathrm{N_2O}` and the
cyanate ion :math:`\mathrm{OCN^-}` with
:class:`~chemistrykit.structure.systems.lewis.LewisStructure`, checks each
atom's octet, and compares their formal charges.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.periodic_table import electronegativity
from chemistrykit.structure.systems.lewis import LewisStructure

candidates = {
    "CO2": {
        "O=C=O": LewisStructure(["O", "C", "O"], {(0, 1): 2, (1, 2): 2}, {0: 2, 2: 2}),
        "O#C-O": LewisStructure(["O", "C", "O"], {(0, 1): 3, (1, 2): 1}, {0: 1, 2: 3}),
    },
    "N2O": {
        "N=N=O": LewisStructure(["N", "N", "O"], {(0, 1): 2, (1, 2): 2}, {0: 2, 2: 2}),
        "N#N-O": LewisStructure(["N", "N", "O"], {(0, 1): 3, (1, 2): 1}, {0: 1, 2: 3}),
        "N-N#O": LewisStructure(["N", "N", "O"], {(0, 1): 1, (1, 2): 3}, {0: 3, 2: 1}),
    },
    "OCN-": {
        "O=C=N": LewisStructure(["O", "C", "N"], {(0, 1): 2, (1, 2): 2}, {0: 2, 2: 2}),
        "O-C#N": LewisStructure(["O", "C", "N"], {(0, 1): 1, (1, 2): 3}, {0: 3, 2: 1}),
        "O#C-N": LewisStructure(["O", "C", "N"], {(0, 1): 3, (1, 2): 1}, {0: 1, 2: 3}),
    },
}


def octet_count(structure, atom):
    """Electrons around an atom in Lewis's picture: lone-pair electrons plus every shared pair."""
    shared = sum(2 * order for pair, order in structure.bond_orders.items() if atom in pair)
    return 2 * structure.lone_pairs.get(atom, 0) + shared


for molecule, structures in candidates.items():
    print(molecule)
    for label, s in structures.items():
        assert all(octet_count(s, a) == 8 for a in range(3)), label
        fc = s.formal_charges()
        print(f"  {label:6s} octets OK, formal charges {[f'{fc[a]:+.0f}' for a in range(3)]}, total {s.total_formal_charge():+.0f}")

# %%
# Every candidate's formal charges sum to the species' net charge (0 for
# CO2 and N2O, -1 for cyanate), the standard consistency check.
# The preferred structure has the smallest total of absolute formal
# charges. Ties are broken by putting negative formal charge on the more
# electronegative atom (the smaller sum of FC times electronegativity).
# For cyanate, O=C=N and O-C#N tie on the first rule, and O-C#N wins
# because its -1 sits on oxygen:


def score(structure):
    fc = structure.formal_charges()
    return (sum(abs(v) for v in fc.values()), sum(fc[a] * electronegativity(sym) for a, sym in enumerate(structure.symbols)))


for molecule, structures in candidates.items():
    best = min(structures, key=lambda label: score(structures[label]))
    print(f"{molecule}: preferred Lewis structure {best}  (sum |FC| = {score(structures[best])[0]:.0f})")
assert min(candidates["OCN-"], key=lambda label: score(candidates["OCN-"][label])) == "O-C#N"

# %%
# Formal charge on each atom for every candidate structure:

fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
for ax, (molecule, structures) in zip(axes, candidates.items()):
    width = 0.8 / len(structures)
    for k, (label, s) in enumerate(structures.items()):
        fc = s.formal_charges()
        ax.bar([a + (k - (len(structures) - 1) / 2) * width for a in range(3)], [fc[a] for a in range(3)], width=width, label=label)
    ax.set_xticks(range(3), [f"{sym}{a + 1}" for a, sym in enumerate(next(iter(structures.values())).symbols)])
    ax.axhline(0.0, color="gray", linewidth=0.8)
    ax.set_title(molecule)
    ax.legend(fontsize=8)
axes[0].set_ylabel("formal charge")
fig.suptitle("Lewis structures ranked by formal charge (# = triple bond)")
fig.tight_layout()
plt.show()
