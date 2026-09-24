r"""
Huckel's 4n+2 aromaticity rule from the computed pi spectrum
==============================================================

Erich Huckel's 1931 pi-electron theory reduces a planar conjugated ring to
a small matrix eigenvalue problem: one 2p_z orbital per carbon, Coulomb
integral :math:`\alpha` on the diagonal, and resonance integral
:math:`\beta` only between bonded neighbours
(:class:`~chemistrykit.quantum.systems.huckel.HuckelSystem`). Filling the
resulting levels two electrons at a time shows why rings with
:math:`4n+2` pi electrons (benzene, 6) form a closed shell while rings with
:math:`4n` (cyclobutadiene, 4) leave two electrons unpaired in a
degenerate pair. :func:`~chemistrykit.quantum.systems.huckel.is_aromatic_by_huckel_rule`
checks the rule against the *computed* degeneracies, not only by counting
electrons.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.quantum.systems.huckel import HuckelSystem, is_aromatic_by_huckel_rule
from chemistrykit.quantum.visualizers.quantum_plots import plot_huckel_levels

# alpha = 0, beta = -1: energies are in units of |beta|, bonding levels negative.
benzene = HuckelSystem.cyclic_polyene(6)
cyclobutadiene = HuckelSystem.cyclic_polyene(4)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))
plot_huckel_levels(benzene.solve(), n_pi_electrons=6, ax=ax1)
ax1.set_title("Benzene, 6 pi electrons: closed shell")
plot_huckel_levels(cyclobutadiene.solve(), n_pi_electrons=4, ax=ax2)
ax2.set_title("Cyclobutadiene, 4 pi electrons: open shell")
for ax in (ax1, ax2):
    ax.set_ylabel(r"orbital energy $(E-\alpha)/|\beta|$")
fig.tight_layout()

# %%
# The same test applied to a family of rings and ring ions. The cyclopropenyl
# cation (2 electrons), cyclopentadienyl anion (6) and tropylium cation (6)
# are all 4n+2 and closed shell; cyclobutadiene and cyclooctatetraene (4n)
# are not.

cases = [
    ("cyclopropenyl cation", 3, 2),
    ("cyclobutadiene", 4, 4),
    ("cyclopentadienyl anion", 5, 6),
    ("benzene", 6, 6),
    ("tropylium cation", 7, 6),
    ("cyclooctatetraene", 8, 8),
    ("cyclooctatetraene dianion", 8, 10),
]
names, verdicts = [], []
for name, n_atoms, n_pi in cases:
    energies = HuckelSystem.cyclic_polyene(n_atoms).solve().energies
    aromatic = is_aromatic_by_huckel_rule(n_pi, energies)
    names.append(f"{name}\n({n_pi} pi e-)")
    verdicts.append(aromatic)
    print(f"{name:27s} ring={n_atoms} pi electrons={n_pi:2d} aromatic={aromatic}")

fig2, ax3 = plt.subplots(figsize=(9, 4))
ax3.bar(range(len(cases)), [1] * len(cases), color=["seagreen" if v else "lightgray" for v in verdicts])
ax3.set_xticks(range(len(cases)))
ax3.set_xticklabels(names, fontsize=8)
ax3.set_yticks([])
ax3.set_title("Huckel 4n+2 rule: green = aromatic (closed-shell 4n+2 ring)")
fig2.tight_layout()

plt.show()
