r"""
Standard reduction potentials, redox-couple balancing, and cell potentials
=============================================================================

:data:`~chemistrykit.electrochem.systems.standard_potentials.STANDARD_REDUCTION_POTENTIALS`
is a small curated table of standard reduction potentials.
:func:`~chemistrykit.electrochem.systems.standard_potentials.cell_potential`
combines any two half-reactions into a standard cell potential, and
:func:`~chemistrykit.electrochem.systems.standard_potentials.balance_redox_reaction`
finds the electron-balancing multiples needed to write out the full
mass-balanced overall reaction.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.electrochem.systems.standard_potentials import (
    STANDARD_REDUCTION_POTENTIALS,
    balance_redox_reaction,
    cell_potential,
    is_spontaneous,
)

# %%
# The classic Daniell cell, Cu cathode / Zn anode.
cu = STANDARD_REDUCTION_POTENTIALS["Cu2+/Cu"]
zn = STANDARD_REDUCTION_POTENTIALS["Zn2+/Zn"]
E_daniell = cell_potential(cu, zn)
print(f"Daniell cell: E = {E_daniell:.2f} V, spontaneous = {is_spontaneous(E_daniell)}")

# %%
# Electron balancing: MnO4-/Mn2+ (5 e-) paired against Fe3+/Fe2+ (1 e-)
# needs the iron half-reaction scaled 5-fold to match electron counts,
# even though the cell *potential* itself needs no such scaling.
mno4 = STANDARD_REDUCTION_POTENTIALS["MnO4-/Mn2+"]
fe = STANDARD_REDUCTION_POTENTIALS["Fe3+/Fe2+"]
cathode_mult, anode_mult, n_total = balance_redox_reaction(mno4, fe)
E_permanganate = cell_potential(mno4, fe)
print(f"\nMnO4-/Fe2+ titration reaction: E = {E_permanganate:.2f} V")
print(f"Balancing: {cathode_mult}x MnO4- half-reaction + {anode_mult}x Fe2+ half-reaction, n_total = {n_total}")

# %%
# A bar chart of the tabulated reduction potentials, sorted -- the
# "electrochemical series". A half-reaction near the top (most negative)
# is a strong reducing agent in its reduced form; one near the bottom
# (most positive) is a strong oxidizing agent in its oxidized form.
names = list(STANDARD_REDUCTION_POTENTIALS.keys())
values = np.array([STANDARD_REDUCTION_POTENTIALS[name].E_standard for name in names])
order = np.argsort(values)

fig, ax = plt.subplots(figsize=(6, 7))
ax.barh(range(len(names)), values[order], color=np.where(values[order] >= 0, "crimson", "steelblue"))
ax.set_yticks(range(len(names)))
ax.set_yticklabels([names[i] for i in order], fontsize=8)
ax.axvline(0.0, color="black", linewidth=0.8)
ax.set_xlabel(r"$E^\circ$ (V vs. SHE)")
ax.set_title("The electrochemical series")
fig.tight_layout()
plt.show()
