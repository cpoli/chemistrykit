r"""
Volta's pile: stacking zinc-silver cells in series
====================================================

Volta's 1800 pile was a stack of identical zinc | brine | silver (or
copper) pairs. Each pair contributes the difference between the two
metals' reduction potentials -- the "contact" driving force Volta
insisted on against Galvani's animal electricity -- and stacking pairs
in series adds their voltages. This example computes the single-pair
voltage with :func:`~chemistrykit.electrochem.systems.standard_potentials.cell_potential`,
shows how the pile voltage grows linearly with the number of pairs, and
compares Volta's zinc-silver pairing with zinc-copper and a frog-leg-free
"same metal" pair that produces nothing at all.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.electrochem.systems.standard_potentials import STANDARD_REDUCTION_POTENTIALS as T
from chemistrykit.electrochem.systems.standard_potentials import cell_potential

# %%
# One pair: silver is the cathode (more positive reduction potential),
# zinc the anode. Two pieces of the *same* metal give exactly zero --
# the dissimilarity of the metals is what matters.
pairs = {
    "Zn | Ag (Volta)": cell_potential(T["Ag+/Ag"], T["Zn2+/Zn"]),
    "Zn | Cu": cell_potential(T["Cu2+/Cu"], T["Zn2+/Zn"]),
    "Zn | Zn": cell_potential(T["Zn2+/Zn"], T["Zn2+/Zn"]),
}
for name, E in pairs.items():
    print(f"{name:16s}: {E:.2f} V per pair")

# %%
# Stacking N identical pairs in series multiplies the voltage by N.
n_pairs = np.arange(1, 61)
fig, ax = plt.subplots()
for name, E in pairs.items():
    ax.plot(n_pairs, n_pairs * E, label=name)
ax.set_xlabel("Number of pairs in the pile")
ax.set_ylabel("Open-circuit pile voltage (V)")
ax.set_title("Volta's pile: series stacking of two-metal cells")
ax.legend()
fig.tight_layout()
plt.show()
