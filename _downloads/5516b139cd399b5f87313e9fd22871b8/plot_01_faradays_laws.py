r"""
Faraday's laws of electrolysis
================================

Faraday's first law: the mass liberated at an electrode is proportional
to the charge passed, :math:`m = QM/(nF)`. His second law: for a fixed
charge, different substances are liberated in proportion to their
equivalent weights :math:`M/n`. This example verifies both with
:func:`~chemistrykit.electrochem.systems.electrolysis.faradays_law_mass`
and :func:`~chemistrykit.electrochem.systems.electrolysis.mass_from_charge`,
recovering the single universal constant -- the Faraday -- that links
charge to chemical change.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.constants import FARADAY
from chemistrykit.electrochem.systems.electrolysis import faradays_law_mass, mass_from_charge

# %%
# First law: copper electrorefining, Cu2+ + 2e- -> Cu (M = 63.546 g/mol).
# Mass deposited in one hour is a straight line through the origin in
# the current.
currents = np.linspace(1.0, 20.0, 8)
mass_deposited = faradays_law_mass(currents, 3600.0, 63.546, n=2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(currents, mass_deposited, "o-")
ax1.set_xlabel("Current (A)")
ax1.set_ylabel("Mass Cu deposited in 1 h (g)")
ax1.set_title("First law: m proportional to Q")

# %%
# Second law: the same charge (here one faraday) through several
# electrolytes connected in series liberates each substance in proportion
# to its equivalent weight M/n.
species = {"H2": (2.016, 2), "O2": (31.998, 4), "Cu": (63.546, 2), "Ag": (107.868, 1), "Al": (26.982, 3)}
masses = {name: mass_from_charge(FARADAY, M, n) for name, (M, n) in species.items()}
equivalents = {name: M / n for name, (M, n) in species.items()}
for name in species:
    print(f"{name:3s}: mass per faraday = {masses[name]:7.3f} g, equivalent weight M/n = {equivalents[name]:7.3f} g/mol")

ax2.scatter(list(equivalents.values()), list(masses.values()))
for name in species:
    ax2.annotate(name, (equivalents[name], masses[name]), textcoords="offset points", xytext=(5, -10))
ax2.set_xlabel("Equivalent weight M/n (g/mol)")
ax2.set_ylabel("Mass liberated by 1 F (g)")
ax2.set_title("Second law: m proportional to M/n")
fig.tight_layout()
plt.show()
