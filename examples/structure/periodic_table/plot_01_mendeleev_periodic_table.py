r"""
Mendeleev's periodic table: periodicity and the atomic-weight inversions
========================================================================

Mendeleev ordered the elements by atomic weight but let chemical
behavior overrule the order where the two disagreed. Most famously he
placed tellurium (127.6) before iodine (126.9) so that iodine fell with
the halogens. Moseley later showed that the true ordering quantity is
the atomic number :math:`Z`. Left: the standard atomic weights of all 118
elements in :data:`chemistrykit.periodic_table.PERIODIC_TABLE` against
:math:`Z`. The points where weight *decreases* as :math:`Z` increases are
exactly the "inversions" an atomic-weight ordering gets wrong. Ar/K,
Co/Ni and Te/I are the chemically real ones; those past uranium reflect
only the convention of listing the longest-lived isotope's mass number
for elements with no stable isotope. Right: the
periodicity itself, Pauling electronegativity rising across each period
and collapsing at each alkali metal.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.periodic_table import PAULING_ELECTRONEGATIVITY, PERIODIC_TABLE_BY_NUMBER, molar_mass

Z = np.arange(1, 119)
mass = np.array([PERIODIC_TABLE_BY_NUMBER[z].atomic_mass for z in Z])
inversions = [z for z in Z[:-1] if mass[z] < mass[z - 1]]

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
axes[0].plot(Z, mass, ".", color="steelblue")
for z in inversions:
    a, b = PERIODIC_TABLE_BY_NUMBER[z], PERIODIC_TABLE_BY_NUMBER[z + 1]
    axes[0].plot([z, z + 1], [a.atomic_mass, b.atomic_mass], "o-", color="firebrick")
    if z < 90:  # past U the "weights" are isotope mass numbers
        axes[0].annotate(f"{a.symbol}/{b.symbol}", (z, a.atomic_mass), textcoords="offset points", xytext=(-30, 8))
axes[0].set_xlabel("atomic number Z")
axes[0].set_ylabel("standard atomic weight (u)")
axes[0].set_title("Atomic weight vs. Z: inversions in red")

en = [(z, PAULING_ELECTRONEGATIVITY[PERIODIC_TABLE_BY_NUMBER[z].symbol]) for z in Z if PERIODIC_TABLE_BY_NUMBER[z].symbol in PAULING_ELECTRONEGATIVITY]
zs, chis = np.array(en).T
axes[1].plot(zs, chis, "o-", color="darkorange", markersize=3)
for alkali in (3, 11, 19, 37, 55, 87):
    axes[1].axvline(alkali, color="gray", linewidth=0.6)
axes[1].set_xlabel("atomic number Z")
axes[1].set_ylabel("Pauling electronegativity")
axes[1].set_title("Periodicity: each period starts at an alkali metal (lines)")
fig.tight_layout()

# %%
print("inversions:", ", ".join(f"{PERIODIC_TABLE_BY_NUMBER[z].symbol}/{PERIODIC_TABLE_BY_NUMBER[z + 1].symbol}" for z in inversions))

# %%
# With the full table, molar masses of heavy-element compounds come
# straight from their formulas:

for f in ("UF6", "K2[PtCl6]", "HgCl2", "Gd2O3", "Pb(NO3)2"):
    print(f"{f:10s} {molar_mass(f):8.2f} g/mol")

plt.show()
