r"""
Kapustinskii equation: lattice energy without a crystal structure
===================================================================

Kapustinskii (1956) replaced the structure-specific Madelung constant with
a near-universal constant per ion, so
:class:`~chemistrykit.crystal.systems.lattice_energy.KapustinskiiLatticeEnergy`
needs only the number of ions per formula unit, the charges, and the
ionic radii. On rock-salt halides it tracks the structure-based
Born-Lande value; it also handles salts of other structures (fluorite
CaF2) where no Madelung constant is supplied at all.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.crystal.systems.lattice_energy import BornLandeLatticeEnergy, KapustinskiiLatticeEnergy
from chemistrykit.crystal.utils.reference_data import SHANNON_IONIC_RADII_PM as R
from chemistrykit.crystal.utils.reference_data import average_born_exponent

rock_salt = {
    # name: (cation, anion, r0 in pm, cation config, anion config)
    "LiF": ("Li+", "F-", 201.0, "He", "Ne"),
    "NaCl": ("Na+", "Cl-", 282.0, "Ne", "Ar"),
    "KCl": ("K+", "Cl-", 315.0, "Ar", "Ar"),
    "KBr": ("K+", "Br-", 330.0, "Ar", "Kr"),
    "RbI": ("Rb+", "I-", 367.0, "Kr", "Xe"),
}
kap, bl = [], []
for name, (cat, an, r0, c1, c2) in rock_salt.items():
    U_k = KapustinskiiLatticeEnergy(2, 1, 1, R[cat], R[an]).lattice_energy() / 1000.0
    U_b = BornLandeLatticeEnergy(1, 1, r0 * 1e-12, average_born_exponent(c1, c2)).lattice_energy() / 1000.0
    kap.append(U_k)
    bl.append(U_b)
    print(f"{name:5s} Kapustinskii {U_k:7.1f}   Born-Lande {U_b:7.1f} kJ/mol   ({(U_k - U_b) / U_b:+.1%})")
    assert abs(U_k - U_b) / abs(U_b) < 0.08

# %%
# Structures Kapustinskii needs no Madelung constant for: doubly charged
# MgO (rock salt) and three-ion CaF2 (fluorite). Approximate experimental
# values: MgO about -3800 kJ/mol, CaF2 about -2630 kJ/mol.
U_mgo = KapustinskiiLatticeEnergy(2, 2, 2, R["Mg2+"], R["O2-"]).lattice_energy() / 1000.0
U_caf2 = KapustinskiiLatticeEnergy(3, 2, 1, R["Ca2+"], R["F-"]).lattice_energy() / 1000.0
print(f"\nMgO  (Kapustinskii): {U_mgo:7.1f} kJ/mol")
print(f"CaF2 (Kapustinskii): {U_caf2:7.1f} kJ/mol")

# %%
fig, ax = plt.subplots(figsize=(5.5, 5))
ax.scatter(bl, kap)
for name, x, y in zip(rock_salt, bl, kap, strict=True):
    ax.annotate(name, (x, y), textcoords="offset points", xytext=(5, -10))
lims = [min(bl + kap) - 50, max(bl + kap) + 50]
ax.plot(lims, lims, "k--", lw=0.8, label="perfect agreement")
ax.set_xlabel("Born-Lande U (kJ/mol, needs structure)")
ax.set_ylabel("Kapustinskii U (kJ/mol, radii only)")
ax.set_title("Kapustinskii's structure-free estimate")
ax.legend()
plt.tight_layout()
plt.show()
