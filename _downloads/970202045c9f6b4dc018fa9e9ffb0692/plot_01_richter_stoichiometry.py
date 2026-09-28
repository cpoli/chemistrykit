r"""
Richter's stoichiometry: equivalent weights, balanced equations, limiting reagents
==================================================================================

Jeremias Richter, who coined the word *stoichiometry*, tabulated how many
parts by weight of each base neutralize a fixed weight of an acid, and
found the proportions were fixed and mutually consistent. Today those
numbers follow from a balanced equation and molar masses.
:func:`~chemistrykit.stoichiometry.balance_equation` finds the smallest
integer coefficients that conserve every element (and charge) as the
exact null space of the composition matrix. Left: Richter-style
neutralization equivalents, the grams of each base that neutralize
1000 g of sulfuric acid. Right: product yield as the N2:H2 feed ratio
varies. It rises while H2 is limiting and falls once N2 is, with the
kink at the stoichiometric 1:3 ratio
(:func:`~chemistrykit.stoichiometry.limiting_reagent`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.periodic_table import molar_mass
from chemistrykit.stoichiometry import balance_equation, limiting_reagent, theoretical_yield

for equation in ("C3H8 + O2 -> CO2 + H2O", "KMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2", "MnO4- + Fe2+ + H+ -> Mn2+ + Fe3+ + H2O"):
    print(balance_equation(equation))

bases = ["NaOH", "KOH", "Ca(OH)2", "Mg(OH)2", "NH3"]
products = {"NaOH": "Na2SO4 + H2O", "KOH": "K2SO4 + H2O", "Ca(OH)2": "CaSO4 + H2O", "Mg(OH)2": "MgSO4 + H2O", "NH3": "(NH4)2SO4"}
grams = []
for base in bases:
    eq = balance_equation(f"H2SO4 + {base} -> {products[base]}")
    grams.append(1000.0 / molar_mass("H2SO4") / eq.reactants["H2SO4"] * eq.reactants[base] * molar_mass(base))

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
axes[0].bar(bases, grams, color="steelblue")
axes[0].set_ylabel("g of base per 1000 g H2SO4")
axes[0].set_title("Neutralization equivalents, Richter-style")

total = 4.0  # mol of N2 + H2 fed
f_N2 = np.linspace(0.01, 0.99, 197)
nh3 = [theoretical_yield("N2 + H2 -> NH3", {"N2": f * total, "H2": (1 - f) * total}, "NH3") for f in f_N2]
axes[1].plot(f_N2, nh3, color="darkorange")
axes[1].axvline(0.25, color="gray", linestyle="--", label="stoichiometric N2:H2 = 1:3")
axes[1].set_xlabel("mole fraction of N2 in the feed")
axes[1].set_ylabel("theoretical NH3 (mol)")
axes[1].set_title("Yield is set by the limiting reagent")
axes[1].legend()
fig.tight_layout()

# %%
for f in (0.15, 0.20, 0.30, 0.40):
    moles = {"N2": f * total, "H2": (1 - f) * total}
    print(f"x(N2) = {f:.2f}: limiting reagent {limiting_reagent('N2 + H2 -> NH3', moles)}")

plt.show()
