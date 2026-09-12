r"""
Formal charge and oxidation state from a Lewis structure
===========================================================

Computes formal charges and oxidation states -- two different, both
genuinely useful, electron-bookkeeping conventions -- for a handful of
classic examples: water, the ammonium and hydronium cations, hydrogen
peroxide, and carbon dioxide, using only connectivity, bond orders, lone
pairs, and Pauling electronegativities (from
:mod:`chemistrykit.periodic_table`).
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.structure.systems.lewis import LewisStructure

structures = {
    "water": LewisStructure(symbols=["O", "H", "H"], bond_orders={(0, 1): 1, (0, 2): 1}, lone_pairs={0: 2}),
    "ammonium (NH4+)": LewisStructure(symbols=["N", "H", "H", "H", "H"], bond_orders={(0, k): 1 for k in (1, 2, 3, 4)}),
    "hydronium (H3O+)": LewisStructure(symbols=["O", "H", "H", "H"], bond_orders={(0, k): 1 for k in (1, 2, 3)}, lone_pairs={0: 1}),
    "hydrogen peroxide": LewisStructure(symbols=["H", "O", "O", "H"], bond_orders={(0, 1): 1, (1, 2): 1, (2, 3): 1}, lone_pairs={1: 2, 2: 2}),
    "carbon dioxide": LewisStructure(symbols=["O", "C", "O"], bond_orders={(0, 1): 2, (1, 2): 2}, lone_pairs={0: 2, 2: 2}),
}

for name, structure in structures.items():
    charges = structure.formal_charges()
    states = structure.oxidation_states()
    print(f"{name}:")
    for i, symbol in enumerate(structure.symbols):
        print(f"    atom {i} ({symbol}): formal charge = {charges[i]:+.1f}, oxidation state = {states[i]:+.1f}")
    print(f"    total formal charge = {structure.total_formal_charge():+.1f}")

# %%
# A useful self-consistency check on any proposed Lewis structure: the
# formal charges must sum to the species' actual net charge (0 for water
# and hydrogen peroxide, +1 for ammonium and hydronium):

for name, expected_charge in [("water", 0.0), ("ammonium (NH4+)", 1.0), ("hydronium (H3O+)", 1.0), ("hydrogen peroxide", 0.0)]:
    total = structures[name].total_formal_charge()
    print(f"{name}: total formal charge {total:+.1f} (expected {expected_charge:+.1f}) -> {'OK' if abs(total - expected_charge) < 1e-9 else 'MISMATCH'}")

# %%
# Water and hydrogen peroxide's central-oxygen oxidation states -2 and -1
# bracket O2's 0 -- the textbook progression of oxygen oxidation states:

fig, ax = plt.subplots(figsize=(6, 4))
labels = ["water (O)", "H2O2 (O)"]
oxidation = [structures["water"].oxidation_states()[0], structures["hydrogen peroxide"].oxidation_states()[1]]
ax.bar(labels, oxidation, color=["steelblue", "darkorange"])
ax.axhline(0.0, color="gray", linewidth=0.8)
ax.set_ylabel("oxygen oxidation state")
ax.set_title("Oxygen oxidation state depends on its bonding environment")
fig.tight_layout()
plt.show()
