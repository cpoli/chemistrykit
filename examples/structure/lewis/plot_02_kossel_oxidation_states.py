r"""
Kossel's ionic bond: oxidation states as the fully ionic limit
==============================================================

Walther Kossel (1916) described bonding as the complete transfer of
electrons from one atom to another, leaving ions with noble-gas
configurations held together by electrostatic attraction. Oxidation
states are the bookkeeping version of that picture: give every bonding
pair entirely to the more electronegative atom (sharing only between
identical atoms) and compare the electron count with the free atom.

For a salt such as NaCl or MgO this recovers the real ionic charges. For
covalent molecules it gives the charges the atoms *would* carry if the
bonds were fully ionic. The example below follows oxygen through
:math:`\mathrm{H_2O}`, :math:`\mathrm{H_2O_2}`, :math:`\mathrm{O_2}` and
:math:`\mathrm{OF_2}`, where it ranges from -2 to +2 depending on which
partner is more electronegative, using
:meth:`~chemistrykit.structure.systems.lewis.LewisStructure.oxidation_states`.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.structure.systems.lewis import LewisStructure

ionic = {
    "NaCl": LewisStructure(["Na", "Cl"], {(0, 1): 1}, {1: 3}),
    "KF": LewisStructure(["K", "F"], {(0, 1): 1}, {1: 3}),
    "MgO": LewisStructure(["Mg", "O"], {(0, 1): 2}, {1: 2}),
    "CaO": LewisStructure(["Ca", "O"], {(0, 1): 2}, {1: 2}),
}
print("Salts: the ionic limit gives each ion's real charge")
for name, s in ionic.items():
    ox = s.oxidation_states()
    print(f"  {name:5s}: " + ", ".join(f"{sym} {ox[i]:+.0f}" for i, sym in enumerate(s.symbols)))
assert ionic["NaCl"].oxidation_states() == {0: 1.0, 1: -1.0}
assert ionic["MgO"].oxidation_states() == {0: 2.0, 1: -2.0}

# %%
# Oxygen's oxidation state in four molecules. Bonds to hydrogen give the
# electrons to oxygen, an O-O bond splits them evenly, and bonds to
# fluorine (the only element more electronegative than oxygen) give them
# away:

oxygen_compounds = {
    "H2O": (LewisStructure(["O", "H", "H"], {(0, 1): 1, (0, 2): 1}, {0: 2}), 0),
    "H2O2": (LewisStructure(["H", "O", "O", "H"], {(0, 1): 1, (1, 2): 1, (2, 3): 1}, {1: 2, 2: 2}), 1),
    "O2": (LewisStructure(["O", "O"], {(0, 1): 2}, {0: 2, 1: 2}), 0),
    "OF2": (LewisStructure(["O", "F", "F"], {(0, 1): 1, (0, 2): 1}, {0: 2, 1: 3, 2: 3}), 0),
}
oxygen_state = {}
for name, (s, o_index) in oxygen_compounds.items():
    oxygen_state[name] = s.oxidation_states()[o_index]
    print(f"{name:5s}: O oxidation state {oxygen_state[name]:+.0f}, formal charge on O {s.formal_charges()[o_index]:+.0f}")
assert [oxygen_state[k] for k in oxygen_compounds] == [-2.0, -1.0, 0.0, 2.0]

# %%
# The formal charge on oxygen is zero in all four, because formal charge
# is Lewis's even-sharing limit; the oxidation state is Kossel's ionic
# limit, and only it follows the partner's electronegativity.

fig, ax = plt.subplots(figsize=(6, 4))
names = list(oxygen_state)
ax.bar(names, [oxygen_state[n] for n in names], color=["C0", "C1", "C2", "C3"])
ax.axhline(0.0, color="gray", linewidth=0.8)
ax.set_ylabel("oxidation state of O")
ax.set_title("Kossel's ionic limit: oxygen from -2 to +2")
fig.tight_layout()
plt.show()
