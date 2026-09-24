r"""
Bethe's crystal-field splitting: d orbitals in an octahedral field
==================================================================

Hans Bethe asked in 1929 what happens to a free ion's degenerate levels
when it sits in a crystal, surrounded by charges with only the symmetry of
a point group. Group theory alone gives the answer. Under a rotation by
angle :math:`\alpha`, the :math:`2l+1` orbitals of angular momentum
:math:`l` have the character

.. math::

    \chi_l(\alpha) = \frac{\sin\left[(l + \tfrac{1}{2})\alpha\right]}{\sin(\alpha/2)},
    \qquad \chi_l(0) = 2l+1,

and, because d orbitals are even under inversion, an improper operation
:math:`iR` has the same character as :math:`R`. Reducing this
representation with the :math:`O_h` character table gives which
degeneracies survive. Five d orbitals split into :math:`e_g` (2) +
:math:`t_{2g}` (3): the central result of crystal-field theory.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.structure.systems.point_group import get_character_table

oh = get_character_table("Oh")

# Rotation angle of each O_h class, and whether it is improper (i times a
# rotation). The improper classes are i = i*E, S4 = i*C4, S6 = i*C3,
# sigma_h = i*C2, sigma_d = i*C2'.
proper_angle = {
    "E": 0.0,
    "8C3": 2 * np.pi / 3,
    "6C2": np.pi,
    "6C4": np.pi / 2,
    "3C2(=C4^2)": np.pi,
    "i": 0.0,
    "6S4": np.pi / 2,
    "8S6": 2 * np.pi / 3,
    "3sigma_h": np.pi,
    "6sigma_d": np.pi,
}
improper = {"i", "6S4", "8S6", "3sigma_h", "6sigma_d"}


def orbital_characters(l):
    """Characters of the (2l+1) orbitals of angular momentum l over the O_h classes."""
    parity = (-1) ** l  # inversion multiplies an l-orbital by (-1)^l
    chars = []
    for op in oh.operations:
        a = proper_angle[op]
        chi = 2 * l + 1 if np.isclose(a, 0.0) else np.sin((l + 0.5) * a) / np.sin(a / 2)
        chars.append(round(chi * (parity if op in improper else 1)))
    return chars


for l, name in [(1, "p"), (2, "d"), (3, "f")]:
    chars = orbital_characters(l)
    print(f"{name} orbitals (l={l}): characters {chars}")
    terms = [(str(m) if m > 1 else "") + irrep for irrep, m in oh.reduce(chars).items()]
    print("    -> " + " + ".join(terms))

d_split = oh.reduce(orbital_characters(2))
assert d_split == {"Eg": 1, "T2g": 1}

# %%
# The energy-level diagram. The splitting :math:`\Delta_o` is not fixed by
# symmetry, but the barycentre rule is: the 2 :math:`e_g` orbitals rise by
# :math:`\tfrac{3}{5}\Delta_o` and the 3 :math:`t_{2g}` orbitals fall by
# :math:`\tfrac{2}{5}\Delta_o`, so the average energy is unchanged.

levels = {"Eg": 0.6, "T2g": -0.4}
dims = {irrep: int(oh.character(irrep, "E")) for irrep in d_split}
assert np.isclose(sum(dims[k] * levels[k] for k in dims), 0.0)

fig, ax = plt.subplots(figsize=(6, 4.5))
for k in range(5):
    ax.hlines(0.0, 0.1 + 0.16 * k, 0.24 + 0.16 * k, color="gray", linewidth=3)
ax.text(0.47, 0.06, "free ion: 5 d orbitals", ha="center")
for irrep, x0 in [("Eg", 1.4), ("T2g", 1.3)]:
    for k in range(dims[irrep]):
        ax.hlines(levels[irrep], x0 + 0.2 * k, x0 + 0.14 + 0.2 * k, color="C0" if irrep == "Eg" else "C3", linewidth=3)
        ax.plot([0.9, x0], [0.0, levels[irrep]], color="lightgray", linestyle=":", linewidth=1)
ax.text(1.95, 0.6, r"$e_g$ ($+\frac{3}{5}\Delta_o$)", va="center")
ax.text(1.95, -0.4, r"$t_{2g}$ ($-\frac{2}{5}\Delta_o$)", va="center")
ax.set_xlim(0, 2.6)
ax.set_ylim(-0.7, 0.9)
ax.set_xticks([])
ax.set_ylabel(r"energy / $\Delta_o$")
ax.set_title(r"Bethe (1929): d orbitals in an $O_h$ field")
fig.tight_layout()
plt.show()
