r"""
Madelung's lattice sum: the electrostatic energy of rock salt
===============================================================

Madelung (1918) asked for the electrostatic energy of one ion in the
field of every other ion of an infinite crystal. For rock salt it is
:math:`U=-Me^2/(4\pi\varepsilon_0r_0)` per ion pair, with the Madelung
constant :math:`M=\sum_j \pm r_0/r_j` summed shell by shell: 6 opposite
charges at :math:`r_0`, 12 like charges at :math:`\sqrt2r_0`, 8 opposite
charges at :math:`\sqrt3r_0`, and so on.
:func:`~chemistrykit.crystal.systems.madelung.madelung_constant_nacl`
returns the converged value.
"""

# %%
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.constants import ELEMENTARY_CHARGE, NA, VACUUM_PERMITTIVITY
from chemistrykit.crystal.systems.madelung import madelung_constant_nacl

# Group the ions around a central Na+ into coordination shells (r in units of r0).
shells = defaultdict(lambda: [0, 0])
n = 4
for i in range(-n, n + 1):
    for j in range(-n, n + 1):
        for k in range(-n, n + 1):
            r2 = i * i + j * j + k * k
            if 0 < r2 <= 12:
                shells[r2][0] += 1
                shells[r2][1] = 1 if (i + j + k) % 2 else -1  # +1 attractive (opposite charge)

print(f"{'shell':>5s} {'r/r0':>7s} {'ions':>5s} {'charge':>8s} {'term':>9s}")
r_values, terms = [], []
for r2 in sorted(shells):
    count, sign = shells[r2]
    term = sign * count / np.sqrt(r2)
    r_values.append(np.sqrt(r2))
    terms.append(term)
    print(f"{len(terms):5d} {np.sqrt(r2):7.3f} {count:5d} {'opposite' if sign > 0 else 'like':>8s} {term:+9.4f}")

# %%
# The shell terms are large and alternate in sign; only the full lattice
# sum gives Madelung's constant, and with it the Coulomb energy of the
# crystal:
M = madelung_constant_nacl()
r0 = 282e-12
U_pair = -M * ELEMENTARY_CHARGE**2 / (4.0 * np.pi * VACUUM_PERMITTIVITY * r0)
print(f"\nMadelung constant of NaCl: M = {M:.4f}")
print(f"Electrostatic energy per ion pair: {U_pair / ELEMENTARY_CHARGE:.2f} eV")
print(f"                          per mole: {U_pair * NA / 1000.0:.0f} kJ/mol")
print(f"(an isolated Na+Cl- pair at the same distance: {-1.0 * ELEMENTARY_CHARGE / (4.0 * np.pi * VACUUM_PERMITTIVITY * r0):.2f} eV)")
assert 1.74 < M < 1.75

# %%
fig, ax = plt.subplots(figsize=(8, 4))
colors = ["C0" if t > 0 else "C3" for t in terms]
ax.bar(r_values, terms, width=0.08, color=colors)
ax.axhline(M, color="k", ls="--", lw=1, label=f"Madelung constant M = {M:.4f}")
ax.axhline(0, color="gray", lw=0.5)
ax.set_xlabel(r"shell radius $r/r_0$")
ax.set_ylabel(r"shell term $\pm N_{shell}\,r_0/r$")
ax.set_title("Madelung's lattice sum for NaCl, shell by shell (blue: opposite charge)")
ax.legend()
plt.tight_layout()
plt.show()
