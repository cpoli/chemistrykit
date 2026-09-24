r"""
Born-Lande equation: lattice energies of the alkali halides
============================================================

Born and Lande (1918) balanced the Madelung attraction :math:`-A/r`
against a short-range repulsion :math:`B/r^n`; setting :math:`dU/dr=0` at
the observed spacing :math:`r_0` gives
:math:`U=-\frac{N_AM|z_+z_-|e^2}{4\pi\varepsilon_0r_0}(1-\frac{1}{n})`,
which
:class:`~chemistrykit.crystal.systems.lattice_energy.BornLandeLatticeEnergy`
evaluates. This example shows the energy curve for NaCl and compares the
equation with Born-Haber values for six rock-salt halides.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.constants import ELEMENTARY_CHARGE, NA, VACUUM_PERMITTIVITY
from chemistrykit.crystal.systems.lattice_energy import BornLandeLatticeEnergy
from chemistrykit.crystal.systems.madelung import madelung_constant_nacl
from chemistrykit.crystal.utils.reference_data import average_born_exponent

M = madelung_constant_nacl()
r0 = 282e-12  # NaCl nearest-neighbor distance, m
n = average_born_exponent("Ne", "Ar")
A = NA * M * ELEMENTARY_CHARGE**2 / (4.0 * np.pi * VACUUM_PERMITTIVITY)  # J m / mol
B = A * r0 ** (n - 1) / n  # chosen so that dU/dr = 0 at r0

r = np.linspace(200e-12, 600e-12, 400)
U_attr = -A / r
U_rep = B / r**n
U = U_attr + U_rep
U_born_lande = BornLandeLatticeEnergy(1, 1, r0=r0, born_exponent=n).lattice_energy()
print(f"Minimum of the Born-Lande curve:  {U.min() / 1000.0:8.1f} kJ/mol at r = {r[np.argmin(U)] * 1e12:.0f} pm")
print(f"Born-Lande equation at r0=282 pm: {U_born_lande / 1000.0:8.1f} kJ/mol")
assert np.isclose(-A / r0 + B / r0**n, U_born_lande)

# %%
# Born-Lande vs. experiment (approximate Born-Haber lattice energies, kJ/mol):
salts = {
    # name: (r0 in pm, cation config, anion config, experimental U)
    "LiF": (201.0, "He", "Ne", -1037.0),
    "NaCl": (282.0, "Ne", "Ar", -787.0),
    "NaBr": (299.0, "Ne", "Kr", -747.0),
    "KCl": (315.0, "Ar", "Ar", -715.0),
    "KBr": (330.0, "Ar", "Kr", -682.0),
    "KI": (353.0, "Ar", "Xe", -649.0),
}
names, calc, expt = [], [], []
for name, (r0_pm, c1, c2, U_exp) in salts.items():
    model = BornLandeLatticeEnergy(1, 1, r0=r0_pm * 1e-12, born_exponent=average_born_exponent(c1, c2))
    U_calc = model.lattice_energy() / 1000.0
    names.append(name)
    calc.append(U_calc)
    expt.append(U_exp)
    print(f"{name:5s} Born-Lande {U_calc:7.1f}  experiment {U_exp:7.1f}  ({(U_calc - U_exp) / U_exp:+.1%})")
    assert abs(U_calc - U_exp) / abs(U_exp) < 0.05

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(r * 1e12, U_attr / 1000.0, "--", label=r"Madelung attraction $-A/r$")
axes[0].plot(r * 1e12, U_rep / 1000.0, ":", label=r"Born repulsion $B/r^n$")
axes[0].plot(r * 1e12, U / 1000.0, "k-", label="total")
axes[0].axvline(282.0, color="gray", lw=0.8)
axes[0].set_ylim(-1500, 1000)
axes[0].set_xlabel("r (pm)")
axes[0].set_ylabel("U (kJ/mol)")
axes[0].set_title("NaCl: Born-Lande energy curve")
axes[0].legend()
x = np.arange(len(names))
axes[1].bar(x - 0.2, np.abs(calc), 0.4, label="Born-Lande")
axes[1].bar(x + 0.2, np.abs(expt), 0.4, label="Born-Haber (expt.)")
axes[1].set_xticks(x, names)
axes[1].set_ylabel("|U| (kJ/mol)")
axes[1].set_title("Born-Lande vs. experiment")
axes[1].legend()
plt.tight_layout()
plt.show()
