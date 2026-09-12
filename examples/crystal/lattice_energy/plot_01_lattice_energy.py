r"""
Born-Lande vs. Kapustinskii lattice energy
==============================================

:class:`~chemistrykit.crystal.systems.lattice_energy.BornLandeLatticeEnergy`
and
:class:`~chemistrykit.crystal.systems.lattice_energy.KapustinskiiLatticeEnergy`
are two independently-derived estimates of the same physical quantity,
sharing the
:class:`~chemistrykit.crystal.core.base_system.LatticeEnergyModel`
interface: Born-Lande needs a Madelung constant and Born exponent (crystal-
structure-specific inputs), while Kapustinskii needs only ionic radii and
charges.
"""

# %%
from chemistrykit.crystal.systems.lattice_energy import BornLandeLatticeEnergy, KapustinskiiLatticeEnergy
from chemistrykit.crystal.systems.madelung import madelung_constant_nacl
from chemistrykit.crystal.utils.reference_data import SHANNON_IONIC_RADII_PM, average_born_exponent

M = madelung_constant_nacl()
n_born = average_born_exponent("Ne", "Ar")  # Na+ (Ne-like), Cl- (Ar-like)
r0 = 282e-12  # m, NaCl nearest-neighbor distance

born_lande = BornLandeLatticeEnergy(cation_charge=1, anion_charge=1, r0=r0, born_exponent=n_born, madelung_constant=M)
kapustinskii = KapustinskiiLatticeEnergy(
    n_ions=2,
    cation_charge=1,
    anion_charge=1,
    r_cation_pm=SHANNON_IONIC_RADII_PM["Na+"],
    r_anion_pm=SHANNON_IONIC_RADII_PM["Cl-"],
)

print(f"Madelung constant (NaCl, Evjen method): {M:.6f}")
print(f"Born exponent (Na+/Cl- average):        {n_born:.1f}")
print(f"Born-Lande lattice energy:    {born_lande.lattice_energy() / 1000.0:8.1f} kJ/mol")
print(f"Kapustinskii lattice energy:  {kapustinskii.lattice_energy() / 1000.0:8.1f} kJ/mol")
print("Experimental (Born-Haber cycle) reference: approximately -787 kJ/mol")

# %%
# Both approximations agree with each other, and with experiment, to
# within about 10%:
agreement = abs(born_lande.lattice_energy() - kapustinskii.lattice_energy()) / abs(born_lande.lattice_energy())
print(f"\nBorn-Lande vs. Kapustinskii relative difference: {agreement:.2%}")
assert agreement < 0.10

# %%
# A doubly-charged salt (MgO) is bound much more strongly, since lattice
# energy scales with the product of ionic charges:
mgo_kapustinskii = KapustinskiiLatticeEnergy(
    n_ions=2,
    cation_charge=2,
    anion_charge=2,
    r_cation_pm=SHANNON_IONIC_RADII_PM["Mg2+"],
    r_anion_pm=SHANNON_IONIC_RADII_PM["O2-"],
)
print(f"\nMgO (Kapustinskii): {mgo_kapustinskii.lattice_energy() / 1000.0:.1f} kJ/mol (compare NaCl above)")
