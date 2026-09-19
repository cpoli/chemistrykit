r"""
Triatomic normal modes via the Wilson GF-matrix method: CO2 and H2O
======================================================================

A genuine (not tabulated) normal-mode calculation for a linear (CO2) and
a bent (H2O) triatomic, using
:class:`~chemistrykit.spectro.systems.vibrational.TriatomicNormalModes`:
a numerically built Wilson B-matrix, mass-weighted `G` matrix, and simple
diagonal valence-force-field `F` matrix are combined into a generalized
eigenvalue problem and solved for the true normal-mode wavenumbers.
"""

# %%
import scipy.constants as sc

from chemistrykit.spectro.systems.vibrational import TriatomicNormalModes

co2 = TriatomicNormalModes.linear(
    mass_terminal=15.999 * sc.atomic_mass,
    mass_central=12.011 * sc.atomic_mass,
    bond_length=116.3e-12,
    k_r=1600.0,
    k_theta=0.85e-18,
)
co2_result = co2.solve()
print("CO2 (linear) normal-mode wavenumbers (cm^-1):")
print(f"  bend (doubly degenerate): {co2_result.wavenumbers[0]:.1f}  (real: 667)")
print(f"  symmetric stretch (IR-inactive): {co2_result.wavenumbers[1]:.1f}  (real: 1388)")
print(f"  antisymmetric stretch: {co2_result.wavenumbers[2]:.1f}  (real: 2349)")

# %%
water = TriatomicNormalModes.bent(
    mass_terminal=1.008 * sc.atomic_mass,
    mass_central=15.999 * sc.atomic_mass,
    bond_length=95.8e-12,
    bond_angle_degrees=104.5,
    k_r=770.0,
    k_theta=0.7e-18,
)
water_result = water.solve()
print("\nH2O (bent) normal-mode wavenumbers (cm^-1):")
print(f"  bend: {water_result.wavenumbers[0]:.1f}  (real: 1595)")
print(f"  symmetric stretch: {water_result.wavenumbers[1]:.1f}  (real: 3657)")
print(f"  antisymmetric stretch: {water_result.wavenumbers[2]:.1f}  (real: 3756)")

# %%
# These come from a *simple diagonal* valence force field (no
# stretch-stretch or stretch-bend interaction constants) with
# literature-typical, but not precisely fitted, force constants -- so
# they land in the right ballpark and reproduce the correct qualitative
# ordering (bend well below both stretches) without matching experiment
# to spectroscopic precision. Stiffening the bend force constant raises
# only the bend frequency, as expected for a diagonal force field:

soft_bend = TriatomicNormalModes.bent(
    mass_terminal=1.008 * sc.atomic_mass,
    mass_central=15.999 * sc.atomic_mass,
    bond_length=95.8e-12,
    bond_angle_degrees=104.5,
    k_r=770.0,
    k_theta=0.3e-18,
).solve()
stiff_bend = TriatomicNormalModes.bent(
    mass_terminal=1.008 * sc.atomic_mass,
    mass_central=15.999 * sc.atomic_mass,
    bond_length=95.8e-12,
    bond_angle_degrees=104.5,
    k_r=770.0,
    k_theta=1.2e-18,
).solve()
print(f"\nBend wavenumber with soft k_theta: {soft_bend.wavenumbers[0]:.1f} cm^-1")
print(f"Bend wavenumber with stiff k_theta: {stiff_bend.wavenumbers[0]:.1f} cm^-1")
print(f"Stretch wavenumbers unchanged: {soft_bend.wavenumbers[1:]} vs {stiff_bend.wavenumbers[1:]}")
