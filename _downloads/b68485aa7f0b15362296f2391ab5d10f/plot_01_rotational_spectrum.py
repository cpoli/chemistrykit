r"""
HCl rotational spectrum, relative intensities, and the H->D isotope shift
============================================================================

Builds an HCl-like rigid rotor from :mod:`chemistrykit.quantum`, predicts
its rotational-spectrum line positions and relative (Boltzmann-weighted)
intensities, and shows how substituting deuterium for hydrogen shifts
every line -- the standard way isotopic substitution is used to confirm
a rotational-spectrum assignment.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

from chemistrykit.quantum.systems.rigid_rotor import RigidRotor
from chemistrykit.spectro.systems.rotational import isotope_shift_ratio, rotational_line_wavenumbers, rotational_spectrum
from chemistrykit.spectro.visualizers.spectro_plots import plot_broadened_spectrum, plot_stick_spectrum

rotor_hcl = RigidRotor.from_diatomic(mass1=1.008 * sc.atomic_mass, mass2=34.97 * sc.atomic_mass, bond_length=127.5e-12)
rotor_dcl = RigidRotor.from_diatomic(mass1=2.014 * sc.atomic_mass, mass2=34.97 * sc.atomic_mass, bond_length=127.5e-12)

lines_hcl = rotational_line_wavenumbers(rotor_hcl, J_max=8)
lines_dcl = rotational_line_wavenumbers(rotor_dcl, J_max=8)
print(f"HCl lines (cm^-1): {np.round(lines_hcl, 2)}")
print(f"DCl lines (cm^-1): {np.round(lines_dcl, 2)}")

ratio = isotope_shift_ratio(mass1a=1.008 * sc.atomic_mass, mass2=34.97 * sc.atomic_mass, mass1b=2.014 * sc.atomic_mass)
print(f"B(DCl)/B(HCl) = {ratio:.4f} -- close to the reduced-mass-ratio prediction of ~0.5")

# %%
# The full spectrum, with relative intensities from the Boltzmann
# population of each initial J level:

spectrum = rotational_spectrum(rotor_hcl, J_max=15, temperature=300.0)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
plot_stick_spectrum(spectrum, ax=ax1, color="steelblue")
ax1.set_xlabel("wavenumber (cm^-1)")
ax1.set_title("HCl rotational spectrum (stick, 300 K)")

x = np.linspace(0, spectrum.positions[-1] + 20.0, 2000)
plot_broadened_spectrum(spectrum, x, ax=ax2, shape="lorentzian", fwhm=3.0, color="crimson")
ax2.set_xlabel("wavenumber (cm^-1)")
ax2.set_title("Lorentzian-broadened spectrum")

fig.tight_layout()
plt.show()
