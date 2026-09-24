r"""
Raman scattering vs. infrared absorption: the mutual exclusion rule for CO2
=============================================================================

Raman and Krishnan (1928) found that a small fraction of scattered light
is shifted in frequency by a molecule's vibrational wavenumbers. A mode
is Raman-active if it changes the molecule's *polarizability*. It is
infrared-active if it changes the *dipole moment*. For a
centrosymmetric molecule such as CO2, no mode is active in both spectra.

This example takes CO2's normal modes from
:class:`~chemistrykit.spectro.systems.vibrational.TriatomicNormalModes`
and evaluates each mode's dipole and polarizability derivatives with two
simple models. The dipole comes from fixed partial charges (O
:math:`-q`, C :math:`+2q`). The polarizability comes from bond
polarizabilities that grow with bond length. The result is an infrared
spectrum and a Raman (Stokes-shift) spectrum with no lines in common.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

from chemistrykit.spectro.core.base_system import Spectrum
from chemistrykit.spectro.systems.vibrational import TriatomicNormalModes
from chemistrykit.spectro.visualizers.spectro_plots import plot_broadened_spectrum

bond_length = 116.3e-12
co2 = TriatomicNormalModes.linear(
    mass_terminal=15.999 * sc.atomic_mass,
    mass_central=12.011 * sc.atomic_mass,
    bond_length=bond_length,
    k_r=1600.0,
    k_theta=0.85e-18,
)
result = co2.solve()
L = result.internal_coordinate_vectors  # columns: modes; rows: (r1, r2, bend)
names = ["bend", "symmetric stretch", "antisymmetric stretch"]

# %%
# Dipole derivative: with charges -q, +2q, -q, stretching the two bonds
# changes the axial dipole by q(dr1 - dr2). The linearized bend
# coordinate S=(x0-2x1+x2)/r changes the transverse dipole by -q r S.
# Polarizability derivative: the isotropic sum of bond polarizabilities
# changes by alpha'(dr1 + dr2), while bending changes it only at second
# order.

q, alpha_prime = 1.0, 1.0  # arbitrary units: only zero vs. nonzero matters
ir_activity = (q * (L[0] - L[1])) ** 2 + (q * bond_length * L[2]) ** 2
raman_activity = (alpha_prime * (L[0] + L[1])) ** 2
ir_activity /= ir_activity.max()
raman_activity /= raman_activity.max()

for name, nu, ir, ra in zip(names, result.wavenumbers, ir_activity, raman_activity, strict=True):
    print(f"{name:>22s}: {nu:7.1f} cm^-1  IR {'active' if ir > 1e-6 else 'silent':>6s}  Raman {'active' if ra > 1e-6 else 'silent':>6s}")

# %%
# No mode appears in both spectra: the symmetric stretch is seen only in
# Raman, and the bend and antisymmetric stretch only in the infrared.

ir = Spectrum(positions=result.wavenumbers, intensities=np.where(ir_activity > 1e-6, ir_activity, 0.0))
raman = Spectrum(positions=result.wavenumbers, intensities=np.where(raman_activity > 1e-6, raman_activity, 0.0))
x = np.linspace(400.0, 2700.0, 3000)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
plot_broadened_spectrum(ir, x, ax=ax1, shape="lorentzian", fwhm=25.0, color="crimson")
ax1.set_title("Infrared absorption (dipole change)")
plot_broadened_spectrum(raman, x, ax=ax2, shape="lorentzian", fwhm=25.0, color="darkgreen")
ax2.set_title("Raman scattering (polarizability change), Stokes shift")
ax2.set_xlabel(r"wavenumber / Raman shift (cm$^{-1}$)")
for ax in (ax1, ax2):
    for nu, name in zip(result.wavenumbers, names, strict=True):
        ax.axvline(nu, color="0.7", lw=0.6, ls=":")
        ax.text(nu - 15.0, 0.45, name, transform=ax.get_xaxis_transform(), rotation=90, fontsize=8, ha="right", va="center")
fig.suptitle("CO2: mutual exclusion of Raman and infrared activity")
fig.tight_layout()
plt.show()
