r"""
Herzberg's linear-molecule rules: CO2's doubly degenerate bend and 3N-5 vibrations
===================================================================================

Herzberg's *Infrared and Raman Spectra of Polyatomic Molecules* (1945)
set out the rules for counting and labelling vibrations of polyatomic
molecules. A nonlinear molecule of N atoms has :math:`3N-6` vibrations.
A linear molecule has :math:`3N-5`, because rotation about its axis
does not count as a degree of freedom. For CO2 this means four vibrations: the
symmetric stretch, the antisymmetric stretch, and a *doubly degenerate*
bend (in-plane and out-of-plane bending have the same frequency by
cylindrical symmetry). This example compares linear CO2 with bent H2O
using
:class:`~chemistrykit.spectro.systems.vibrational.TriatomicNormalModes`,
whose ``is_linear`` flag marks the bend as a degenerate pair.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

from chemistrykit.spectro.systems.vibrational import TriatomicNormalModes

u = sc.atomic_mass
co2 = TriatomicNormalModes.linear(mass_terminal=15.999 * u, mass_central=12.011 * u, bond_length=116.3e-12, k_r=1600.0, k_theta=0.85e-18).solve()
h2o = TriatomicNormalModes.bent(
    mass_terminal=1.008 * u, mass_central=15.999 * u, bond_length=95.8e-12, bond_angle_degrees=104.5, k_r=770.0, k_theta=0.7e-18
).solve()

n_atoms = 3
for name, res in (("CO2", co2), ("H2O", h2o)):
    degeneracies = np.array([2 if (res.is_linear and k == 0) else 1 for k in range(3)])
    rule = "3N-5" if res.is_linear else "3N-6"
    expected = 3 * n_atoms - (5 if res.is_linear else 6)
    print(f"{name}: is_linear={res.is_linear}; vibrations counted with degeneracy = {degeneracies.sum()} ({rule} = {expected})")
    print(f"   wavenumbers (cm^-1): {np.round(res.wavenumbers, 1)}, degeneracies {degeneracies.tolist()}")

# %%
# In the linear molecule the bend coordinate is the transverse
# displacement :math:`(x_0-2x_1+x_2)/r`. The out-of-plane (y) bend gives
# the same frequency, so a stick for the bend carries twice the weight.

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax, name, res in ((axes[0], "CO2 (linear, 3N-5 = 4)", co2), (axes[1], "H2O (bent, 3N-6 = 3)", h2o)):
    degeneracies = [2 if (res.is_linear and k == 0) else 1 for k in range(3)]
    ax.vlines(res.wavenumbers, 0.0, degeneracies, color="tab:purple", lw=3)
    for nu, g in zip(res.wavenumbers, degeneracies, strict=True):
        ax.annotate("doubly degenerate bend" if g == 2 else "", (nu, g), xytext=(5, -15), textcoords="offset points", fontsize=8)
    ax.set_ylim(0.0, 2.5)
    ax.set_ylabel("number of modes")
    ax.set_xlabel(r"wavenumber (cm$^{-1}$)")
    ax.set_title(name)
fig.suptitle("Counting vibrations: a linear molecule's degenerate bend")
fig.tight_layout()
plt.show()
