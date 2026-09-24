r"""
Bragg's law: indexing powder XRD peaks
========================================

Each line of a diffraction pattern is a family of :math:`(hkl)` planes
that satisfies Bragg's law :math:`n\lambda=2d\sin\theta`.
:func:`~chemistrykit.crystal.systems.xrd.powder_xrd_peaks` computes each
reflection's structure factor
(:func:`~chemistrykit.crystal.systems.xrd.structure_factor`) directly
from atomic positions and Miller indices, so BCC's and FCC's systematic
absences (reflections with zero intensity by symmetry) fall out
automatically rather than being hardcoded as a rule.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.crystal.systems.xrd import powder_xrd_peaks
from chemistrykit.crystal.visualizers.crystal_plots import plot_xrd_pattern

wavelength = 154.18  # Cu-Kalpha, pm

iron_bcc = powder_xrd_peaks("BCC", a=286.65, wavelength=wavelength, hkl_max=3)
copper_fcc = powder_xrd_peaks("FCC", a=361.5, wavelength=wavelength, hkl_max=3)
polonium_sc = powder_xrd_peaks("SC", a=336.0, wavelength=wavelength, hkl_max=2)

print("Fe (BCC) allowed reflections (h+k+l even only):")
for peak in iron_bcc:
    print(f"  {peak.hkl}  2theta={peak.two_theta:6.2f} deg  d={peak.d_spacing:6.2f} pm")

print("\nCu (FCC) allowed reflections (unmixed parity only):")
for peak in copper_fcc:
    print(f"  {peak.hkl}  2theta={peak.two_theta:6.2f} deg  d={peak.d_spacing:6.2f} pm")

print("\nalpha-Po (SC) allowed reflections (no absences):")
for peak in polonium_sc:
    print(f"  {peak.hkl}  2theta={peak.two_theta:6.2f} deg  d={peak.d_spacing:6.2f} pm")

# %%
# BCC's (100) reflection (h+k+l=1, odd) is absent -- the first line is
# (110) instead. FCC's (100) and (110) are both absent (mixed parity) --
# the first line is (111).
bcc_hkls = {p.hkl for p in iron_bcc}
fcc_hkls = {p.hkl for p in copper_fcc}
print(f"\n(100) present for BCC? {(1, 0, 0) in bcc_hkls}")
print(f"(110) present for BCC? {(1, 1, 0) in bcc_hkls}")
print(f"(100) present for FCC? {(1, 0, 0) in fcc_hkls}")
print(f"(111) present for FCC? {(1, 1, 1) in fcc_hkls}")

# %%
# Every listed peak obeys Bragg's law exactly -- recover the wavelength
# from each peak's d-spacing and angle, and the NaCl-style inverse problem
# (lattice constant from a measured angle) that the Braggs solved:

for peak in copper_fcc:
    lam = 2.0 * peak.d_spacing * np.sin(np.radians(peak.two_theta / 2.0))
    assert abs(lam - wavelength) < 1e-9
first = copper_fcc[0]
h, k, l = first.hkl
a_recovered = wavelength / (2.0 * np.sin(np.radians(first.two_theta / 2.0))) * np.sqrt(h * h + k * k + l * l)
print(f"\nCu lattice constant recovered from the (111) angle via Bragg's law: {a_recovered:.2f} pm")

# %%
fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
plot_xrd_pattern(iron_bcc, ax=axes[0])
axes[0].set_title("Fe (BCC)")
plot_xrd_pattern(copper_fcc, ax=axes[1])
axes[1].set_title("Cu (FCC)")
plot_xrd_pattern(polonium_sc, ax=axes[2])
axes[2].set_title("alpha-Po (SC)")
plt.tight_layout()
plt.show()
