r"""
Classifying the 7 crystal systems
====================================

:func:`~chemistrykit.crystal.systems.crystal_systems.classify_crystal_system`
classifies a unit cell purely from the equalities/inequalities among its
six lattice parameters :math:`(a,b,c,\alpha,\beta,\gamma)`, and
:func:`~chemistrykit.crystal.systems.crystal_systems.unit_cell_volume`
gives the general-cell volume formula that reduces to familiar shortcuts
(:math:`abc`, :math:`\frac{\sqrt3}{2}a^2c`) in special cases.
"""

# %%
import numpy as np

from chemistrykit.crystal.systems.crystal_systems import classify_crystal_system, unit_cell_volume

cells = {
    "NaCl (rock salt)": (5.640, 5.640, 5.640, 90.0, 90.0, 90.0),
    "TiO2 (rutile)": (4.593, 4.593, 2.959, 90.0, 90.0, 90.0),
    "alpha-quartz": (4.913, 4.913, 5.405, 90.0, 90.0, 120.0),
    "calcite (rhombohedral setting)": (6.375, 6.375, 6.375, 46.08, 46.08, 46.08),
    "gypsum": (5.679, 15.202, 6.522, 90.0, 113.83, 90.0),
    "K2Cr2O7 (potassium dichromate)": (7.418, 7.318, 13.379, 98.6, 90.4, 95.5),
}

# %%
for name, (a, b, c, alpha, beta, gamma) in cells.items():
    system = classify_crystal_system(a, b, c, alpha, beta, gamma)
    V = unit_cell_volume(a, b, c, alpha, beta, gamma)
    print(f"{name:35s} -> {system:14s} V = {V:8.2f} A^3")

# %%
# A right-angle (orthorhombic/tetragonal/cubic) cell's volume is exactly
# a*b*c -- the general formula's simplest special case:
V_general = unit_cell_volume(4.593, 4.593, 2.959, 90.0, 90.0, 90.0)
V_shortcut = 4.593 * 4.593 * 2.959
print(f"\nRutile: general formula = {V_general:.4f}, a*b*c shortcut = {V_shortcut:.4f}")
assert np.isclose(V_general, V_shortcut)

# %%
# The hexagonal shortcut (sqrt(3)/2 * a^2 * c) matches the general
# formula too, for quartz's hexagonal cell:
a, c = 4.913, 5.405
V_general = unit_cell_volume(a, a, c, 90.0, 90.0, 120.0)
V_shortcut = (np.sqrt(3.0) / 2.0) * a**2 * c
print(f"Quartz:  general formula = {V_general:.4f}, hexagonal shortcut = {V_shortcut:.4f}")
assert np.isclose(V_general, V_shortcut)
