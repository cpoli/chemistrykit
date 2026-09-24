r"""
Baeyer's strain theory: angle strain in planar cycloalkane rings
================================================================

Adolf von Baeyer (1885) combined van't Hoff's tetrahedral carbon with the
assumption that rings are flat. A planar ring of :math:`n` carbons has
interior angles :math:`180^\circ(n-2)/n`, so each bond must bend away
from the tetrahedral 109.47 degrees by

.. math::

    \delta(n) = \tfrac{1}{2}\left[109.47^\circ - \frac{180^\circ\,(n-2)}{n}\right].

The theory correctly explains why three- and four-membered rings are
strained and five- and six-membered rings are common. It wrongly predicts
growing strain for larger rings: real rings pucker. This example computes
Baeyer's :math:`\delta(n)` with
:func:`~chemistrykit.structure.systems.ring_strain.baeyer_angle_strain`,
compares it with measured strain energies, and builds a chair cyclohexane
whose angles are exactly tetrahedral.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.structure.core.base_system import Molecule, angle_between
from chemistrykit.structure.systems.ring_strain import TETRAHEDRAL_ANGLE, baeyer_angle_strain, chair_cyclohexane_coordinates, planar_ring_angle
from chemistrykit.structure.visualizers.structure_plots import plot_molecule_3d

sizes = np.arange(3, 9)
# Total ring strain energies of the cycloalkanes (CH2)n from heats of
# combustion, kJ/mol (standard textbook values).
measured_strain = {3: 115.0, 4: 110.0, 5: 26.0, 6: 0.0, 7: 26.0, 8: 40.0}

for n in sizes:
    d = baeyer_angle_strain(int(n))
    angle = planar_ring_angle(int(n))
    strain = measured_strain[int(n)]
    print(f"C{n}H{2 * n}: planar angle {angle:6.2f} deg, Baeyer strain {d:+6.2f} deg/bond, measured ring strain {strain:5.0f} kJ/mol")

# %%
# Baeyer's cyclohexane, if planar, would have 120 degree angles and a
# strain of -5.26 degrees per bond. The chair form removes it entirely:
# every C-C-C angle is tetrahedral.

chair = chair_cyclohexane_coordinates(bond_length=1.54)
angles = [angle_between(chair[k - 1] - chair[k], chair[(k + 1) % 6] - chair[k]) for k in range(6)]
print(f"\nchair cyclohexane C-C-C angles: {np.round(angles, 4)} (tetrahedral = {TETRAHEDRAL_ANGLE:.4f})")
assert np.allclose(angles, TETRAHEDRAL_ANGLE)

# %%
fig = plt.figure(figsize=(12, 4.5))
ax1 = fig.add_subplot(1, 2, 1)
ax1.bar(sizes - 0.2, [baeyer_angle_strain(int(n)) for n in sizes], width=0.4, color="C0", label="Baeyer angle strain (deg/bond)")
ax1.axhline(0.0, color="gray", linewidth=0.8)
ax1.set_xlabel("ring size n")
ax1.set_ylabel("Baeyer strain (degrees per bond)", color="C0")
ax1b = ax1.twinx()
ax1b.bar(sizes + 0.2, [measured_strain[int(n)] for n in sizes], width=0.4, color="C1", label="measured strain energy")
ax1b.set_ylabel("measured ring strain (kJ/mol)", color="C1")
ax1.set_title("Baeyer's planar-ring prediction vs experiment")

ax2 = fig.add_subplot(1, 2, 2, projection="3d")
ring = Molecule(symbols=["C"] * 6, coordinates=chair, bonds=[(k, (k + 1) % 6) for k in range(6)])
plot_molecule_3d(ring, ax=ax2)
ax2.set_title("Chair cyclohexane: puckered, no angle strain")
fig.tight_layout()
plt.show()
