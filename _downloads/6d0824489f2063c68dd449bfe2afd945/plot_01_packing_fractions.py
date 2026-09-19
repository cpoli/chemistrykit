r"""
Packing efficiency of SC, BCC, FCC, and HCP
==============================================

The four lattices implemented here
(:mod:`chemistrykit.crystal.systems.packing`) share the
:class:`~chemistrykit.crystal.core.base_system.LatticePacking` interface:
each derives its atomic packing factor from pure geometry (touching
spheres, atoms-per-cell, and a lattice-constant-to-radius relation).
FCC and ideal HCP -- both close-packed, differing only in stacking
sequence -- turn out to share the same maximum packing fraction.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.crystal.systems.packing import (
    BodyCenteredCubicPacking,
    FaceCenteredCubicPacking,
    HexagonalClosePacking,
    SimpleCubicPacking,
)
from chemistrykit.crystal.visualizers.crystal_plots import plot_packing_fractions

lattices = {
    "SC": SimpleCubicPacking(),
    "BCC": BodyCenteredCubicPacking(),
    "FCC": FaceCenteredCubicPacking(),
    "HCP (ideal)": HexagonalClosePacking(),
}

for name, lattice in lattices.items():
    pf = lattice.packing_fraction()
    print(f"{name:12s} packing fraction = {pf:.4f}   coordination = {lattice.coordination_number:2d}   atoms/cell = {lattice.atoms_per_cell}")

# %%
# FCC and ideal HCP are both close-packed and share the same packing
# fraction, even though the underlying models (cubic vs. hexagonal cell,
# 4 vs. 2 atoms per cell) look nothing alike:
fcc_pf = lattices["FCC"].packing_fraction()
hcp_pf = lattices["HCP (ideal)"].packing_fraction()
print(f"\nFCC packing fraction:  {fcc_pf:.9f}")
print(f"HCP packing fraction:  {hcp_pf:.9f}")
print(f"Match: {abs(fcc_pf - hcp_pf) < 1e-9}")

# %%
# A real HCP metal's actual c/a ratio deviates from the ideal value,
# which lowers its packing fraction below the FCC/ideal-HCP maximum:
zinc_hcp = HexagonalClosePacking(c_over_a=1.856)  # zinc's anomalous c/a
print(f"\nZinc-like HCP (c/a=1.856) packing fraction: {zinc_hcp.packing_fraction():.4f} (below the ideal {hcp_pf:.4f})")

# %%
names = list(lattices.keys())
fractions = [lattices[name].packing_fraction() for name in names]
ax = plot_packing_fractions(names, fractions)
plt.tight_layout()
plt.show()
