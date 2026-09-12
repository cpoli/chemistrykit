r"""Hard-sphere packing of the simple cubic, BCC, FCC, and (ideal) HCP lattices.

Each model below places identical touching spheres of radius `r` on a
Bravais lattice of constant `a` and derives, from pure geometry, how many
spheres share a conventional unit cell, how many nearest neighbors each
sphere touches, and what fraction of space the spheres fill (the atomic
packing factor, APF). See Ashcroft & Mermin, *Solid State Physics*
(1976), Ch. 4, or West, *Solid State Chemistry and its Applications*, 2nd
ed. (2014), Ch. 1.2, for the standard derivations reproduced here.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.crystal.core.base_system import LatticePacking

__all__ = [
    "SimpleCubicPacking",
    "BodyCenteredCubicPacking",
    "FaceCenteredCubicPacking",
    "HexagonalClosePacking",
]


class SimpleCubicPacking(LatticePacking):
    r"""Simple cubic (SC) packing: one lattice point per cell, spheres touching along the cell edge.

    Touching along the edge gives :math:`r = a/2`, one sphere per cell,
    and 6 nearest neighbors (Ashcroft & Mermin, Ch. 4).

    Examples
    --------
    The exact, textbook atomic packing factor is :math:`\pi/6\approx0.5236`
    (the least efficient of the four lattices here):

    >>> bool(round(float(SimpleCubicPacking().packing_fraction()), 6) == round(np.pi / 6.0, 6))
    True
    """

    coordination_number = 6
    atoms_per_cell = 1.0

    def unit_cell_volume(self, a: float):
        a = np.asarray(a, dtype=np.float64)
        return a**3

    def atomic_radius(self, a: float):
        return np.asarray(a, dtype=np.float64) / 2.0


class BodyCenteredCubicPacking(LatticePacking):
    r"""Body-centered cubic (BCC) packing: 2 lattice points per cell, spheres touching along the body diagonal.

    The body diagonal has length :math:`\sqrt3 a` and spans 4 sphere
    radii (corner sphere - center sphere - opposite corner sphere), so
    :math:`r=\sqrt3a/4`; 2 atoms per conventional cell (8 corners x 1/8 +
    1 body center), 8 nearest neighbors (Ashcroft & Mermin, Ch. 4).

    Examples
    --------
    The exact atomic packing factor is :math:`\sqrt3\pi/8\approx0.6802`:

    >>> bool(round(float(BodyCenteredCubicPacking().packing_fraction()), 6) == round(float(np.sqrt(3.0) * np.pi / 8.0), 6))
    True
    """

    coordination_number = 8
    atoms_per_cell = 2.0

    def unit_cell_volume(self, a: float):
        a = np.asarray(a, dtype=np.float64)
        return a**3

    def atomic_radius(self, a: float):
        return np.sqrt(3.0) * np.asarray(a, dtype=np.float64) / 4.0


class FaceCenteredCubicPacking(LatticePacking):
    r"""Face-centered cubic (FCC) packing: 4 lattice points per cell, spheres touching along a face diagonal.

    The face diagonal has length :math:`\sqrt2a` and spans 4 sphere
    radii, so :math:`r=\sqrt2a/4`; 4 atoms per conventional cell (8
    corners x 1/8 + 6 faces x 1/2), 12 nearest neighbors -- the
    close-packed maximum for identical spheres (Ashcroft & Mermin, Ch. 4).

    Examples
    --------
    The exact atomic packing factor, :math:`\pi/(3\sqrt2)\approx0.7405`,
    is the well-known maximum density for identical spheres on a regular
    lattice:

    >>> bool(round(float(FaceCenteredCubicPacking().packing_fraction()), 6) == round(float(np.pi / (3.0 * np.sqrt(2.0))), 6))
    True
    """

    coordination_number = 12
    atoms_per_cell = 4.0

    def unit_cell_volume(self, a: float):
        a = np.asarray(a, dtype=np.float64)
        return a**3

    def atomic_radius(self, a: float):
        return np.sqrt(2.0) * np.asarray(a, dtype=np.float64) / 4.0


class HexagonalClosePacking(LatticePacking):
    r"""Hexagonal close-packed (HCP) packing on the primitive hexagonal cell.

    The primitive hexagonal cell (basal rhombus of side `a`, height `c`)
    has volume :math:`V=\frac{\sqrt3}{2}a^2c` and contains 2 lattice
    points (the ABAB stacking's 2-atom basis); each sphere touches 12
    nearest neighbors (6 in-plane, 3 in the layer above, 3 below), the
    same coordination as FCC -- both are close-packed arrangements of
    identical spheres, differing only in stacking sequence (ABAB... vs.
    ABCABC...), which is why they share the same maximum packing
    fraction.

    Parameters
    ----------
    c_over_a : float, optional
        The stacking-axis-to-basal-plane-edge ratio :math:`c/a`. Defaults
        to the *ideal* close-packing value :math:`\sqrt{8/3}\approx1.633`
        (Ashcroft & Mermin, Ch. 4), derived purely from geometry (equal
        touching spheres stacked ABAB); real HCP metals deviate from this
        (e.g. Mg: 1.624, Zn: 1.856 -- Zn's anomalously large ratio means
        its actual packing fraction is *below* the ideal value below, an
        approximation this model does not otherwise account for), so
        passing the real `c_over_a` gives a more accurate, but no longer
        exactly :math:`\pi/(3\sqrt2)`, packing fraction.

    Examples
    --------
    At the ideal :math:`c/a`, the packing fraction exactly matches FCC's
    :math:`\pi/(3\sqrt2)\approx0.7405` -- verified numerically here rather
    than assumed, since it is not obvious by inspection that two
    differently-stacked lattices share an atomic packing factor:

    >>> from chemistrykit.crystal.systems.packing import FaceCenteredCubicPacking
    >>> hcp_pf = HexagonalClosePacking().packing_fraction()
    >>> fcc_pf = FaceCenteredCubicPacking().packing_fraction()
    >>> round(float(hcp_pf), 9) == round(float(fcc_pf), 9)
    True

    A non-ideal ``c_over_a`` (e.g. zinc's 1.856) gives a *lower* packing
    fraction, since the spheres no longer touch snugly in the stacking
    direction:

    >>> HexagonalClosePacking(c_over_a=1.856).packing_fraction() < HexagonalClosePacking().packing_fraction()
    True
    """

    coordination_number = 12
    atoms_per_cell = 2.0

    def __init__(self, c_over_a: float | None = None):
        self.c_over_a = float(c_over_a) if c_over_a is not None else float(np.sqrt(8.0 / 3.0))

    def unit_cell_volume(self, a: float):
        a = np.asarray(a, dtype=np.float64)
        c = self.c_over_a * a
        return (np.sqrt(3.0) / 2.0) * a**2 * c

    def atomic_radius(self, a: float):
        return np.asarray(a, dtype=np.float64) / 2.0
