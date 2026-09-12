r"""Abstract base classes for crystal-structure models, and a note on why exactly these two ABCs (and no others).

This domain has two genuine "several interchangeable models of one
physical quantity, compared side by side" shapes -- the same shape that
:class:`chemistrykit.thermo.core.base_system.EquationOfState`,
:class:`chemistrykit.surface.core.base_system.AdsorptionIsotherm`, and
:class:`chemistrykit.polymer.core.base_system.PolymerChainModel` already
capture in their own domains:

* :class:`LatticePacking` -- the four common close-packed/simple hard-
  sphere lattice types (simple cubic, body-centered cubic, face-centered
  cubic, hexagonal close-packed;
  :mod:`chemistrykit.crystal.systems.packing`) each answer the same
  question -- "given identical touching spheres on this lattice, what
  fraction of space do they fill, and how many nearest neighbors does
  each have?" -- via a different atoms-per-cell count and
  radius-to-lattice-constant relation.
* :class:`LatticeEnergyModel` -- the Born-Lande and Kapustinskii
  equations (:mod:`chemistrykit.crystal.systems.lattice_energy`) are two
  independently-derived ways to estimate the same physical quantity (the
  lattice energy of an ionic solid) from different inputs (crystal
  structure + Born exponent, vs. only ionic radii and charges), so they
  are usefully compared side by side on the same interface.

Every other model in this domain is a self-contained closed-form
relationship with no swappable sibling, so -- following
``chemistrykit.electrochem``/``chemistrykit.photochem``/``chemistrykit.surface``'s
precedent of not forcing an ABC where none is warranted -- each stays a
plain function (plus, where useful, a small result dataclass) in its own
``systems/`` module:

* Crystal-system classification and unit-cell volume
  (:mod:`chemistrykit.crystal.systems.crystal_systems`) is a single
  classification function, not a family of swappable models.
* The Madelung constant (:mod:`chemistrykit.crystal.systems.madelung`) is
  a single numerical lattice summation, not several competing physical
  models of the same thing.
* Bragg's law and powder-XRD structure factors
  (:mod:`chemistrykit.crystal.systems.xrd`) are exact geometric/
  diffraction relationships.
* Schottky/Frenkel defect equilibrium
  (:mod:`chemistrykit.crystal.systems.defects`) are two distinct defect
  *mechanisms* (vacancy pairs vs. vacancy-interstitial pairs) with
  different formulas and different inputs (one vs. two site counts) --
  not two competing estimates of the same quantity -- so, like
  :mod:`chemistrykit.solutions.systems.solubility`'s Ksp/common-ion
  functions, they stay as plain functions rather than being forced onto
  a shared interface.

**On overlap with** :mod:`chemistrykit.structure.systems.point_group`:
that module detects a *molecule's* point-group symmetry by geometrically
testing candidate rotation/reflection operations against a finite set of
atomic coordinates. Classifying a *crystal* into one of the 7 crystal
systems, by contrast, is a lookup on the equalities/inequalities among
the six lattice parameters :math:`(a,b,c,\alpha,\beta,\gamma)` -- there is
no analogous "does this coordinate set map onto itself under this
operation" test to share, and a space group (the crystallographic
analogue of a point group, combining point-symmetry operations with
lattice translations) is a substantially larger object that this domain
does not attempt to compute. The two modules are conceptually related
(both ultimately classify structures by symmetry) but do not share
reusable code.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

__all__ = ["LatticePacking", "LatticeEnergyModel"]


class LatticePacking(ABC):
    """Common interface for a hard-sphere packing model of a simple Bravais lattice.

    Concrete subclasses (:class:`chemistrykit.crystal.systems.packing.SimpleCubicPacking`,
    :class:`~chemistrykit.crystal.systems.packing.BodyCenteredCubicPacking`,
    :class:`~chemistrykit.crystal.systems.packing.FaceCenteredCubicPacking`,
    :class:`~chemistrykit.crystal.systems.packing.HexagonalClosePacking`)
    set the class attributes `coordination_number` and `atoms_per_cell`,
    and implement :meth:`unit_cell_volume` and :meth:`atomic_radius`;
    :meth:`packing_fraction` is then available for every subclass for
    free, exactly as :meth:`chemistrykit.surface.core.base_system.AdsorptionIsotherm.fractional_coverage`
    is built once atop each concrete isotherm's ``loading``.
    """

    #: int: Number of nearest neighbors each sphere touches.
    coordination_number: int

    #: float: Number of lattice points (spheres) per conventional unit cell.
    atoms_per_cell: float

    @abstractmethod
    def unit_cell_volume(self, a: float) -> float:
        """Return the (conventional) unit-cell volume for lattice constant `a`.

        Parameters
        ----------
        a : float or array-like of float
            Cubic (or, for HCP, basal) lattice constant.

        Returns
        -------
        float or ndarray
        """

    @abstractmethod
    def atomic_radius(self, a: float):
        """Return the touching-sphere radius implied by lattice constant `a`.

        Parameters
        ----------
        a : float or array-like of float

        Returns
        -------
        float or ndarray
        """

    def packing_fraction(self, a: float = 1.0):
        r"""Return the fraction of the unit cell's volume occupied by touching spheres.

        .. math::

            \text{APF} = \frac{Z \cdot \frac{4}{3}\pi r^3}{V_{cell}}

        with `Z` = :attr:`atoms_per_cell` and `r` from :meth:`atomic_radius`
        -- independent of `a` for any of these lattices (both the sphere
        volume and the cell volume scale as :math:`a^3`), so the default
        `a=1.0` is only a placeholder unit.

        Parameters
        ----------
        a : float or array-like of float, default 1.0
            Lattice constant (cancels out of the result).

        Returns
        -------
        float or ndarray
        """
        a = np.asarray(a, dtype=np.float64)
        r = self.atomic_radius(a)
        sphere_volume = (4.0 / 3.0) * np.pi * r**3
        result = self.atoms_per_cell * sphere_volume / self.unit_cell_volume(a)
        return float(result) if result.ndim == 0 else result


class LatticeEnergyModel(ABC):
    """Common interface for a model estimating an ionic crystal's lattice energy.

    Concrete subclasses (:class:`chemistrykit.crystal.systems.lattice_energy.BornLandeLatticeEnergy`,
    :class:`~chemistrykit.crystal.systems.lattice_energy.KapustinskiiLatticeEnergy`)
    implement :meth:`lattice_energy`.
    """

    @abstractmethod
    def lattice_energy(self) -> float:
        """Return the (negative, exothermic-formation-convention) lattice energy.

        Returns
        -------
        float
            Lattice energy, in J/mol (negative for a stable ionic solid).
        """
