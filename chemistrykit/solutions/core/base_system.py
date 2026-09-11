"""Abstract base classes for solution-chemistry models, and result containers.

Two ABCs cover the two shapes of model in this subpackage, mirroring
:mod:`chemistrykit.kinetics.core.base_system`'s ``RateLaw``/``ReactionNetwork``
split:

* :class:`WeakElectrolyte` -- a single weak acid or weak base equilibrium
  with a closed-form (cubic, solved numerically) ion concentration
  (:mod:`chemistrykit.solutions.systems.acid_base`).
* :class:`Titration` -- a full pH-vs-titrant-volume curve, where each
  point is itself found by root-finding a charge-balance equation
  (:mod:`chemistrykit.solutions.systems.titration`). :meth:`Titration.curve`
  and :meth:`Titration.find_equivalence_point` are shared machinery
  provided here once for every concrete titration type, the same way
  ``ReactionNetwork.integrate`` is shared machinery for every concrete
  reaction network.

Solubility equilibria (:mod:`chemistrykit.solutions.systems.solubility`)
and activity-coefficient models
(:mod:`chemistrykit.solutions.systems.activity`) don't share a common
polymorphic interface the way the two families above do, so (following
``chemistrykit.kinetics.systems.arrhenius``/``enzyme``'s precedent) they
are implemented directly as functions in their own ``systems/`` modules.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

__all__ = ["WeakElectrolyte", "TitrationResult", "Titration"]


class WeakElectrolyte(ABC):
    r"""Common interface for a single weak acid or weak base dissociation equilibrium.

    A concrete subclass represents *either* a weak acid (in which case
    :meth:`_ion_concentration` returns :math:`[H^+]`) *or* a weak base
    (returning :math:`[OH^-]`) -- never both -- so
    :meth:`percent_dissociation`'s formula
    :math:`\alpha = ([\text{ion}] - K_w/[\text{ion}])/C` applies unchanged
    to either case (for an acid, :math:`[\text{ion}] - K_w/[\text{ion}] =
    [H^+] - [OH^-] = [A^-]`, the dissociated fraction's concentration; for
    a base, symmetrically, :math:`[OH^-] - [H^+] = [BH^+]`).
    """

    Kw: float = 1.0e-14
    """float: Water autoionization constant, default 1.0e-14 (25 degC)."""

    @property
    @abstractmethod
    def total_concentration(self) -> float:
        """float: Total (analytical) concentration of the weak electrolyte, in mol/L."""

    @property
    @abstractmethod
    def equilibrium_constant(self) -> float:
        """float: Ka (for a weak acid) or Kb (for a weak base)."""

    @abstractmethod
    def _ion_concentration(self) -> float:
        """Return :math:`[H^+]` (weak acid) or :math:`[OH^-]` (weak base), in mol/L."""

    def percent_dissociation(self) -> float:
        r"""Percentage of the weak electrolyte that has dissociated, :math:`100\alpha`.

        Returns
        -------
        float
        """
        ion = self._ion_concentration()
        alpha = (ion - self.Kw / ion) / self.total_concentration
        return 100.0 * alpha


@dataclass
class TitrationResult:
    """Container for the output of a :meth:`Titration.curve` call."""

    Vb: np.ndarray
    """ndarray: Volume(s) of titrant added, in L (or any consistent volume unit)."""

    pH: np.ndarray
    """ndarray: pH at each volume in `Vb`."""


class Titration(ABC):
    """Common base for a strong/weak acid-base titration-curve model.

    Concrete subclasses implement :meth:`pH_at`; :meth:`curve` and
    :meth:`find_equivalence_point` are then available on every subclass
    for free, mirroring how
    :meth:`chemistrykit.kinetics.core.base_system.ReactionNetwork.integrate`
    is shared machinery built once atop each subclass's ``rhs``.
    """

    @abstractmethod
    def pH_at(self, Vb: np.ndarray) -> np.ndarray:
        """Return the pH at each titrant volume in `Vb`.

        Parameters
        ----------
        Vb : ndarray
            Volume(s) of titrant added, in L.

        Returns
        -------
        ndarray
        """

    def curve(self, Vb) -> TitrationResult:
        """Compute the full titration curve over a range of titrant volumes.

        Parameters
        ----------
        Vb : array-like of float
            Volumes of titrant added, in L.

        Returns
        -------
        TitrationResult
        """
        Vb = np.atleast_1d(np.asarray(Vb, dtype=np.float64))
        pH = np.asarray(self.pH_at(Vb), dtype=np.float64)
        return TitrationResult(Vb=Vb, pH=pH)

    def find_equivalence_point(self, Vb) -> float:
        r"""Numerically locate the equivalence point as the point of steepest pH ascent.

        The equivalence point is where :math:`d(\text{pH})/dV_b` is
        largest in magnitude -- the inflection point of the titration
        curve (Harris, *Quantitative Chemical Analysis*, 9th ed., Ch.
        11). This is a generic numerical detector shared by every
        :class:`Titration` subclass; compare its result against a
        subclass's closed-form equivalence volume (stoichiometric
        ``Ca*Va/Cb``) as a correctness check.

        Parameters
        ----------
        Vb : array-like of float
            A sufficiently fine grid of titrant volumes spanning the
            equivalence point, in L.

        Returns
        -------
        float
            The titrant volume of steepest ascent.
        """
        result = self.curve(Vb)
        dpH_dV = np.gradient(result.pH, result.Vb)
        idx = int(np.argmax(np.abs(dpH_dV)))
        return float(result.Vb[idx])
