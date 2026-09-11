r"""Abstract base class for a molecular partition function, and its thermodynamic-functions container.

A molecular partition function :math:`q(T)` (translational, rotational,
vibrational, ...; :mod:`chemistrykit.statmech.systems.partition_functions`)
is the single object every one of statistical thermodynamics'
macroscopic quantities is built from (McQuarrie, *Statistical
Mechanics*, Ch. 3-6; Atkins & de Paula, *Physical Chemistry*, 11th ed.,
Ch. 13-14):

.. math::

    U - U(0) = k_BT^2\left(\frac{\partial \ln q}{\partial T}\right)_V,
    \qquad
    A = -k_BT\ln q, \qquad S = \frac{U-A}{T}

:class:`PartitionFunction` declares this common interface
(``value``/``internal_energy``/``entropy``/``heat_capacity_v``); each
:mod:`~chemistrykit.statmech.systems.partition_functions` subclass
implements the specific closed-form result for its own mode (the
formulas genuinely differ in structure between translational,
rotational, and vibrational partition functions -- e.g. only the
translational mode needs the :math:`N!` indistinguishability correction
that gives the Sackur-Tetrode equation -- so, mirroring
:class:`chemistrykit.md.core.base_system.PairPotential`, the derived
thermodynamic functions are declared abstract here rather than derived
generically by numerical differentiation).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

from chemistrykit.constants import NA

__all__ = ["PartitionFunction", "ThermodynamicFunctions"]


class PartitionFunction(ABC):
    """Common interface for a single-mode molecular partition function.

    Every method accepts `N`, the number of independent, non-interacting
    copies of this mode contributing to an extensive thermodynamic
    quantity -- `N` moles' worth of molecules if `N` is left at its
    default of Avogadro's number (giving molar quantities directly), a
    single molecule's contribution if ``N=1``, or any other ensemble
    size.
    """

    @abstractmethod
    def value(self, T):
        """Return the partition function q(T) for a single molecule/mode.

        Parameters
        ----------
        T : float
            Absolute temperature, in K.

        Returns
        -------
        float
        """

    @abstractmethod
    def internal_energy(self, T, N: float = NA):
        r"""Return :math:`U(T)-U(0) = Nk_BT^2(\partial\ln q/\partial T)_V`.

        Parameters
        ----------
        T : float
            Absolute temperature, in K.
        N : float, default :data:`chemistrykit.constants.NA`
            Number of molecules/modes.

        Returns
        -------
        float
            Energy, in J (relative to this mode's own ground state).
        """

    @abstractmethod
    def entropy(self, T, N: float = NA):
        """Return the entropy S(T) contributed by this mode.

        Parameters
        ----------
        T : float
        N : float, default :data:`chemistrykit.constants.NA`

        Returns
        -------
        float
            Entropy, in J/K.
        """

    @abstractmethod
    def heat_capacity_v(self, T, N: float = NA):
        """Return the constant-volume heat capacity Cv(T) contributed by this mode.

        Parameters
        ----------
        T : float
        N : float, default :data:`chemistrykit.constants.NA`

        Returns
        -------
        float
            Heat capacity, in J/K.
        """

    def helmholtz_free_energy(self, T, N: float = NA):
        r"""Return :math:`A=U-TS` for this mode, from :meth:`internal_energy` and :meth:`entropy`.

        Parameters
        ----------
        T : float
        N : float, default :data:`chemistrykit.constants.NA`

        Returns
        -------
        float
            Helmholtz free energy, in J.
        """
        return self.internal_energy(T, N) - T * self.entropy(T, N)


@dataclass
class ThermodynamicFunctions:
    """A snapshot of one mode's (or a combined molecule's) thermodynamic functions at a given T.

    Returned by
    :meth:`chemistrykit.statmech.systems.partition_functions.IdealGasMolecule.thermodynamic_functions`.
    """

    T: float
    """float: Temperature, in K."""

    q: float
    """float: The (combined) partition function value at `T`."""

    U: float
    """float: Internal energy (relative to the ground state), in J."""

    S: float
    """float: Entropy, in J/K."""

    Cv: float
    """float: Constant-volume heat capacity, in J/K."""

    A: float
    """float: Helmholtz free energy, ``U - T*S``, in J."""
