r"""The rigid-rotor model of molecular rotation.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 7.5 (rotational
motion) and Ch. 12.2 (rotational spectroscopy), or Levine, *Quantum
Chemistry*, 7th ed., Ch. 6.4.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import HBAR
from chemistrykit.quantum.core.base_system import QuantumSystem

__all__ = ["RigidRotor"]


class RigidRotor(QuantumSystem):
    r"""A linear rigid rotor of moment of inertia `I`.

    .. math::

        E_J = J(J+1)\frac{\hbar^2}{2I}, \qquad g_J = 2J+1, \qquad J=0,1,2,\dots

    (Atkins & de Paula, *Physical Chemistry*, 11th ed., eq. 7.20-7.21.)
    Each level is :math:`(2J+1)`-fold degenerate (one state per value of
    the projection quantum number :math:`M_J=-J,\dots,J`), and allowed
    microwave (pure rotational) transitions follow the selection rule
    :math:`\Delta J=\pm1` (Atkins & de Paula, *Physical Chemistry*, 11th
    ed., Ch. 12.2), giving evenly spaced absorption lines at
    :math:`2B(J+1)` where :math:`B=\hbar^2/(2I)` is the rotational
    constant.

    Parameters
    ----------
    moment_of_inertia : float
        Moment of inertia `I`, in kg m^2.

    Examples
    --------
    >>> rotor = RigidRotor(moment_of_inertia=1.45e-46)  # ~HCl-like
    >>> rotor.degeneracy(0)
    1
    >>> rotor.degeneracy(2)
    5
    """

    def __init__(self, moment_of_inertia: float):
        if moment_of_inertia <= 0:
            raise ValueError("moment_of_inertia must be positive")
        self.moment_of_inertia = float(moment_of_inertia)

    @classmethod
    def from_diatomic(cls, mass1: float, mass2: float, bond_length: float) -> RigidRotor:
        r"""Build a :class:`RigidRotor` from a diatomic's atomic masses and bond length.

        Uses the reduced mass :math:`\mu=m_1m_2/(m_1+m_2)` and
        :math:`I=\mu r^2` (Atkins & de Paula, *Physical Chemistry*, 11th
        ed., eq. 7.19).

        Parameters
        ----------
        mass1, mass2 : float
            Atomic masses, in kg.
        bond_length : float
            Bond length `r`, in m.

        Returns
        -------
        RigidRotor

        Examples
        --------
        >>> import scipy.constants as sc
        >>> rotor = RigidRotor.from_diatomic(mass1=1.008 * sc.atomic_mass, mass2=34.97 * sc.atomic_mass, bond_length=127.5e-12)  # HCl
        >>> bool(1.0e-47 < rotor.moment_of_inertia < 5.0e-47)
        True
        """
        mu = mass1 * mass2 / (mass1 + mass2)
        return cls(moment_of_inertia=mu * bond_length**2)

    @property
    def rotational_constant(self) -> float:
        r"""float: The rotational constant :math:`B=\hbar^2/(2I)`, in J."""
        return HBAR**2 / (2.0 * self.moment_of_inertia)

    def energy(self, J):
        r"""Return :math:`E_J=J(J+1)\hbar^2/(2I)`.

        Parameters
        ----------
        J : int or array-like of int
            Rotational quantum number(s), each >= 0.

        Returns
        -------
        float or ndarray
            Energy, in J.
        """
        J = np.asarray(J, dtype=np.float64)
        if np.any(J < 0):
            raise ValueError("J must be >= 0")
        return J * (J + 1.0) * self.rotational_constant

    def degeneracy(self, J: int) -> int:
        """Return the degeneracy :math:`g_J=2J+1` of level `J`.

        Parameters
        ----------
        J : int

        Returns
        -------
        int
        """
        if J < 0:
            raise ValueError("J must be >= 0")
        return 2 * J + 1

    def transition_energy(self, J: int) -> float:
        r"""Return the :math:`J\to J+1` absorption transition energy, :math:`2B(J+1)`.

        Parameters
        ----------
        J : int
            Initial (lower) rotational quantum number, >= 0.

        Returns
        -------
        float
            Energy, in J.

        Examples
        --------
        Successive rotational-spectrum lines are evenly spaced by
        :math:`2B` -- the textbook rigid-rotor selection-rule result:

        >>> rotor = RigidRotor(moment_of_inertia=1.45e-46)
        >>> spacing1 = rotor.transition_energy(1) - rotor.transition_energy(0)
        >>> spacing2 = rotor.transition_energy(2) - rotor.transition_energy(1)
        >>> bool(abs(spacing1 - spacing2) < 1e-30)
        True
        """
        if J < 0:
            raise ValueError("J must be >= 0")
        return float(self.energy(J + 1) - self.energy(J))
