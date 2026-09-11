"""Abstract base class for equation-of-state (EOS) models.

An equation of state relates pressure, molar volume, and temperature for
a pure substance -- ``P = P(Vm, T)`` or, inverted, ``Vm = Vm(P, T)``.
Every EOS in :mod:`chemistrykit.thermo.systems.equations_of_state`
(:class:`IdealGas`, :class:`VanDerWaals`, :class:`RedlichKwong`)
implements this common interface so they can be swapped in and compared
directly (e.g. plotting P-V isotherms for all three side by side),
mirroring the role
:class:`chemistrykit.kinetics.core.base_system.RateLaw` plays for the
integrated rate laws -- a small ABC capturing the one shape of model this
domain's ``systems/equations_of_state`` module has, the same way
``chemistrykit.kinetics`` has two ABCs (``RateLaw``, ``ReactionNetwork``)
for its two shapes of model. Reaction/phase equilibrium and mixture
properties elsewhere in :mod:`chemistrykit.thermo` don't share a common
polymorphic interface the way the three EOS models do, so (following
``chemistrykit.kinetics.systems.arrhenius``/``enzyme``'s precedent) they
are implemented directly as functions and small result dataclasses in
their own ``systems/`` modules rather than forced into an ABC here.

Solving :math:`V_m(P, T)` for a cubic EOS (van der Waals, Redlich-Kwong)
requires finding roots of a cubic polynomial in :math:`V_m`; below the
critical temperature there can be up to three real positive roots
(liquid, unstable, and vapor branches) -- see the ``branch`` parameter on
each concrete class's ``molar_volume``, and
:func:`chemistrykit.thermo.utils.cubic_roots.real_positive_roots` for the
shared root-selection numerics.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

from chemistrykit.constants import R as _R

__all__ = ["EquationOfState"]


class EquationOfState(ABC):
    """Common interface for a pure-substance PVT equation of state.

    Concrete subclasses implement :meth:`pressure` (an explicit formula)
    and :meth:`molar_volume` (generally a root-find, for any EOS cubic or
    higher in :math:`V_m`); :meth:`compressibility_factor` then follows
    for every subclass for free.
    """

    R: float = _R

    @abstractmethod
    def pressure(self, Vm, T):
        """Return the pressure at molar volume `Vm` and temperature `T`.

        Parameters
        ----------
        Vm : float or array-like of float
            Molar volume, in m^3/mol.
        T : float or array-like of float
            Absolute temperature, in K.

        Returns
        -------
        float or ndarray
            Pressure, in Pa.
        """

    @abstractmethod
    def molar_volume(self, P, T, **kwargs):
        """Return the molar volume at pressure `P` and temperature `T`.

        Parameters
        ----------
        P : float
            Pressure, in Pa.
        T : float
            Absolute temperature, in K.

        Returns
        -------
        float
            Molar volume, in m^3/mol.
        """

    def compressibility_factor(self, P, T, **kwargs):
        r"""Return the compressibility factor :math:`Z = PV_m/(RT)`.

        :math:`Z = 1` exactly for an ideal gas; deviations from 1
        quantify real-gas non-ideality (Atkins & de Paula, *Physical
        Chemistry*, 11th ed., Ch. 1.3).

        Parameters
        ----------
        P : float
            Pressure, in Pa.
        T : float
            Absolute temperature, in K.
        **kwargs
            Forwarded to :meth:`molar_volume` (e.g. ``branch`` for a
            cubic EOS).

        Returns
        -------
        float
        """
        Vm = self.molar_volume(P, T, **kwargs)
        return P * Vm / (self.R * T)
