"""Abstract base class for adsorption-isotherm models, and a note on why no other ABC is added here.

Langmuir, Freundlich, and BET are three interchangeable models of the
same physical quantity -- how much of a gas (or solute) is adsorbed on a
surface at equilibrium, as a function of pressure (or concentration) --
each with its own standard linearization for fitting parameters from
data (see each concrete ``systems/`` module). This is exactly the "several
swappable models of one physical relationship, compared side by side"
shape that :class:`chemistrykit.thermo.core.base_system.EquationOfState`
already captures for pure-substance PVT equations of state, and
:class:`chemistrykit.electrochem.core.base_system.BatteryDischargeModel`
captures for battery discharge curves, so :class:`AdsorptionIsotherm` is
this domain's one ABC.

Langmuir-Hinshelwood surface-reaction kinetics
(:mod:`chemistrykit.surface.systems.langmuir_hinshelwood`) *consumes* a
Langmuir coverage rather than providing an isotherm interface of its own,
and the turnover-number/rate-enhancement catalysis model
(:mod:`chemistrykit.surface.systems.catalysis`) is a small, self-contained
set of algebraic relationships built directly on
:func:`chemistrykit.kinetics.systems.arrhenius.arrhenius_rate_constant` --
neither shares polymorphic behavior worth capturing in an ABC, so
(following ``chemistrykit.electrochem``/``chemistrykit.photochem``'s
precedent) both stay as plain functions and small result dataclasses in
their own ``systems/`` modules.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

__all__ = ["AdsorptionIsotherm"]


class AdsorptionIsotherm(ABC):
    """Common interface for an equilibrium adsorption isotherm.

    Concrete subclasses (:class:`chemistrykit.surface.systems.langmuir.LangmuirIsotherm`,
    :class:`chemistrykit.surface.systems.freundlich.FreundlichIsotherm`,
    :class:`chemistrykit.surface.systems.bet.BETIsotherm`) implement
    :meth:`loading`, the amount adsorbed at a given equilibrium pressure;
    :meth:`fractional_coverage` (relative to a model-specific saturation
    loading) is then available for every subclass for free.
    """

    #: float: Saturation (monolayer) loading used to normalize
    #: :meth:`fractional_coverage`. Concrete subclasses set this in
    #: ``__init__`` (Langmuir/BET: a genuine monolayer capacity;
    #: Freundlich has no saturation loading, so its subclass overrides
    #: :meth:`fractional_coverage` to raise instead).
    saturation_loading: float = 1.0

    @abstractmethod
    def loading(self, P):
        """Return the amount adsorbed at equilibrium pressure(s) `P`.

        Parameters
        ----------
        P : float or array-like of float
            Equilibrium (partial) pressure of the adsorbate.

        Returns
        -------
        float or ndarray
            Amount adsorbed, in whatever units the model's capacity
            parameter was given (e.g. mol/g, cm^3(STP)/g).
        """

    def fractional_coverage(self, P):
        """Return the loading as a fraction of the saturation (monolayer) loading.

        Parameters
        ----------
        P : float or array-like of float

        Returns
        -------
        float or ndarray
            In :math:`[0, 1]` for a genuine monolayer model.
        """
        return np.asarray(self.loading(P), dtype=np.float64) / self.saturation_loading
