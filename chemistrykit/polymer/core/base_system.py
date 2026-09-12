r"""Abstract base class for polymer chain-statistics models, and a note on why no other ABC is added here.

An ideal (random-walk, theta-solvent) chain and a real chain under good-
or poor-solvent conditions are described by the same two headline
quantities -- the mean-square end-to-end distance and the mean-square
radius of gyration -- differing only in how each scales with the number
of segments :math:`n` (Rubinstein & Colby, *Polymer Physics*, 2003, Ch.
2-3): ideally, :math:`\langle R^2\rangle = nb^2` *exactly* (a simple
random walk); more generally (real chains, any solvent quality),
:math:`R \sim bn^\nu` with a solvent-quality-dependent Flory exponent
:math:`\nu`. This is exactly the "several interchangeable models of one
physical relationship, compared side by side" shape that
:class:`chemistrykit.thermo.core.base_system.EquationOfState` already
captures for pure-substance PVT equations of state, and
:class:`chemistrykit.surface.core.base_system.AdsorptionIsotherm`
captures for adsorption models, so :class:`PolymerChainModel` is this
domain's one ABC.

Molecular-weight distributions (Mn, Mw, PDI, the Flory-Schulz
distribution), step-growth kinetics (the Carothers equation), and
chain-growth/free-radical polymerization kinetics are each a
self-contained closed-form relationship or reaction-network calculation
with no shared polymorphic interface among them, so -- following
``chemistrykit.electrochem``/``chemistrykit.photochem``'s precedent of
not forcing an ABC where the underlying models are algebraically
independent -- they stay as plain functions and small result dataclasses
in their own ``systems/`` modules rather than folded into an ABC here.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

__all__ = ["PolymerChainModel"]


class PolymerChainModel(ABC):
    """Common interface for a polymer chain's characteristic size as a function of chain length.

    Concrete subclasses (:class:`chemistrykit.polymer.systems.chain_statistics.IdealChain`,
    :class:`chemistrykit.polymer.systems.chain_statistics.RealChain`)
    implement :meth:`mean_square_end_to_end` and
    :meth:`mean_square_radius_of_gyration`; :meth:`end_to_end_distance`
    and :meth:`radius_of_gyration` (their square roots) are then
    available for every subclass for free.
    """

    @abstractmethod
    def mean_square_end_to_end(self, n, b):
        r"""Return :math:`\langle R^2\rangle` for a chain of `n` segments of length `b`.

        Parameters
        ----------
        n : float or array-like of float
            Number of chain segments (statistical/Kuhn segments).
        b : float
            Segment (Kuhn) length.

        Returns
        -------
        float or ndarray
        """

    @abstractmethod
    def mean_square_radius_of_gyration(self, n, b):
        r"""Return :math:`\langle R_g^2\rangle` for a chain of `n` segments of length `b`.

        Parameters
        ----------
        n : float or array-like of float
        b : float

        Returns
        -------
        float or ndarray
        """

    def end_to_end_distance(self, n, b):
        r"""Return :math:`R = \sqrt{\langle R^2\rangle}`.

        Parameters
        ----------
        n : float or array-like of float
        b : float

        Returns
        -------
        float or ndarray
        """
        return np.sqrt(self.mean_square_end_to_end(n, b))

    def radius_of_gyration(self, n, b):
        r"""Return :math:`R_g = \sqrt{\langle R_g^2\rangle}`.

        Parameters
        ----------
        n : float or array-like of float
        b : float

        Returns
        -------
        float or ndarray
        """
        return np.sqrt(self.mean_square_radius_of_gyration(n, b))
