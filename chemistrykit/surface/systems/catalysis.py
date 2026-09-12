r"""Turnover frequency/number, and rate enhancement from a catalyst's lowered activation energy.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 22.5-22.6. A
catalyst does not change a reaction's thermodynamics (:math:`\Delta G`,
:math:`K_{eq}`); it opens an alternative mechanistic pathway with a
*lower activation energy* :math:`E_a`, which -- via the Arrhenius
equation -- increases the rate constant by a factor of
:math:`\exp(\Delta E_a/RT)`, :math:`\Delta E_a` being the reduction in
activation energy. This module reuses
:func:`chemistrykit.kinetics.systems.arrhenius.arrhenius_rate_constant`
directly for that calculation rather than reimplementing the Arrhenius
equation, per this package's convention of not duplicating shared
kinetics machinery across domains.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from chemistrykit.constants import R
from chemistrykit.kinetics.systems.arrhenius import arrhenius_rate_constant

__all__ = ["turnover_frequency", "turnover_number", "CatalyticRateComparison", "compare_catalyzed_rate"]


def turnover_frequency(rate: float, active_site_concentration: float) -> float:
    r"""Turnover frequency :math:`\text{TOF} = \text{rate}/[\text{active sites}]`.

    The rate of product formation per catalytic active site per unit
    time -- the standard measure of intrinsic catalytic activity,
    independent of how much catalyst is used (Atkins & de Paula,
    *Physical Chemistry*, 11th ed., Ch. 22.6).

    Parameters
    ----------
    rate : float
        Reaction (product-formation) rate, in concentration/time.
    active_site_concentration : float
        Concentration (or amount) of catalytically active sites, in the
        same amount units as `rate`'s numerator.

    Returns
    -------
    float
        TOF, in 1/time.

    Examples
    --------
    >>> round(turnover_frequency(rate=5.0e-3, active_site_concentration=2.0e-6), 2)
    2500.0
    """
    return rate / active_site_concentration


def turnover_number(moles_converted: float, moles_catalyst: float) -> float:
    r"""Turnover number :math:`\text{TON} = n_{\text{converted}}/n_{\text{catalyst}}`.

    The total number of catalytic cycles a given amount of catalyst has
    performed over the course of a reaction -- a dimensionless measure of
    catalyst *durability/efficiency*, as distinct from :func:`turnover_frequency`'s
    measure of instantaneous *activity* (Atkins & de Paula, *Physical
    Chemistry*, 11th ed., Ch. 22.6).

    Parameters
    ----------
    moles_converted : float
        Total moles of substrate converted to product.
    moles_catalyst : float
        Moles of catalyst (active sites) used.

    Returns
    -------
    float
        Dimensionless.

    Examples
    --------
    >>> round(turnover_number(moles_converted=0.5, moles_catalyst=1.0e-4), 1)
    5000.0
    """
    return moles_converted / moles_catalyst


@dataclass
class CatalyticRateComparison:
    """Result of comparing a catalyzed and uncatalyzed rate constant at a given temperature."""

    k_uncatalyzed: float
    """float: Uncatalyzed rate constant at temperature `T`."""

    k_catalyzed: float
    """float: Catalyzed rate constant at temperature `T`."""

    rate_enhancement: float
    """float: :math:`k_{catalyzed}/k_{uncatalyzed}`."""

    delta_Ea: float
    """float: Activation-energy reduction :math:`E_{a,uncat}-E_{a,cat}`, in J/mol."""


def compare_catalyzed_rate(
    Ea_uncatalyzed: float,
    Ea_catalyzed: float,
    T: float,
    A_uncatalyzed: float,
    A_catalyzed: Optional[float] = None,
    R_gas: float = R,
) -> CatalyticRateComparison:
    r"""Compare catalyzed vs. uncatalyzed Arrhenius rate constants at temperature `T`.

    Both rate constants are evaluated with
    :func:`chemistrykit.kinetics.systems.arrhenius.arrhenius_rate_constant`;
    when the two mechanisms share the same pre-exponential factor `A`
    (the common simplifying assumption -- see Atkins & de Paula, *Physical
    Chemistry*, 11th ed., Ch. 22.5), the rate enhancement reduces to a
    pure exponential in the activation-energy reduction:

    .. math::

        \frac{k_{cat}}{k_{uncat}} = \exp\!\left(\frac{E_{a,uncat}-E_{a,cat}}{RT}\right)

    Parameters
    ----------
    Ea_uncatalyzed, Ea_catalyzed : float
        Activation energies of the uncatalyzed and catalyzed pathways, in
        J/mol (``Ea_catalyzed < Ea_uncatalyzed`` for a genuine catalyst).
    T : float
        Absolute temperature, in K.
    A_uncatalyzed : float
        Pre-exponential factor of the uncatalyzed pathway.
    A_catalyzed : float, optional
        Pre-exponential factor of the catalyzed pathway. Defaults to
        `A_uncatalyzed` (the common simplifying assumption that the
        catalyst changes only :math:`E_a`, not the attempt frequency).
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.

    Returns
    -------
    CatalyticRateComparison

    Examples
    --------
    With equal pre-exponential factors, the rate enhancement is exactly
    :math:`\exp(\Delta E_a/RT)`:

    >>> import numpy as np
    >>> comparison = compare_catalyzed_rate(Ea_uncatalyzed=80e3, Ea_catalyzed=50e3, T=298.15, A_uncatalyzed=1e13)
    >>> expected = np.exp((80e3 - 50e3) / (8.31446261815324 * 298.15))
    >>> round(float(comparison.rate_enhancement / expected), 6)
    1.0

    A catalyst that lowers :math:`E_a` by 30 kJ/mol at room temperature
    speeds the reaction up by many orders of magnitude:

    >>> comparison.rate_enhancement > 1.0e5
    True
    """
    if A_catalyzed is None:
        A_catalyzed = A_uncatalyzed
    k_un = float(arrhenius_rate_constant(A_uncatalyzed, Ea_uncatalyzed, T, R_gas))
    k_cat = float(arrhenius_rate_constant(A_catalyzed, Ea_catalyzed, T, R_gas))
    return CatalyticRateComparison(
        k_uncatalyzed=k_un,
        k_catalyzed=k_cat,
        rate_enhancement=k_cat / k_un,
        delta_Ea=Ea_uncatalyzed - Ea_catalyzed,
    )
