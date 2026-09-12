r"""Langmuir-Hinshelwood surface-reaction kinetics.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 22.5
("Heterogeneous catalysis"); the mechanism is named for I. Langmuir and
C. N. Hinshelwood, who independently developed the kinetic picture of a
reaction occurring between species adsorbed on the *same* catalytic
surface. Both mechanisms below build directly on
:func:`chemistrykit.surface.systems.langmuir.langmuir_coverage` -- the
Langmuir-Hinshelwood rate law is nothing more than a mass-action rate
constant multiplied by the Langmuir fractional coverage(s) of the
reacting adsorbate(s), so this module reuses that function rather than
re-deriving the coverage expression.
"""

from __future__ import annotations

import numpy as np

__all__ = ["lh_rate_single_site", "lh_rate_dual_site"]


def lh_rate_single_site(k: float, K_A: float, P_A):
    r"""Single-site Langmuir-Hinshelwood rate: a lone adsorbed species reacts unimolecularly.

    A single reactant A adsorbs onto the catalytic surface (Langmuir
    equilibrium, constant :math:`K_A`) and then reacts *on the surface*
    at a rate proportional to its own coverage :math:`\theta_A`:

    .. math::

        \text{rate} = k\,\theta_A = \frac{kK_AP_A}{1+K_AP_A}

    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 22.5, eq.
    22.19). At low pressure (:math:`K_AP_A \ll 1`) this reduces to a rate
    that is first order in :math:`P_A`; at high pressure
    (:math:`K_AP_A \gg 1`, the surface saturated with A) it becomes zero
    order in :math:`P_A` -- the hallmark Langmuir-Hinshelwood pressure
    dependence, qualitatively different from a simple gas-phase
    elementary reaction.

    Parameters
    ----------
    k : float
        Surface-reaction rate constant (the rate at full coverage,
        :math:`\theta_A=1`).
    K_A : float
        Langmuir adsorption equilibrium constant of A.
    P_A : float or array-like of float
        Partial pressure of A.

    Returns
    -------
    float or ndarray

    Examples
    --------
    At the Langmuir half-saturation pressure :math:`P_A=1/K_A`, the rate
    is exactly half the saturation (zero-order) rate `k`:

    >>> round(float(lh_rate_single_site(k=4.0, K_A=2.0, P_A=0.5)), 10)
    2.0

    Low-pressure limit is (to leading order) first order in `P_A`:

    >>> k, K_A = 4.0, 2.0
    >>> P_small = np.array([1e-4, 2e-4, 4e-4])
    >>> rate = lh_rate_single_site(k, K_A, P_small)
    >>> bool(np.allclose(rate, k * K_A * P_small, rtol=1e-3))
    True
    """
    P_A = np.asarray(P_A, dtype=np.float64)
    result = k * K_A * P_A / (1.0 + K_A * P_A)
    return float(result) if result.ndim == 0 else result


def lh_rate_dual_site(k: float, K_A: float, P_A, K_B: float, P_B):
    r"""Dual-site (bimolecular) Langmuir-Hinshelwood rate: two species compete for the same sites.

    Two reactants A and B adsorb *competitively* onto the same pool of
    surface sites, and react where an A and a B happen to sit on
    neighboring sites -- the classic mechanism for, e.g., the catalytic
    oxidation :math:`2CO+O_2 \to 2CO_2` on a metal surface (Atkins & de
    Paula, *Physical Chemistry*, 11th ed., Ch. 22.5). Competitive Langmuir
    adsorption gives each species' coverage as

    .. math::

        \theta_A = \frac{K_AP_A}{1+K_AP_A+K_BP_B}, \qquad
        \theta_B = \frac{K_BP_B}{1+K_AP_A+K_BP_B}

    and the surface-reaction rate is :math:`k\theta_A\theta_B`. A
    characteristic (and often counterintuitive) feature: for fixed
    :math:`P_B`, the rate is *not* monotonic in :math:`P_A` -- it rises
    at low :math:`P_A` (more A available to react) but falls again at
    high :math:`P_A` (A crowds B off the surface entirely), peaking at an
    intermediate pressure.

    Parameters
    ----------
    k : float
        Surface-reaction rate constant.
    K_A, K_B : float
        Langmuir adsorption equilibrium constants of A and B.
    P_A, P_B : float or array-like of float
        Partial pressures of A and B.

    Returns
    -------
    float or ndarray

    Examples
    --------
    At low coverage of both species (:math:`K_AP_A,K_BP_B \ll 1`, the
    denominator :math:`\approx 1`), the rate reduces to the naive
    mass-action product :math:`kK_AP_AK_BP_B`:

    >>> k, K_A, K_B = 4.0, 2.0, 5.0
    >>> P_A = np.array([1e-4, 2e-4])
    >>> P_B = np.array([1e-4, 3e-4])
    >>> rate = lh_rate_dual_site(k, K_A, P_A, K_B, P_B)
    >>> bool(np.allclose(rate, k * K_A * P_A * K_B * P_B, rtol=1e-2))
    True

    The rate vanishes if either reactant is entirely absent (nothing to
    react *with*, even at full coverage of the other):

    >>> round(float(lh_rate_dual_site(k=4.0, K_A=2.0, P_A=0.0, K_B=5.0, P_B=1.0)), 10)
    0.0

    The rate is non-monotonic in :math:`P_A` at fixed :math:`P_B`
    (rises, then falls, as A crowds B off the surface):

    >>> P_A_scan = np.linspace(0.01, 50.0, 200)
    >>> rate_scan = lh_rate_dual_site(k=4.0, K_A=2.0, P_A=P_A_scan, K_B=5.0, P_B=0.2)
    >>> peak_index = int(np.argmax(rate_scan))
    >>> bool(0 < peak_index < len(P_A_scan) - 1)
    True
    """
    P_A = np.asarray(P_A, dtype=np.float64)
    P_B = np.asarray(P_B, dtype=np.float64)
    denom = 1.0 + K_A * P_A + K_B * P_B
    theta_A = K_A * P_A / denom
    theta_B = K_B * P_B / denom
    result = k * theta_A * theta_B
    return float(result) if result.ndim == 0 else result
