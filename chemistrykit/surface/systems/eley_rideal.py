r"""Eley-Rideal surface-reaction kinetics.

D. D. Eley and E. K. Rideal, "Parahydrogen Conversion on Tungsten,"
*Nature* 146 (1940), 401. In the Eley-Rideal mechanism only reactant A
is adsorbed (Langmuir coverage :math:`\theta_A`); B reacts with it
*directly from the gas phase*, without adsorbing first (Atkins & de
Paula, *Physical Chemistry*, 11th ed.):

.. math::

    \text{rate} = k\,\theta_A P_B = \frac{kK_AP_AP_B}{1+K_AP_A}

Unlike the dual-site Langmuir-Hinshelwood rate
(:func:`chemistrykit.surface.systems.langmuir_hinshelwood.lh_rate_dual_site`),
the rate therefore rises *monotonically* with :math:`P_A` (B never has to
compete with A for sites) and stays first order in :math:`P_B` at every
pressure -- the standard kinetic test for telling the two mechanisms
apart.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.surface.systems.langmuir import langmuir_coverage

__all__ = ["er_rate"]


def er_rate(k: float, K_A: float, P_A, P_B):
    r"""Eley-Rideal rate :math:`k\theta_AP_B`, with :math:`\theta_A` the Langmuir coverage of A.

    Parameters
    ----------
    k : float
        Rate constant for gas-phase B striking adsorbed A.
    K_A : float
        Langmuir adsorption equilibrium constant of A.
    P_A, P_B : float or array-like of float
        Partial pressures of the adsorbed reactant A and gas-phase reactant B.

    Returns
    -------
    float or ndarray

    Examples
    --------
    At the half-saturation pressure of A the rate is :math:`kP_B/2`:

    >>> round(float(er_rate(k=4.0, K_A=2.0, P_A=0.5, P_B=3.0)), 10)
    6.0

    Doubling :math:`P_B` always doubles the rate (first order in B):

    >>> r1 = er_rate(4.0, 2.0, 10.0, 1.0)
    >>> r2 = er_rate(4.0, 2.0, 10.0, 2.0)
    >>> round(float(r2 / r1), 10)
    2.0
    """
    theta_A = np.asarray(langmuir_coverage(K_A, P_A), dtype=np.float64)
    result = k * theta_A * np.asarray(P_B, dtype=np.float64)
    return float(result) if result.ndim == 0 else result
