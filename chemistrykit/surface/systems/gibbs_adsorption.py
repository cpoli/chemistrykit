r"""The Gibbs adsorption equation, and the Szyszkowski surface-tension equation.

J. W. Gibbs, "On the Equilibrium of Heterogeneous Substances," *Trans.
Connecticut Acad. Arts Sci.* 3 (1875-1878), whose theory of surfaces of
discontinuity gives the surface excess :math:`\Gamma` of a dilute solute
directly from how the liquid's surface tension :math:`\gamma` changes
with the solute's concentration :math:`c`:

.. math::

    \Gamma = -\frac{1}{RT}\,\frac{d\gamma}{d\ln c}

(Atkins & de Paula, *Physical Chemistry*, 11th ed.; Adamson &
Gast, *Physical Chemistry of Surfaces*, 6th ed., Ch. III). A solute that
*lowers* the surface tension (a surfactant) is therefore positively
adsorbed at the surface -- measured without ever looking at the surface
directly.

The empirical Szyszkowski (1908) equation for aqueous fatty-acid
solutions,

.. math::

    \gamma = \gamma_0 - RT\,\Gamma_{max}\ln(1 + Kc),

turns, under the Gibbs equation, into exactly a Langmuir isotherm for the
surface excess, :math:`\Gamma = \Gamma_{max}Kc/(1+Kc)` -- the link that
:func:`gibbs_surface_excess` recovers numerically from surface-tension
data alone.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import R

__all__ = ["szyszkowski_surface_tension", "gibbs_surface_excess"]


def szyszkowski_surface_tension(c, gamma0: float, Gamma_max: float, K: float, T: float, R_gas: float = R):
    r"""Szyszkowski surface tension :math:`\gamma = \gamma_0 - RT\Gamma_{max}\ln(1+Kc)`.

    Parameters
    ----------
    c : float or array-like of float
        Bulk solute concentration (mol/m^3 or any unit consistent with `K`).
    gamma0 : float
        Surface tension of the pure solvent, in N/m.
    Gamma_max : float
        Saturation surface excess, in mol/m^2.
    K : float
        Adsorption constant, inverse to `c`'s units.
    T : float
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.

    Returns
    -------
    float or ndarray
        Surface tension, in N/m.

    Examples
    --------
    Pure solvent (``c=0``) has surface tension exactly :math:`\gamma_0`:

    >>> round(float(szyszkowski_surface_tension(0.0, gamma0=0.072, Gamma_max=5e-6, K=1.0, T=298.15)), 10)
    0.072
    """
    c = np.asarray(c, dtype=np.float64)
    result = gamma0 - R_gas * T * Gamma_max * np.log1p(K * c)
    return float(result) if result.ndim == 0 else result


def gibbs_surface_excess(c, gamma, T: float, R_gas: float = R):
    r"""Surface excess from surface-tension data via the Gibbs equation :math:`\Gamma=-(1/RT)\,d\gamma/d\ln c`.

    The derivative is taken numerically (second-order finite differences
    in :math:`\ln c`, :func:`numpy.gradient`), so `c` should be sampled
    densely enough that :math:`\gamma(\ln c)` is smooth.

    Parameters
    ----------
    c : array-like of float
        Strictly positive, increasing bulk concentrations.
    gamma : array-like of float
        Measured surface tensions at those concentrations, in N/m.
    T : float
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.

    Returns
    -------
    ndarray
        Surface excess :math:`\Gamma`, in mol/m^2 (for `gamma` in N/m).

    Examples
    --------
    Applied to Szyszkowski surface-tension data, the Gibbs equation
    recovers a Langmuir isotherm for the surface excess:

    >>> c = np.logspace(-3, 2, 2001)
    >>> gamma = szyszkowski_surface_tension(c, gamma0=0.072, Gamma_max=5e-6, K=2.0, T=298.15)
    >>> Gamma = gibbs_surface_excess(c, gamma, T=298.15)
    >>> bool(np.allclose(Gamma[1:-1], 5e-6 * 2.0 * c[1:-1] / (1 + 2.0 * c[1:-1]), rtol=1e-4))
    True
    """
    c = np.asarray(c, dtype=np.float64)
    gamma = np.asarray(gamma, dtype=np.float64)
    if np.any(c <= 0):
        raise ValueError("concentrations must be strictly positive")
    dgamma_dlnc = np.gradient(gamma, np.log(c))
    return -dgamma_dlnc / (R_gas * T)
