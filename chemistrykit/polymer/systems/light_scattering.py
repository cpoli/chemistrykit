r"""Absolute molar mass from dilute-solution light scattering and osmometry.

P. Debye, *J. Appl. Phys.* 15, 338 (1944). For a dilute solution of
polymer species :math:`i` (mass concentration :math:`c_i`, molar mass
:math:`M_i`), each chain scatters in proportion to the square of its
mass, so the zero-angle excess Rayleigh ratio is
:math:`R_0=K\sum_ic_iM_i` (optical constant :math:`K`), and the Debye
relation

.. math::

    \frac{Kc}{R_0} = \frac{1}{M_w} + 2A_2c + \dots

gives, extrapolated to :math:`c\to0`, the *weight*-average molar mass
:math:`M_w=\sum c_iM_i/\sum c_i`. Osmotic pressure instead counts
molecules, :math:`\Pi=RT\sum c_i/M_i`, so :math:`\Pi/(RTc)\to1/M_n`.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "rayleigh_ratio_dilute_mixture",
    "debye_Kc_over_R",
    "osmotic_pressure_dilute_mixture",
]


def rayleigh_ratio_dilute_mixture(c_i, M_i, K: float = 1.0) -> float:
    r"""Zero-angle excess Rayleigh ratio of an ideal dilute mixture, :math:`R_0=K\sum c_iM_i`.

    Parameters
    ----------
    c_i : array-like of float
        Mass concentration of each species.
    M_i : array-like of float
        Molar mass of each species.
    K : float, optional
        Optical constant.

    Returns
    -------
    float

    Examples
    --------
    The apparent molar mass :math:`R_0/(Kc)` is the weight average:

    >>> round(rayleigh_ratio_dilute_mixture([1.0, 1.0], [1e4, 1e5]) / 2.0, 6)
    55000.0
    """
    return float(K * np.sum(np.asarray(c_i, dtype=float) * np.asarray(M_i, dtype=float)))


def debye_Kc_over_R(c, Mw: float, A2: float = 0.0):
    r"""Debye relation :math:`Kc/R_0=1/M_w+2A_2c` (zero scattering angle).

    Parameters
    ----------
    c : float or array-like of float
        Total polymer mass concentration.
    Mw : float
        Weight-average molar mass.
    A2 : float, optional
        Second virial coefficient.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> float(debye_Kc_over_R(0.0, Mw=1e5, A2=1e-4))
    1e-05
    """
    return 1.0 / Mw + 2.0 * A2 * np.asarray(c, dtype=float)


def osmotic_pressure_dilute_mixture(c_i, M_i, T: float, R: float = 8.314462618) -> float:
    r"""Van 't Hoff osmotic pressure of an ideal dilute mixture, :math:`\Pi=RT\sum c_i/M_i`.

    Parameters
    ----------
    c_i : array-like of float
        Mass concentrations (kg m^-3).
    M_i : array-like of float
        Molar masses (kg mol^-1).
    T : float
        Temperature (K).
    R : float, optional
        Gas constant.

    Returns
    -------
    float
        Osmotic pressure (Pa).

    Examples
    --------
    The apparent molar mass :math:`RTc/\Pi` is the number average:

    >>> Pi = osmotic_pressure_dilute_mixture([1.0, 1.0], [10.0, 100.0], T=300.0)
    >>> round(8.314462618 * 300.0 * 2.0 / Pi, 6)
    18.181818
    """
    return float(R * T * np.sum(np.asarray(c_i, dtype=float) / np.asarray(M_i, dtype=float)))
