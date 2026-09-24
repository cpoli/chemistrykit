r"""The Debye model of a crystalline solid's heat capacity.

P. Debye, "Zur Theorie der spezifischen Wärmen," *Ann. Phys.* 39,
789-839 (1912), replaced Einstein's single vibrational frequency with a
continuous spectrum of lattice vibrations (sound waves), with density of
states :math:`g(\nu)\propto\nu^2` cut off at a maximum frequency
:math:`\nu_D` chosen so the solid has exactly :math:`3N` modes. Writing
:math:`\Theta_D=h\nu_D/k_B` (the Debye temperature):

.. math::

    C_V = 9Nk_B\left(\frac{T}{\Theta_D}\right)^3
          \int_0^{\Theta_D/T}\frac{x^4e^x}{(e^x-1)^2}\,dx

which tends to the Dulong-Petit value :math:`3Nk_B` at high `T` and to
Debye's :math:`T^3` law, :math:`C_V\to\frac{12\pi^4}{5}Nk_B(T/\Theta_D)^3`,
as :math:`T\to0` (McQuarrie, *Statistical Mechanics*, Ch. 11).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from chemistrykit.constants import K_B, NA

__all__ = ["DebyeSolid"]


def _debye_integral(power: int, upper: float) -> float:
    r"""Return :math:`\int_0^{u}x^{p}e^x/(e^x-1)^2\,dx` (``power=4``) or :math:`\int_0^u x^3/(e^x-1)\,dx` (``power=3``)."""
    if upper <= 0.0:
        return 0.0
    if power == 4:
        integrand = lambda x: x**4 * np.exp(-x) / (-np.expm1(-x)) ** 2 if x > 0 else 0.0
    else:
        integrand = lambda x: x**3 / np.expm1(x) if x > 0 else 0.0
    value, _ = quad(integrand, 0.0, min(upper, 200.0), limit=200)
    return value


class DebyeSolid:
    r"""Debye's continuum model of a monatomic crystal's lattice vibrations.

    Parameters
    ----------
    debye_temperature : float
        The Debye temperature :math:`\Theta_D=h\nu_D/k_B`, in K.

    Examples
    --------
    Far above :math:`\Theta_D` the molar heat capacity approaches the
    Dulong-Petit value :math:`3R`; far below, it follows Debye's
    :math:`T^3` law:

    >>> copper = DebyeSolid(debye_temperature=343.0)
    >>> round(float(copper.heat_capacity_v(5000.0) / (3 * 8.31446261815324)), 3)
    1.0
    >>> ratio = copper.heat_capacity_v(3.0) / copper.low_temperature_heat_capacity(3.0)
    >>> round(float(ratio), 6)
    1.0
    """

    def __init__(self, debye_temperature: float):
        if debye_temperature <= 0:
            raise ValueError("debye_temperature must be positive")
        self.debye_temperature = float(debye_temperature)

    def heat_capacity_v(self, T, N: float = NA):
        """Return the constant-volume heat capacity of `N` atoms (default one mole).

        Parameters
        ----------
        T : float or array_like
            Absolute temperature, in K (must be positive).
        N : float, default :data:`chemistrykit.constants.NA`
            Number of atoms.

        Returns
        -------
        float or numpy.ndarray
            Heat capacity, in J/K.
        """
        T_arr = np.asarray(T, dtype=float)
        if np.any(T_arr <= 0):
            raise ValueError("T must be positive")
        values = np.vectorize(lambda t: 9.0 * N * K_B * (t / self.debye_temperature) ** 3 * _debye_integral(4, self.debye_temperature / t))(T_arr)
        return float(values) if values.ndim == 0 else values

    def internal_energy(self, T, N: float = NA):
        r"""Return the thermal internal energy :math:`U(T)-U(0)` of `N` atoms.

        :math:`U-U_0 = 9Nk_BT(T/\Theta_D)^3\int_0^{\Theta_D/T}x^3/(e^x-1)\,dx`,
        measured from the zero-point energy :math:`U_0=\frac98Nk_B\Theta_D`.

        Parameters
        ----------
        T : float or array_like
            Absolute temperature, in K (must be positive).
        N : float, default :data:`chemistrykit.constants.NA`

        Returns
        -------
        float or numpy.ndarray
            Energy, in J.
        """
        T_arr = np.asarray(T, dtype=float)
        if np.any(T_arr <= 0):
            raise ValueError("T must be positive")
        values = np.vectorize(lambda t: 9.0 * N * K_B * t * (t / self.debye_temperature) ** 3 * _debye_integral(3, self.debye_temperature / t))(T_arr)
        return float(values) if values.ndim == 0 else values

    def low_temperature_heat_capacity(self, T, N: float = NA):
        r"""Return Debye's :math:`T^3` law, :math:`C_V=\frac{12\pi^4}{5}Nk_B(T/\Theta_D)^3`, exact as :math:`T\to0`.

        Parameters
        ----------
        T : float or array_like
            Absolute temperature, in K.
        N : float, default :data:`chemistrykit.constants.NA`

        Returns
        -------
        float or numpy.ndarray
            Heat capacity, in J/K.
        """
        return 12.0 * np.pi**4 / 5.0 * N * K_B * (np.asarray(T, dtype=float) / self.debye_temperature) ** 3
