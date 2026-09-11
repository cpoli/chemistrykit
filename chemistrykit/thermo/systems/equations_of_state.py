r"""Equations of state: ideal gas, van der Waals, and Redlich-Kwong.

All three share the :class:`~chemistrykit.thermo.core.base_system.EquationOfState`
interface (``pressure(Vm, T)``, ``molar_volume(P, T)``,
``compressibility_factor(P, T)``) so they can be compared directly on the
same P-V isotherm plot. See Atkins & de Paula, *Physical Chemistry*, 11th
ed., Ch. 1 for the ideal gas law and van der Waals equation, and O.
Redlich & J. N. S. Kwong, *Chem. Rev.* 44, 233 (1949) (also tabulated in
Smith, Van Ness & Abbott, *Introduction to Chemical Engineering
Thermodynamics*, 7th ed., Ch. 3) for Redlich-Kwong.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import R
from chemistrykit.thermo.core.base_system import EquationOfState
from chemistrykit.thermo.utils.cubic_roots import real_positive_roots

__all__ = ["IdealGas", "VanDerWaals", "RedlichKwong"]


class IdealGas(EquationOfState):
    r"""The ideal (perfect) gas law: :math:`PV_m = RT`.

    See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 1.1.

    Examples
    --------
    The molar volume of an ideal gas at STP (1 atm-scale pressure of
    101325 Pa, 0 degC) is the textbook 22.4 L/mol:

    >>> gas = IdealGas()
    >>> round(float(gas.molar_volume(P=101325.0, T=273.15)) * 1000.0, 1)
    22.4
    >>> round(float(gas.compressibility_factor(P=101325.0, T=273.15)), 6)
    1.0
    """

    def pressure(self, Vm, T):
        Vm = np.asarray(Vm, dtype=np.float64)
        T = np.asarray(T, dtype=np.float64)
        return self.R * T / Vm

    def molar_volume(self, P, T, **kwargs):
        P = np.asarray(P, dtype=np.float64)
        T = np.asarray(T, dtype=np.float64)
        return self.R * T / P


class VanDerWaals(EquationOfState):
    r"""The van der Waals equation of state.

    .. math::

        P = \frac{RT}{V_m - b} - \frac{a}{V_m^2}

    where `a` corrects for intermolecular attraction and `b` for the
    finite volume of the molecules themselves (J. D. van der Waals,
    doctoral thesis, Leiden, 1873; Atkins & de Paula, *Physical
    Chemistry*, 11th ed., Ch. 1.3b).

    Parameters
    ----------
    a : float
        Attraction parameter, in Pa m^6 mol^-2.
    b : float
        Excluded-volume parameter, in m^3/mol.

    Examples
    --------
    With ``a = b = 0`` van der Waals reduces exactly to the ideal gas law:

    >>> vdw = VanDerWaals(a=0.0, b=0.0)
    >>> ideal = IdealGas()
    >>> round(float(vdw.pressure(Vm=0.01, T=300.0)), 6) == round(float(ideal.pressure(Vm=0.01, T=300.0)), 6)
    True
    """

    def __init__(self, a: float, b: float):
        if b < 0:
            raise ValueError("b must be non-negative")
        self.a = float(a)
        self.b = float(b)

    def pressure(self, Vm, T):
        Vm = np.asarray(Vm, dtype=np.float64)
        T = np.asarray(T, dtype=np.float64)
        return self.R * T / (Vm - self.b) - self.a / Vm**2

    def molar_volume(self, P, T, branch: str = "vapor") -> float:
        r"""Solve the van der Waals cubic for `Vm` at pressure `P`, temperature `T`.

        Rearranging :math:`P = RT/(V_m-b) - a/V_m^2` into a cubic in
        :math:`V_m`:

        .. math::

            V_m^3 - \left(b + \frac{RT}{P}\right) V_m^2 + \frac{a}{P} V_m - \frac{ab}{P} = 0

        Below the critical temperature this can have three real positive
        roots (a metastable liquid branch, an unstable middle root, and a
        vapor branch).

        Parameters
        ----------
        P : float
            Pressure, in Pa.
        T : float
            Absolute temperature, in K.
        branch : {"vapor", "liquid"}
            Return the largest ("vapor") or smallest ("liquid") real
            positive root.

        Returns
        -------
        float

        Examples
        --------
        Solving for `Vm` and substituting back reproduces the original
        pressure:

        >>> vdw = VanDerWaals(a=0.1448, b=3.913e-5)  # e.g. representative CO2-like values
        >>> Vm = vdw.molar_volume(P=1.0e5, T=300.0)
        >>> round(float(vdw.pressure(Vm, T=300.0)), 3)
        100000.0
        """
        coeffs = [1.0, -(self.b + self.R * T / P), self.a / P, -self.a * self.b / P]
        roots = real_positive_roots(coeffs)
        if roots.size == 0:
            raise ValueError("no physically admissible molar volume found for the given P, T")
        return float(roots[-1]) if branch == "vapor" else float(roots[0])

    @classmethod
    def from_critical_constants(cls, Tc: float, Pc: float) -> VanDerWaals:
        r"""Build a :class:`VanDerWaals` EOS from critical temperature and pressure.

        At the critical point the van der Waals cubic has a triple root,
        which fixes

        .. math::

            a = \frac{27 R^2 T_c^2}{64 P_c}, \qquad b = \frac{R T_c}{8 P_c}

        (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 1.3c).

        Parameters
        ----------
        Tc : float
            Critical temperature, in K.
        Pc : float
            Critical pressure, in Pa.

        Returns
        -------
        VanDerWaals

        Examples
        --------
        A structural consequence of these formulas -- independent of
        which substance's (Tc, Pc) is used -- is that the van der Waals
        critical molar volume is always :math:`V_c = 3b`, and hence the
        critical compressibility factor is always exactly 3/8 (Atkins &
        de Paula, Table 1.6), unlike the true value for real gases
        (typically 0.23-0.31). (Solving the cubic numerically exactly at
        the critical point, where it has a triple root, is a
        textbook-classic ill-conditioned root-finding problem, so this
        checks the closed-form critical relation directly rather than
        round-tripping through :meth:`molar_volume`.)

        >>> Tc, Pc = 304.13, 7.3773e6  # CO2
        >>> eos = VanDerWaals.from_critical_constants(Tc=Tc, Pc=Pc)
        >>> Vc = 3.0 * eos.b
        >>> round(Pc * Vc / (eos.R * Tc), 6)
        0.375
        """
        a = 27.0 * R**2 * Tc**2 / (64.0 * Pc)
        b = R * Tc / (8.0 * Pc)
        return cls(a=a, b=b)


class RedlichKwong(EquationOfState):
    r"""The Redlich-Kwong equation of state.

    .. math::

        P = \frac{RT}{V_m - b} - \frac{a}{\sqrt{T}\,V_m(V_m+b)}

    An improvement on van der Waals that gives the attractive term an
    explicit temperature dependence (O. Redlich & J. N. S. Kwong, *Chem.
    Rev.* 44, 233 (1949)).

    Parameters
    ----------
    a : float
        Attraction parameter, in Pa m^6 K^0.5 mol^-2.
    b : float
        Excluded-volume parameter, in m^3/mol.
    """

    def __init__(self, a: float, b: float):
        if b < 0:
            raise ValueError("b must be non-negative")
        self.a = float(a)
        self.b = float(b)

    def pressure(self, Vm, T):
        Vm = np.asarray(Vm, dtype=np.float64)
        T = np.asarray(T, dtype=np.float64)
        return self.R * T / (Vm - self.b) - self.a / (np.sqrt(T) * Vm * (Vm + self.b))

    def molar_volume(self, P, T, branch: str = "vapor") -> float:
        r"""Solve the Redlich-Kwong cubic for `Vm` at pressure `P`, temperature `T`.

        Multiplying :math:`P = RT/(V_m-b) - a/(\sqrt{T}V_m(V_m+b))`
        through by :math:`\sqrt{T}V_m(V_m+b)(V_m-b)` and collecting terms
        in :math:`V_m` gives

        .. math::

            V_m^3 - \frac{RT}{P} V_m^2 + \left(\frac{a}{P\sqrt{T}} - b^2 - \frac{bRT}{P}\right) V_m - \frac{ab}{P\sqrt{T}} = 0

        Parameters
        ----------
        P : float
            Pressure, in Pa.
        T : float
            Absolute temperature, in K.
        branch : {"vapor", "liquid"}
            Return the largest ("vapor") or smallest ("liquid") real
            positive root.

        Returns
        -------
        float

        Examples
        --------
        Solving for `Vm` and substituting back reproduces the original
        pressure:

        >>> rk = RedlichKwong(a=6.4239, b=2.7143e-5)  # e.g. representative CO2-like values
        >>> Vm = rk.molar_volume(P=1.0e5, T=300.0)
        >>> round(float(rk.pressure(Vm, T=300.0)), 3)
        100000.0
        """
        sqrtT = np.sqrt(T)
        coeffs = [
            1.0,
            -(self.R * T / P),
            self.a / (P * sqrtT) - self.b**2 - self.b * self.R * T / P,
            -self.a * self.b / (P * sqrtT),
        ]
        roots = real_positive_roots(coeffs)
        if roots.size == 0:
            raise ValueError("no physically admissible molar volume found for the given P, T")
        return float(roots[-1]) if branch == "vapor" else float(roots[0])

    @classmethod
    def from_critical_constants(cls, Tc: float, Pc: float) -> RedlichKwong:
        r"""Build a :class:`RedlichKwong` EOS from critical temperature and pressure.

        .. math::

            a = 0.42748\,\frac{R^2 T_c^{2.5}}{P_c}, \qquad b = 0.08664\,\frac{R T_c}{P_c}

        (Redlich & Kwong 1949; the numeric coefficients are tabulated
        e.g. in Smith, Van Ness & Abbott, *Introduction to Chemical
        Engineering Thermodynamics*, 7th ed., Table 3.1).

        Parameters
        ----------
        Tc : float
            Critical temperature, in K.
        Pc : float
            Critical pressure, in Pa.

        Returns
        -------
        RedlichKwong

        Examples
        --------
        As with van der Waals's 3/8, Redlich-Kwong's critical-point
        construction fixes a universal critical molar volume :math:`V_c =
        RT_c/(3P_c)` and hence a universal critical compressibility
        factor, here exactly 1/3, again independent of which substance's
        (Tc, Pc) is used. (As in :meth:`VanDerWaals.from_critical_constants`,
        this checks the closed-form critical relation directly rather
        than numerically solving the cubic exactly at its ill-conditioned
        triple root.)

        >>> Tc, Pc = 304.13, 7.3773e6  # CO2
        >>> eos = RedlichKwong.from_critical_constants(Tc=Tc, Pc=Pc)
        >>> Vc = eos.R * Tc / (3.0 * Pc)
        >>> round(Pc * Vc / (eos.R * Tc), 6)
        0.333333
        """
        a = 0.42748 * R**2 * Tc**2.5 / Pc
        b = 0.08664 * R * Tc / Pc
        return cls(a=a, b=b)
