r"""Lewis fugacity: fugacity coefficients and vapor-liquid saturation from an equation of state.

G. N. Lewis defined the fugacity :math:`f` of a real gas so that its
chemical potential keeps the ideal-gas form,
:math:`\mu = \mu^\circ + RT\ln(f/P^\circ)`, with :math:`f \to P` as
:math:`P \to 0` (G. N. Lewis, *Proc. Am. Acad. Arts Sci.* 37, 49-69
(1901)). The fugacity coefficient :math:`\phi = f/P` follows from any
pressure-explicit equation of state as (Smith, Van Ness & Abbott,
*Introduction to Chemical Engineering Thermodynamics*, 7th ed., Ch. 11)

.. math::

    \ln\phi = Z - 1 - \ln Z + \frac{1}{RT}\int_{V_m}^{\infty}
        \left(P - \frac{RT}{V'}\right) dV'

and two phases of a pure substance coexist exactly when their fugacities
are equal, which fixes the saturation (vapor) pressure.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq, minimize_scalar

__all__ = ["fugacity_coefficient", "fugacity", "saturation_pressure"]

_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(96)


def fugacity_coefficient(eos, P: float, T: float, branch: str = "vapor") -> float:
    r"""Fugacity coefficient :math:`\phi = f/P` of a pure fluid described by `eos`.

    The residual integral is evaluated in density,
    :math:`\rho = 1/V_m`, as :math:`\int_0^{\rho}(P - RT\rho')/\rho'^2\,d\rho'`
    with 96-point Gauss-Legendre quadrature, so any
    :class:`~chemistrykit.thermo.core.base_system.EquationOfState` works.

    Parameters
    ----------
    eos : EquationOfState
        Pure-fluid equation of state.
    P : float
        Pressure, in Pa.
    T : float
        Absolute temperature, in K.
    branch : {"vapor", "liquid"}, default "vapor"
        Root of a cubic equation of state to use (ignored by
        :class:`~chemistrykit.thermo.systems.equations_of_state.IdealGas`).

    Returns
    -------
    float

    Examples
    --------
    For van der Waals the integral has the closed form
    :math:`\ln\phi = b/(V_m-b) - 2a/(RTV_m) - \ln[P(V_m-b)/(RT)]`:

    >>> import numpy as np
    >>> from chemistrykit.thermo import VanDerWaals
    >>> vdw = VanDerWaals(a=0.3640, b=4.267e-5)
    >>> P, T = 2.0e6, 300.0
    >>> V = vdw.molar_volume(P, T)
    >>> RT = vdw.R * T
    >>> exact = np.exp(vdw.b / (V - vdw.b) - 2 * vdw.a / (RT * V) - np.log(P * (V - vdw.b) / RT))
    >>> bool(np.isclose(fugacity_coefficient(vdw, P, T), exact, rtol=1e-10))
    True
    """
    try:
        Vm = eos.molar_volume(P, T, branch=branch)
    except TypeError:
        Vm = eos.molar_volume(P, T)
    Vm = float(Vm)
    RT = eos.R * T
    Z = P * Vm / RT
    rho_max = 1.0 / Vm
    rho = 0.5 * rho_max * (_GL_NODES + 1.0)
    integrand = (eos.pressure(1.0 / rho, T) - RT * rho) / rho**2
    residual = 0.5 * rho_max * np.sum(_GL_WEIGHTS * integrand) / RT
    return float(np.exp(Z - 1.0 - np.log(Z) + residual))


def fugacity(eos, P: float, T: float, branch: str = "vapor") -> float:
    r"""Fugacity :math:`f = \phi P` of a pure fluid, in Pa.

    Parameters
    ----------
    eos : EquationOfState
        Pure-fluid equation of state.
    P : float
        Pressure, in Pa.
    T : float
        Absolute temperature, in K.
    branch : {"vapor", "liquid"}, default "vapor"
        Root of a cubic equation of state to use.

    Returns
    -------
    float

    Examples
    --------
    An ideal gas's fugacity is its pressure:

    >>> from chemistrykit.thermo import IdealGas
    >>> round(fugacity(IdealGas(), 5.0e5, 300.0), 6)
    500000.0
    """
    return fugacity_coefficient(eos, P, T, branch=branch) * P


def saturation_pressure(eos, T: float) -> float:
    r"""Vapor pressure of a cubic equation of state from Lewis's equal-fugacity condition.

    Finds the pressure at which the liquid and vapor roots of `eos` have
    equal fugacity, :math:`f^L(P_{sat}) = f^V(P_{sat})`. The search is
    bracketed between the isotherm's two spinodal pressures (its local
    minimum and maximum in :math:`P(V_m)`), where three roots exist.

    Parameters
    ----------
    eos : EquationOfState
        A cubic equation of state with an excluded-volume attribute ``b``
        (:class:`~chemistrykit.thermo.VanDerWaals`,
        :class:`~chemistrykit.thermo.RedlichKwong`,
        :class:`~chemistrykit.thermo.PengRobinson`).
    T : float
        Absolute temperature, below the equation's critical temperature, in K.

    Returns
    -------
    float
        Saturation pressure, in Pa.

    Examples
    --------
    The van der Waals coexistence curve at reduced temperature 0.9 lies at
    reduced pressure 0.647, whatever the substance:

    >>> from chemistrykit.thermo import VanDerWaals
    >>> Tc, Pc = 304.13, 7.3773e6
    >>> vdw = VanDerWaals.from_critical_constants(Tc, Pc)
    >>> round(saturation_pressure(vdw, 0.9 * Tc) / Pc, 3)
    0.647
    """
    b = eos.b
    lnV = np.linspace(np.log(1.05 * b), np.log(200.0 * b), 4000)
    P_iso = eos.pressure(np.exp(lnV), T)
    dP = np.diff(P_iso)
    rising = np.where(dP > 0)[0]
    if rising.size == 0:
        raise ValueError("isotherm has no van der Waals loop; T is at or above the critical temperature")
    first_rise = int(rising[0])
    later_fall = np.where((dP < 0) & (np.arange(dP.size) > first_rise))[0]
    i_max = int(later_fall[0]) if later_fall.size else lnV.size - 1

    def refine(i, sign):
        lo, hi = lnV[max(i - 1, 0)], lnV[min(i + 1, lnV.size - 1)]
        res = minimize_scalar(lambda x: sign * eos.pressure(np.exp(x), T), bounds=(lo, hi), method="bounded")
        return float(eos.pressure(np.exp(res.x), T))

    P_low = max(refine(first_rise, 1.0), 0.0)
    P_high = refine(i_max, -1.0)
    span = P_high - P_low
    lo = P_low + 1e-6 * span if P_low > 0 else 1e-6 * P_high
    hi = P_high - 1e-6 * span

    def g(P):
        return np.log(fugacity_coefficient(eos, P, T, "liquid")) - np.log(fugacity_coefficient(eos, P, T, "vapor"))

    return float(brentq(g, lo, hi, xtol=1e-10 * P_high, rtol=1e-12))
