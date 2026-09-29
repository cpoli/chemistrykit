r"""Raoult's/Henry's law mixtures and colligative properties.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 5 ("Simple
mixtures") throughout: Sec. 5.4 for Raoult's and Henry's laws, Sec. 5.5
for the colligative properties (freezing-point depression, boiling-point
elevation, osmotic pressure). The non-ideal activity-coefficient models
(Margules, Wilson, NRTL, UNIQUAC) share one isothermal modified-Raoult's-law
VLE implementation, :class:`BinaryActivityModel`; see Prausnitz,
Lichtenthaler & de Azevedo, *Molecular Thermodynamics of Fluid-Phase
Equilibria*, 3rd ed., Ch. 6-7, for all four.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.optimize import brentq

from chemistrykit.constants import R

__all__ = [
    "raoult_vapor_pressure",
    "BinaryIdealSolution",
    "henry_law_pressure",
    "CRYOSCOPIC_CONSTANTS",
    "freezing_point_depression",
    "boiling_point_elevation",
    "osmotic_pressure",
    "BinaryActivityModel",
    "MargulesSolution",
    "WilsonSolution",
    "NRTLSolution",
    "UNIQUACSolution",
]


def raoult_vapor_pressure(x, P_pure: float):
    r"""Raoult's law: :math:`P_A = x_A P_A^*`.

    The partial vapor pressure of a component in an ideal mixture is
    proportional to its mole fraction, with the pure-component vapor
    pressure as the constant of proportionality (F. M. Raoult, 1887;
    Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 5.4a).

    Parameters
    ----------
    x : float or array-like of float
        Mole fraction of the component in the liquid.
    P_pure : float
        Vapor pressure of the pure component, :math:`P_A^*`.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> round(float(raoult_vapor_pressure(x=0.4, P_pure=100.0)), 6)
    40.0
    """
    x = np.asarray(x, dtype=np.float64)
    return x * P_pure


@dataclass
class BinaryIdealSolution:
    r"""A binary liquid mixture obeying Raoult's law for both components.

    Combining Raoult's law for each component with Dalton's law for the
    vapor phase gives the total vapor pressure as a function of liquid
    composition, and the vapor-phase composition in equilibrium with it
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 5.4a):

    .. math::

        P_{tot}(x_A) = x_A P_A^* + (1-x_A) P_B^*, \qquad
        y_A = \frac{x_A P_A^*}{P_{tot}(x_A)}

    Parameters
    ----------
    P_A_star : float
        Vapor pressure of pure A.
    P_B_star : float
        Vapor pressure of pure B.
    """

    P_A_star: float
    P_B_star: float

    def total_pressure(self, x_A):
        """Total vapor pressure above a liquid of composition `x_A`.

        Parameters
        ----------
        x_A : float or array-like of float
            Liquid-phase mole fraction of A.

        Returns
        -------
        float or ndarray
        """
        x_A = np.asarray(x_A, dtype=np.float64)
        return x_A * self.P_A_star + (1.0 - x_A) * self.P_B_star

    def vapor_composition(self, x_A):
        r"""Vapor-phase mole fraction of A, :math:`y_A`, in equilibrium with liquid composition `x_A`.

        Parameters
        ----------
        x_A : float or array-like of float
            Liquid-phase mole fraction of A.

        Returns
        -------
        float or ndarray

        Examples
        --------
        For a 1:1 mixture of equally volatile components, the vapor has
        the same composition as the liquid:

        >>> solution = BinaryIdealSolution(P_A_star=100.0, P_B_star=100.0)
        >>> round(float(solution.vapor_composition(0.5)), 6)
        0.5

        The more volatile component (higher pure vapor pressure) is
        enriched in the vapor relative to the liquid:

        >>> solution = BinaryIdealSolution(P_A_star=200.0, P_B_star=50.0)
        >>> bool(solution.vapor_composition(0.5) > 0.5)
        True
        """
        x_A = np.asarray(x_A, dtype=np.float64)
        return x_A * self.P_A_star / self.total_pressure(x_A)


def henry_law_pressure(x, K_H: float):
    r"""Henry's law: :math:`P_B = x_B K_H`.

    For a dilute solute B, the partial vapor pressure is proportional to
    its mole fraction with the (empirical, solute- and solvent-specific)
    Henry's law constant :math:`K_H` as the proportionality constant --
    unlike Raoult's law, :math:`K_H \neq P_B^*` in general, because the
    solute's local environment in dilute solution (surrounded by solvent)
    differs from pure solute (W. Henry, *Philos. Trans. R. Soc.* 93, 29
    (1803); Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 5.4b).

    Parameters
    ----------
    x : float or array-like of float
        Mole fraction of the (dilute) solute in solution.
    K_H : float
        Henry's law constant for this solute/solvent pair, same pressure
        units as the desired result.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> round(float(henry_law_pressure(x=0.001, K_H=1.5e5)), 3)
    150.0
    """
    x = np.asarray(x, dtype=np.float64)
    return x * K_H


#: Cryoscopic (Kf) and ebullioscopic (Kb) constants for common solvents,
#: in K kg/mol (i.e. per unit *molality*), plus each solvent's normal
#: freezing/boiling point in K -- commonly tabulated literature values
#: (Atkins & de Paula, *Physical Chemistry*, 11th ed., Table 5.4). Provided
#: for convenience; :func:`freezing_point_depression` and
#: :func:`boiling_point_elevation` also accept `Kf`/`Kb` directly for any
#: solvent not listed here.
CRYOSCOPIC_CONSTANTS = {
    "water": {"Kf": 1.86, "Kb": 0.51, "Tf": 273.15, "Tb": 373.15},
    "benzene": {"Kf": 5.12, "Kb": 2.53, "Tf": 278.65, "Tb": 353.25},
    "camphor": {"Kf": 37.7, "Kb": None, "Tf": 452.65, "Tb": None},
    "acetic_acid": {"Kf": 3.90, "Kb": 3.07, "Tf": 289.75, "Tb": 391.05},
}


def freezing_point_depression(Kf: float, b: float, i: float = 1.0) -> float:
    r"""Freezing-point depression :math:`\Delta T_f = i K_f b`.

    A colligative property: it depends on the *number* of dissolved
    solute particles per kg of solvent, not their identity (Atkins & de
    Paula, *Physical Chemistry*, 11th ed., Ch. 5.5). See also
    :mod:`chemistrykit.solutions`, which uses the same van't Hoff factor
    `i` convention for strong-electrolyte dissociation, but does not
    duplicate this colligative-property machinery -- it lives here in
    `chemistrykit.thermo` only.

    Parameters
    ----------
    Kf : float
        Cryoscopic constant of the solvent, in K kg/mol (see
        :data:`CRYOSCOPIC_CONSTANTS` for common solvents).
    b : float
        Molality of solute, in mol/kg.
    i : float, default 1.0
        Van't Hoff factor: the number of particles each formula unit of
        solute dissociates into (1 for a nonelectrolyte, 2 for e.g. NaCl
        assuming complete dissociation).

    Returns
    -------
    float
        Freezing-point depression, in K (the new freezing point is
        :math:`T_f^* - \Delta T_f`).

    Examples
    --------
    A 1.00 molal aqueous NaCl solution (i=2, ideal complete dissociation):

    >>> round(freezing_point_depression(Kf=1.86, b=1.00, i=2.0), 2)
    3.72
    """
    return i * Kf * b


def boiling_point_elevation(Kb: float, b: float, i: float = 1.0) -> float:
    r"""Boiling-point elevation :math:`\Delta T_b = i K_b b`.

    Parameters
    ----------
    Kb : float
        Ebullioscopic constant of the solvent, in K kg/mol (see
        :data:`CRYOSCOPIC_CONSTANTS` for common solvents).
    b : float
        Molality of solute, in mol/kg.
    i : float, default 1.0
        Van't Hoff factor (see :func:`freezing_point_depression`).

    Returns
    -------
    float
        Boiling-point elevation, in K.

    Examples
    --------
    >>> round(boiling_point_elevation(Kb=0.51, b=1.00, i=2.0), 2)
    1.02
    """
    return i * Kb * b


def osmotic_pressure(M: float, T: float, i: float = 1.0, R_gas: float = R) -> float:
    r"""The van't Hoff equation for osmotic pressure, :math:`\Pi = i M R T`.

    Formally identical to the ideal gas law -- van't Hoff originally
    noted the (coincidental, but pedagogically useful) analogy (J. H.
    van't Hoff, 1887; Atkins & de Paula, *Physical Chemistry*, 11th ed.,
    Ch. 5.5e).

    Parameters
    ----------
    M : float
        Molarity of solute, in mol/m^3 (use ``M_mol_per_L * 1000`` to
        convert from the more commonly tabulated mol/L).
    T : float
        Absolute temperature, in K.
    i : float, default 1.0
        Van't Hoff factor (see :func:`freezing_point_depression`).
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1; with `M` in mol/m^3 this gives
        `Pi` in Pa.

    Returns
    -------
    float
        Osmotic pressure, in Pa.

    Examples
    --------
    A 0.100 mol/L (=100 mol/m^3) nonelectrolyte solution at 298.15 K:

    >>> round(osmotic_pressure(M=100.0, T=298.15) / 1000.0, 2)
    247.9
    """
    return i * M * R_gas * T


class BinaryActivityModel(ABC):
    r"""Isothermal vapor-liquid equilibrium of a binary liquid from an activity-coefficient model.

    Subclasses supply :meth:`activity_coefficients` and the pure-component
    vapor pressures ``P1_star``/``P2_star`` at the temperature of
    interest (e.g. from :class:`chemistrykit.thermo.AntoineEquation`);
    this base turns them into the modified Raoult's law
    :math:`P_i = x_i\gamma_iP_i^*` with an ideal vapor, i.e. a P-x-y
    diagram (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 5.3).
    """

    P1_star: float
    P2_star: float

    @abstractmethod
    def activity_coefficients(self, x1):
        r"""Return :math:`(\gamma_1, \gamma_2)` at liquid mole fraction `x1`."""

    def excess_gibbs(self, x1, T: float, R_gas: float = R):
        r"""Molar excess Gibbs energy :math:`G^E = RT(x_1\ln\gamma_1 + x_2\ln\gamma_2)`, in J/mol.

        Parameters
        ----------
        x1 : float or array-like of float
            Liquid-phase mole fraction of component 1.
        T : float
            Absolute temperature, in K.
        R_gas : float, default :data:`chemistrykit.constants.R`
            Gas constant.

        Returns
        -------
        float or ndarray
        """
        x1 = np.asarray(x1, dtype=np.float64)
        g1, g2 = self.activity_coefficients(x1)
        return R_gas * T * (x1 * np.log(g1) + (1.0 - x1) * np.log(g2))

    def partial_pressures(self, x1):
        r"""Return :math:`(P_1, P_2)` from the modified Raoult's law :math:`P_i = x_i\gamma_iP_i^*`.

        Parameters
        ----------
        x1 : float or array-like of float
            Liquid-phase mole fraction of component 1.

        Returns
        -------
        tuple of (float or ndarray)
        """
        x1 = np.asarray(x1, dtype=np.float64)
        g1, g2 = self.activity_coefficients(x1)
        return x1 * g1 * self.P1_star, (1.0 - x1) * g2 * self.P2_star

    def total_pressure(self, x1):
        """Total (bubble-point) vapor pressure above a liquid of composition `x1`.

        Parameters
        ----------
        x1 : float or array-like of float
            Liquid-phase mole fraction of component 1.

        Returns
        -------
        float or ndarray
        """
        P1, P2 = self.partial_pressures(x1)
        return P1 + P2

    def vapor_composition(self, x1):
        """Vapor-phase mole fraction of component 1 in equilibrium with liquid `x1`.

        Parameters
        ----------
        x1 : float or array-like of float
            Liquid-phase mole fraction of component 1.

        Returns
        -------
        float or ndarray
        """
        P1, P2 = self.partial_pressures(x1)
        return P1 / (P1 + P2)

    def relative_volatility(self, x1):
        r"""Relative volatility :math:`\alpha_{12} = \gamma_1P_1^*/(\gamma_2P_2^*)`.

        Parameters
        ----------
        x1 : float or array-like of float
            Liquid-phase mole fraction of component 1.

        Returns
        -------
        float or ndarray
        """
        g1, g2 = self.activity_coefficients(x1)
        return g1 * self.P1_star / (g2 * self.P2_star)

    def azeotrope(self) -> Optional[tuple[float, float]]:
        r"""Locate an azeotrope, where vapor and liquid compositions coincide.

        At an azeotrope :math:`y_1 = x_1`, i.e. the relative volatility is
        1. One exists in :math:`0 < x_1 < 1` exactly when
        :math:`\ln\alpha_{12}` changes sign across the composition range,
        and is found by root bracketing (:func:`scipy.optimize.brentq`).

        Returns
        -------
        (x1, P) or None
            Azeotropic composition and pressure, or ``None`` if there is
            no sign change (no azeotrope, or an even number of them).
        """
        eps = 1e-9

        def f(x1):
            return float(np.log(self.relative_volatility(x1)))

        if f(eps) * f(1.0 - eps) >= 0.0:
            return None
        x_az = brentq(f, eps, 1.0 - eps, xtol=1e-12)
        return x_az, float(self.total_pressure(x_az))


@dataclass
class MargulesSolution(BinaryActivityModel):
    r"""A non-ideal binary liquid described by the one-parameter (two-suffix) Margules model.

    Margules expanded the logarithms of the activity coefficients as power
    series in mole fraction (M. Margules, *Sitzungsber. Kais. Akad. Wiss.
    Wien, Math.-Naturwiss. Kl.* 104, 1243-1278 (1895)). Truncated at its
    first term, the model has a symmetric excess Gibbs energy
    :math:`G^E/(RT) = A\,x_1x_2`, giving

    .. math::

        \ln\gamma_1 = A\,x_2^2, \qquad \ln\gamma_2 = A\,x_1^2

    and the modified Raoult's law :math:`P_i = x_i\gamma_i P_i^*` (Atkins &
    de Paula, *Physical Chemistry*, 11th ed., Ch. 5.3). :math:`A > 0` gives
    positive deviations from Raoult's law, :math:`A < 0` negative ones,
    and :math:`A = 0` recovers :class:`BinaryIdealSolution`.

    Parameters
    ----------
    A : float
        Dimensionless Margules parameter, :math:`W/(RT)`.
    P1_star : float
        Vapor pressure of pure component 1.
    P2_star : float
        Vapor pressure of pure component 2.

    Examples
    --------
    At infinite dilution component 1 obeys Henry's law with
    :math:`K_H = P_1^* e^{A}`, while near :math:`x_1 = 1` its activity
    coefficient tends to 1 (Raoult's law):

    >>> import numpy as np
    >>> sol = MargulesSolution(A=1.2, P1_star=30.0, P2_star=20.0)
    >>> round(sol.henry_constant_1 / 30.0, 6) == round(float(np.exp(1.2)), 6)
    True
    >>> g1, g2 = sol.activity_coefficients(1.0)
    >>> float(g1), round(float(g2), 6) == round(float(np.exp(1.2)), 6)
    (1.0, True)
    """

    A: float
    P1_star: float
    P2_star: float

    def activity_coefficients(self, x1):
        r"""Return :math:`(\gamma_1, \gamma_2)` at liquid mole fraction `x1`.

        Parameters
        ----------
        x1 : float or array-like of float
            Liquid-phase mole fraction of component 1.

        Returns
        -------
        tuple of (float or ndarray)
        """
        x1 = np.asarray(x1, dtype=np.float64)
        x2 = 1.0 - x1
        return np.exp(self.A * x2**2), np.exp(self.A * x1**2)

    @property
    def henry_constant_1(self) -> float:
        r"""Henry's-law constant of component 1 at infinite dilution, :math:`P_1^* e^{A}`."""
        return float(self.P1_star * np.exp(self.A))


@dataclass
class WilsonSolution(BinaryActivityModel):
    r"""Wilson's local-composition model of a non-ideal binary liquid.

    G. M. Wilson, *J. Am. Chem. Soc.* 86, 127 (1964):

    .. math::

        \ln\gamma_1 = -\ln(x_1+\Lambda_{12}x_2)
            + x_2\left(\frac{\Lambda_{12}}{x_1+\Lambda_{12}x_2}
            - \frac{\Lambda_{21}}{x_2+\Lambda_{21}x_1}\right)

    and symmetrically for :math:`\gamma_2`. Wilson's model cannot predict
    liquid-liquid phase splitting, but is accurate for miscible mixtures
    of polar and non-polar components (e.g. alcohols in hydrocarbons),
    including their azeotropes.

    Parameters
    ----------
    Lambda12, Lambda21 : float
        Positive, dimensionless Wilson parameters at the temperature of
        interest (see :meth:`from_energies`).
    P1_star, P2_star : float
        Pure-component vapor pressures at that temperature.

    Examples
    --------
    With :math:`\Lambda_{12}=\Lambda_{21}=1` the mixture is ideal:

    >>> sol = WilsonSolution(Lambda12=1.0, Lambda21=1.0, P1_star=30.0, P2_star=20.0)
    >>> [round(float(g), 12) for g in sol.activity_coefficients(0.3)]
    [1.0, 1.0]
    """

    Lambda12: float
    Lambda21: float
    P1_star: float
    P2_star: float

    @classmethod
    def from_energies(cls, V1: float, V2: float, lambda12: float, lambda21: float, T: float, P1_star: float, P2_star: float, R_gas: float = R):
        r"""Build from liquid molar volumes and interaction energies, :math:`\Lambda_{12}=(V_2/V_1)e^{-\lambda_{12}/RT}`.

        Parameters
        ----------
        V1, V2 : float
            Pure-liquid molar volumes (any consistent unit).
        lambda12, lambda21 : float
            Interaction-energy differences :math:`\lambda_{12}-\lambda_{11}`
            and :math:`\lambda_{21}-\lambda_{22}`, in J/mol.
        T : float
            Absolute temperature, in K.
        P1_star, P2_star : float
            Pure-component vapor pressures at `T`.
        R_gas : float, default :data:`chemistrykit.constants.R`

        Returns
        -------
        WilsonSolution
        """
        Lambda12 = V2 / V1 * np.exp(-lambda12 / (R_gas * T))
        Lambda21 = V1 / V2 * np.exp(-lambda21 / (R_gas * T))
        return cls(float(Lambda12), float(Lambda21), P1_star, P2_star)

    def activity_coefficients(self, x1):
        r"""Return :math:`(\gamma_1, \gamma_2)` at liquid mole fraction `x1`.

        Parameters
        ----------
        x1 : float or array-like of float

        Returns
        -------
        tuple of (float or ndarray)
        """
        x1 = np.asarray(x1, dtype=np.float64)
        x2 = 1.0 - x1
        d1 = x1 + self.Lambda12 * x2
        d2 = x2 + self.Lambda21 * x1
        bracket = self.Lambda12 / d1 - self.Lambda21 / d2
        return np.exp(-np.log(d1) + x2 * bracket), np.exp(-np.log(d2) - x1 * bracket)


@dataclass
class NRTLSolution(BinaryActivityModel):
    r"""The non-random two-liquid (NRTL) model of a non-ideal binary liquid.

    H. Renon & J. M. Prausnitz, *AIChE J.* 14, 135 (1968). With
    :math:`G_{ij}=\exp(-\alpha\tau_{ij})`:

    .. math::

        \ln\gamma_1 = x_2^2\left[\tau_{21}\left(\frac{G_{21}}{x_1+x_2G_{21}}\right)^2
            + \frac{\tau_{12}G_{12}}{(x_2+x_1G_{12})^2}\right]

    and symmetrically for :math:`\gamma_2`. Unlike Wilson's model, NRTL
    can describe partially miscible (liquid-liquid) systems.

    Parameters
    ----------
    tau12, tau21 : float
        Dimensionless interaction parameters at the temperature of
        interest (commonly :math:`\tau_{ij}=b_{ij}/RT`).
    P1_star, P2_star : float
        Pure-component vapor pressures at that temperature.
    alpha : float, default 0.3
        Non-randomness parameter (0.2-0.47 in Renon & Prausnitz's fits).

    Examples
    --------
    At infinite dilution, :math:`\ln\gamma_1^\infty=\tau_{21}+\tau_{12}e^{-\alpha\tau_{12}}`:

    >>> import numpy as np
    >>> sol = NRTLSolution(tau12=0.5, tau21=1.2, P1_star=30.0, P2_star=20.0)
    >>> g1, _ = sol.activity_coefficients(0.0)
    >>> bool(np.isclose(np.log(g1), 1.2 + 0.5 * np.exp(-0.3 * 0.5)))
    True
    """

    tau12: float
    tau21: float
    P1_star: float
    P2_star: float
    alpha: float = 0.3

    def activity_coefficients(self, x1):
        r"""Return :math:`(\gamma_1, \gamma_2)` at liquid mole fraction `x1`.

        Parameters
        ----------
        x1 : float or array-like of float

        Returns
        -------
        tuple of (float or ndarray)
        """
        x1 = np.asarray(x1, dtype=np.float64)
        x2 = 1.0 - x1
        G12 = np.exp(-self.alpha * self.tau12)
        G21 = np.exp(-self.alpha * self.tau21)
        d1 = x1 + x2 * G21
        d2 = x2 + x1 * G12
        ln_g1 = x2**2 * (self.tau21 * (G21 / d1) ** 2 + self.tau12 * G12 / d2**2)
        ln_g2 = x1**2 * (self.tau12 * (G12 / d2) ** 2 + self.tau21 * G21 / d1**2)
        return np.exp(ln_g1), np.exp(ln_g2)


@dataclass
class UNIQUACSolution(BinaryActivityModel):
    r"""The universal quasi-chemical (UNIQUAC) model of a non-ideal binary liquid.

    D. S. Abrams & J. M. Prausnitz, *AIChE J.* 21, 116 (1975).
    :math:`\ln\gamma_i` is the sum of a combinatorial (size/shape) part,
    fixed by each molecule's van der Waals volume `r` and surface area
    `q`, and a residual (energetic) part:

    .. math::

        \ln\gamma_i^C = \ln\frac{\Phi_i}{x_i} + \frac z2 q_i\ln\frac{\theta_i}{\Phi_i}
            + l_i - \frac{\Phi_i}{x_i}\sum_j x_jl_j, \qquad
        \ln\gamma_i^R = q_i\left[1-\ln\sum_j\theta_j\tau_{ji}
            - \sum_j\frac{\theta_j\tau_{ij}}{\sum_k\theta_k\tau_{kj}}\right]

    with :math:`\Phi_i=r_ix_i/\sum_jr_jx_j`,
    :math:`\theta_i=q_ix_i/\sum_jq_jx_j`,
    :math:`l_i=\frac z2(r_i-q_i)-(r_i-1)`, :math:`\tau_{ii}=1` and
    coordination number :math:`z=10`.

    Parameters
    ----------
    r1, r2, q1, q2 : float
        Relative van der Waals volumes and surface areas (tabulated, e.g.
        water r=0.92, q=1.40; ethanol r=2.1055, q=1.972).
    tau12, tau21 : float
        Dimensionless energy parameters, :math:`\tau_{ij}=\exp(-a_{ij}/T)`.
    P1_star, P2_star : float
        Pure-component vapor pressures at the temperature of interest.
    z : float, default 10.0
        Lattice coordination number.

    Examples
    --------
    Identical-size molecules with :math:`\tau_{12}=\tau_{21}=1` form an ideal solution:

    >>> sol = UNIQUACSolution(r1=1.0, r2=1.0, q1=1.0, q2=1.0, tau12=1.0, tau21=1.0, P1_star=30.0, P2_star=20.0)
    >>> [round(float(g), 12) for g in sol.activity_coefficients(0.3)]
    [1.0, 1.0]
    """

    r1: float
    r2: float
    q1: float
    q2: float
    tau12: float
    tau21: float
    P1_star: float
    P2_star: float
    z: float = 10.0

    def activity_coefficients(self, x1):
        r"""Return :math:`(\gamma_1, \gamma_2)` at liquid mole fraction `x1`.

        Parameters
        ----------
        x1 : float or array-like of float

        Returns
        -------
        tuple of (float or ndarray)
        """
        x1 = np.asarray(x1, dtype=np.float64)
        x2 = 1.0 - x1
        r1, r2, q1, q2, z = self.r1, self.r2, self.q1, self.q2, self.z
        rx = r1 * x1 + r2 * x2
        qx = q1 * x1 + q2 * x2
        l1 = z / 2.0 * (r1 - q1) - (r1 - 1.0)
        l2 = z / 2.0 * (r2 - q2) - (r2 - 1.0)
        lx = x1 * l1 + x2 * l2
        theta1, theta2 = q1 * x1 / qx, q2 * x2 / qx
        t12, t21 = self.tau12, self.tau21
        s1 = theta1 + theta2 * t21
        s2 = theta1 * t12 + theta2
        out = []
        for r, q, l, s_own, resid in (
            (r1, q1, l1, s1, theta1 / s1 + theta2 * t12 / s2),
            (r2, q2, l2, s2, theta1 * t21 / s1 + theta2 / s2),
        ):
            phi_over_x = r / rx
            theta_over_phi = (q / qx) / phi_over_x
            ln_comb = np.log(phi_over_x) + z / 2.0 * q * np.log(theta_over_phi) + l - phi_over_x * lx
            ln_res = q * (1.0 - np.log(s_own) - resid)
            out.append(np.exp(ln_comb + ln_res))
        return out[0], out[1]
