r"""Raoult's/Henry's law mixtures and colligative properties.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 5 ("Simple
mixtures") throughout: Sec. 5.4 for Raoult's and Henry's laws, Sec. 5.5
for the colligative properties (freezing-point depression, boiling-point
elevation, osmotic pressure).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.constants import R

__all__ = [
    "raoult_vapor_pressure",
    "BinaryIdealSolution",
    "henry_law_pressure",
    "CRYOSCOPIC_CONSTANTS",
    "freezing_point_depression",
    "boiling_point_elevation",
    "osmotic_pressure",
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
