r"""Clausius-Clapeyron phase boundaries and the Gibbs phase rule.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 4 ("Physical
transformations of pure substances") for both topics.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.constants import R

__all__ = ["ClausiusClapeyron", "gibbs_phase_rule"]


@dataclass
class ClausiusClapeyron:
    r"""The (integrated) Clausius-Clapeyron equation for a liquid-vapor (or solid-vapor) boundary.

    Assuming the enthalpy of vaporization :math:`\Delta H_{vap}` is
    constant over the temperature range of interest, that the vapor
    behaves as an ideal gas, and that the molar volume of the condensed
    phase is negligible compared to the vapor's (Atkins & de Paula,
    *Physical Chemistry*, 11th ed., Ch. 4.1, eq. 4.11), the
    Clausius-Clapeyron differential equation :math:`dP/dT =
    \Delta H_{vap}P/(RT^2)` integrates to

    .. math::

        \ln\frac{P(T)}{P_{ref}} = -\frac{\Delta H_{vap}}{R}
            \left(\frac{1}{T} - \frac{1}{T_{ref}}\right)

    Parameters
    ----------
    delta_h_vap : float
        Enthalpy of vaporization, in J/mol (assumed temperature-independent
        -- an approximation whose quality degrades far from `T_ref`).
    T_ref : float
        A reference temperature at which the vapor pressure `P_ref` is known, in K.
    P_ref : float
        Vapor pressure at `T_ref`, in Pa.

    Examples
    --------
    Water: normal boiling point 373.15 K at 101325 Pa, using the
    (literature) enthalpy of vaporization at the normal boiling point,
    40700 J/mol (Atkins & de Paula, Table 4.1), to estimate the vapor
    pressure at 363.15 K (90 degC). The constant-:math:`\Delta H_{vap}`
    approximation gives 70.6 kPa, close to the accepted steam-table value
    of about 70.1 kPa -- the ~1% discrepancy is exactly the approximation
    error flagged above (:math:`\Delta H_{vap}` actually decreases
    somewhat as T rises toward the critical point):

    >>> water = ClausiusClapeyron(delta_h_vap=40700.0, T_ref=373.15, P_ref=101325.0)
    >>> round(float(water.pressure(363.15)) / 1000.0, 1)
    70.6
    """

    delta_h_vap: float
    T_ref: float
    P_ref: float
    R_gas: float = R

    def pressure(self, T):
        """Return the vapor pressure at temperature(s) `T`.

        Parameters
        ----------
        T : float or array-like of float
            Absolute temperature(s), in K.

        Returns
        -------
        float or ndarray
            Vapor pressure, in Pa (same units as `P_ref`).
        """
        T = np.asarray(T, dtype=np.float64)
        return self.P_ref * np.exp(-self.delta_h_vap / self.R_gas * (1.0 / T - 1.0 / self.T_ref))

    def boiling_point(self, P):
        """Return the temperature at which the vapor pressure equals `P` (invert :meth:`pressure`).

        Parameters
        ----------
        P : float or array-like of float
            Target vapor pressure, in Pa.

        Returns
        -------
        float or ndarray
            Temperature, in K.

        Examples
        --------
        Inverting :meth:`pressure` recovers the reference point exactly:

        >>> water = ClausiusClapeyron(delta_h_vap=40700.0, T_ref=373.15, P_ref=101325.0)
        >>> round(float(water.boiling_point(101325.0)), 2)
        373.15
        """
        P = np.asarray(P, dtype=np.float64)
        inv_T = 1.0 / self.T_ref - (self.R_gas / self.delta_h_vap) * np.log(P / self.P_ref)
        return 1.0 / inv_T

    @classmethod
    def from_two_points(cls, T1: float, P1: float, T2: float, P2: float) -> ClausiusClapeyron:
        r"""Determine :math:`\Delta H_{vap}` from two (T, P) points on the phase boundary.

        Solving the integrated equation for the slope between two known
        points:

        .. math::

            \Delta H_{vap} = \frac{-R \ln(P_2/P_1)}{1/T_2 - 1/T_1}

        Parameters
        ----------
        T1, P1 : float
            First reference point (K, Pa).
        T2, P2 : float
            Second point (K, Pa).

        Returns
        -------
        ClausiusClapeyron
            With `T_ref`, `P_ref` set to ``(T1, P1)``.

        Examples
        --------
        Recovering a known enthalpy of vaporization from two synthetic
        points generated with it:

        >>> model = ClausiusClapeyron(delta_h_vap=35000.0, T_ref=300.0, P_ref=1.0e4)
        >>> T2 = 320.0
        >>> P2 = model.pressure(T2)
        >>> fit = ClausiusClapeyron.from_two_points(300.0, 1.0e4, T2, P2)
        >>> round(fit.delta_h_vap, 2)
        35000.0
        """
        delta_h_vap = -R * np.log(P2 / P1) / (1.0 / T2 - 1.0 / T1)
        return cls(delta_h_vap=float(delta_h_vap), T_ref=T1, P_ref=P1)


def gibbs_phase_rule(n_components: int, n_phases: int, reactions: int = 0) -> int:
    r"""The Gibbs phase rule: :math:`F = C - P + 2 - r`.

    `F` is the number of intensive degrees of freedom (e.g. temperature,
    pressure, composition variables) that can be varied independently
    while the system remains in the same set of coexisting phases at
    equilibrium; `C` is the number of independent components, `P` the
    number of phases in equilibrium, and `r` the number of independent
    reaction equilibria constraining the composition (0 for a
    non-reactive system) (Atkins & de Paula, *Physical Chemistry*, 11th
    ed., Ch. 4.5; J. W. Gibbs, *Trans. Connecticut Acad.* 3, 108 (1876)).

    Parameters
    ----------
    n_components : int
        Number of independent chemical components, `C`.
    n_phases : int
        Number of phases in equilibrium, `P`.
    reactions : int, default 0
        Number of independent reaction-equilibrium constraints, `r`.

    Returns
    -------
    int
        Degrees of freedom, `F`.

    Examples
    --------
    A pure substance (C=1) at its triple point (three phases coexisting)
    has zero degrees of freedom -- the triple point is a single fixed
    (T, P):

    >>> gibbs_phase_rule(n_components=1, n_phases=3)
    0

    A pure substance with two phases in equilibrium (e.g. liquid-vapor)
    has one degree of freedom -- the phase boundary is a curve, P(T):

    >>> gibbs_phase_rule(n_components=1, n_phases=2)
    1
    """
    return n_components - n_phases + 2 - reactions
