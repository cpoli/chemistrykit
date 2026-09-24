r"""Redox and complexometric (EDTA) titration curves.

Acid-base titration curves already exist in
:mod:`chemistrykit.solutions.systems.titration`
(:class:`~chemistrykit.solutions.systems.titration.StrongAcidStrongBaseTitration`,
:class:`~chemistrykit.solutions.systems.titration.WeakAcidStrongBaseTitration`,
:class:`~chemistrykit.solutions.systems.titration.WeakBaseStrongAcidTitration`)
and are reused directly rather than duplicated here (see this domain's
examples gallery, which exercises all three titration *types* -- acid-
base, redox, complexometric -- side by side). This module adds the two
that are new to analytical chemistry: potentiometric redox titrations and
EDTA (complexometric) titrations, both implementing
:class:`chemistrykit.analytical.core.base_system.TitrationCurve`.

See Harris, *Quantitative Chemical Analysis*, 9th ed., Ch. 16 (redox
titrations) and Ch. 11-12 (EDTA titrations), or Skoog, West, Holler &
Crouch, *Fundamentals of Analytical Chemistry*, 9th ed., Ch. 16 and 17,
throughout.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.analytical.core.base_system import TitrationCurve
from chemistrykit.constants import FARADAY, R

__all__ = ["RedoxTitration", "EDTATitration", "GranPlotResult", "gran_plot"]


class RedoxTitration(TitrationCurve):
    r"""A potentiometric redox titration: a reduced analyte titrated by an oxidizing titrant.

    Models the general reaction :math:`n_2\,\text{Red}_1 + n_1\,\text{Ox}_2
    \rightarrow n_2\,\text{Ox}_1 + n_1\,\text{Red}_2` in the standard
    large-equilibrium-constant ("quantitative reaction") approximation
    (Harris, *Quantitative Chemical Analysis*, 9th ed., Ch. 16), under
    which electron-balance stoichiometry alone (not a numerical
    equilibrium solve) fixes the composition at every titrant volume:

    * **Before the equivalence volume** :math:`V_{eq}=n_1C_{analyte}V_{analyte}/(n_2C_{titrant})`,
      the analyte couple sets the potential via the Nernst equation, using
      the fraction of analyte already oxidized (fixed by the moles of
      titrant added, via the electron balance :math:`n_1x=n_2C_{titrant}V`):

      .. math::

          E = E^\circ_1 - \frac{RT}{n_1F}\ln\frac{[\text{Red}_1]}{[\text{Ox}_1]}

    * **At the equivalence point**, the classical weighted-average result
      (from adding the two Nernst equations):

      .. math::

          E_{eq} = \frac{n_1E^\circ_1+n_2E^\circ_2}{n_1+n_2}

    * **Past the equivalence volume**, the titrant couple sets the
      potential from the (now excess) unreacted titrant:

      .. math::

          E = E^\circ_2 - \frac{RT}{n_2F}\ln\frac{[\text{Red}_2]}{[\text{Ox}_2]}

    For the common case :math:`n_1=n_2` (e.g. Fe2+/Ce4+, both
    1-electron couples) this reduces to the textbook symmetric titration
    curve with :math:`E_{eq}=(E^\circ_1+E^\circ_2)/2`.

    Parameters
    ----------
    E1_standard : float
        Standard reduction potential of the analyte couple (Ox1/Red1), in V.
    n1 : int
        Electrons transferred in the analyte half-reaction.
    E2_standard : float
        Standard reduction potential of the titrant couple (Ox2/Red2), in V.
    n2 : int
        Electrons transferred in the titrant half-reaction.
    C_analyte : float
        Initial analyte (Red1) concentration, in mol/L.
    V_analyte : float
        Initial analyte volume, in L.
    C_titrant : float
        Titrant (Ox2) concentration, in mol/L.
    T : float, default :data:`chemistrykit.constants.STANDARD_TEMPERATURE`
        Absolute temperature, in K.

    Examples
    --------
    Fe2+ (E=0.771 V) titrated by Ce4+ (E=1.72 V), both 1-electron
    couples: at the equivalence volume, `E` is exactly the arithmetic
    mean of the two standard potentials -- the textbook symmetric-curve
    result:

    >>> titration = RedoxTitration(E1_standard=0.771, n1=1, E2_standard=1.72, n2=1, C_analyte=0.10, V_analyte=0.050, C_titrant=0.10)
    >>> V_eq = titration.equivalence_volume()
    >>> round(float(titration.response_at(np.array([V_eq]))[0]), 4) == round((0.771 + 1.72) / 2.0, 4)
    True
    """

    def __init__(self, E1_standard: float, n1: int, E2_standard: float, n2: int, C_analyte: float, V_analyte: float, C_titrant: float, T: float = 298.15):
        self.E1_standard = float(E1_standard)
        self.n1 = int(n1)
        self.E2_standard = float(E2_standard)
        self.n2 = int(n2)
        self.C_analyte = float(C_analyte)
        self.V_analyte = float(V_analyte)
        self.C_titrant = float(C_titrant)
        self.T = float(T)

    def equivalence_volume(self) -> float:
        """Return the exact (stoichiometric, electron-balance) equivalence volume.

        Returns
        -------
        float
        """
        return self.n1 * self.C_analyte * self.V_analyte / (self.n2 * self.C_titrant)

    def response_at(self, V):
        V = np.atleast_1d(np.asarray(V, dtype=np.float64))
        V_eq = self.equivalence_volume()
        RT_n1F = R * self.T / (self.n1 * FARADAY)
        RT_n2F = R * self.T / (self.n2 * FARADAY)
        a = self.C_analyte * self.V_analyte
        E = np.empty_like(V)
        for i, v in enumerate(V):
            # The isclose check must come first: a `v` that is a tiny
            # floating-point epsilon *below* V_eq (e.g. from an
            # independently-computed volume array, or plain round-off)
            # would otherwise fall into the `v < V_eq` branch below, where
            # `x` clamps to within 1e-12 of `a` regardless of how close
            # `v` actually is to V_eq -- silently returning a physically
            # meaningless clamped potential instead of the correct
            # equivalence-point value, for a `v` indistinguishable from
            # V_eq at any reasonable tolerance.
            if np.isclose(v, V_eq, rtol=1e-9, atol=1e-15):
                E[i] = (self.n1 * self.E1_standard + self.n2 * self.E2_standard) / (self.n1 + self.n2)
            elif v < V_eq:
                x = (self.n2 / self.n1) * self.C_titrant * v  # mol analyte oxidized
                x = min(x, a * (1.0 - 1e-12))
                E[i] = self.E1_standard - RT_n1F * np.log((a - x) / x)
            else:
                E[i] = self.E2_standard + RT_n2F * np.log((v - V_eq) / V_eq)
        return E


class EDTATitration(TitrationCurve):
    r"""A complexometric titration of a metal ion M by EDTA, via the conditional formation constant.

    Models the 1:1 complexation :math:`M+Y\rightleftharpoons MY` exactly
    (no large-K approximation needed, unlike :class:`RedoxTitration`),
    using the *conditional* formation constant :math:`K_f'` -- the
    effective (pH- and side-reaction-corrected) formation constant at the
    titration's fixed pH (Harris, *Quantitative Chemical Analysis*, 9th
    ed., Ch. 12.2; folding pH-dependence into a single conditional
    constant, rather than modeling EDTA's stepwise protonation
    equilibria explicitly, is itself the standard simplification this
    class makes, flagged here as an approximation).

    Mass balance on total metal :math:`C_M` and total EDTA added
    :math:`C_Y` (both diluted by the growing total volume), with
    :math:`x=[MY]`, :math:`[M]=C_M-x`, :math:`[Y]=C_Y-x`, and
    :math:`K_f'=x/[(C_M-x)(C_Y-x)]`, gives the exact quadratic

    .. math::

        K_f'x^2-\left[K_f'(C_M+C_Y)+1\right]x+K_f'C_MC_Y=0

    solved here in closed form for the physical root :math:`0\le x\le
    \min(C_M,C_Y)`, then :math:`pM=-\log_{10}[M]`.

    Parameters
    ----------
    C_metal : float
        Initial metal-ion concentration, in mol/L.
    V_metal : float
        Initial metal-ion solution volume, in L.
    K_conditional : float
        Conditional formation constant :math:`K_f'` at the titration pH.
    C_edta : float
        Titrant (EDTA) concentration, in mol/L.

    Examples
    --------
    At the equivalence volume (:math:`C_M=C_Y`), `pM` approaches the
    textbook large-:math:`K_f'` approximation :math:`pM\approx
    \frac12\log_{10}(K_f'/C_{M,eq})`
    (Harris, *Quantitative Chemical Analysis*, 9th ed., eq. 12.2) as
    :math:`K_f'` grows -- verified here rather than assumed, since the
    approximation formula is itself only a limit of the exact quadratic
    solved above:

    >>> import numpy as np
    >>> titration = EDTATitration(C_metal=0.010, V_metal=0.050, K_conditional=1e12, C_edta=0.010)
    >>> V_eq = titration.equivalence_volume()
    >>> pM_exact = float(titration.response_at(np.array([V_eq]))[0])
    >>> C_M_eq = titration.C_metal * titration.V_metal / (titration.V_metal + V_eq)
    >>> pM_approx = 0.5 * np.log10(1e12 / C_M_eq)
    >>> bool(abs(pM_exact - pM_approx) < 0.01)
    True
    """

    def __init__(self, C_metal: float, V_metal: float, K_conditional: float, C_edta: float):
        self.C_metal = float(C_metal)
        self.V_metal = float(V_metal)
        self.K_conditional = float(K_conditional)
        self.C_edta = float(C_edta)

    def equivalence_volume(self) -> float:
        """Return the exact (stoichiometric, 1:1) equivalence volume.

        Returns
        -------
        float
        """
        return self.C_metal * self.V_metal / self.C_edta

    def _free_metal_concentration(self, C_M: float, C_Y: float) -> float:
        K = self.K_conditional
        a, b, c = K, -(K * (C_M + C_Y) + 1.0), K * C_M * C_Y
        disc = b * b - 4.0 * a * c
        sqrt_disc = np.sqrt(max(disc, 0.0))
        x_minus = (-b - sqrt_disc) / (2.0 * a)
        x_plus = (-b + sqrt_disc) / (2.0 * a)
        x_max = min(C_M, C_Y)
        for x in (x_minus, x_plus):
            if -1e-15 <= x <= x_max + 1e-15:
                return max(C_M - x, 1e-300)
        raise ValueError("no physical root found for the M/Y mass-action quadratic")

    def response_at(self, V):
        V = np.atleast_1d(np.asarray(V, dtype=np.float64))
        pM = np.empty_like(V)
        for i, v in enumerate(V):
            V_tot = self.V_metal + v
            C_M = self.C_metal * self.V_metal / V_tot
            C_Y = self.C_edta * v / V_tot
            free_M = self._free_metal_concentration(C_M, C_Y)
            pM[i] = -np.log10(free_M)
        return pM


@dataclass
class GranPlotResult:
    """Result of a :func:`gran_plot` call."""

    V: np.ndarray
    """ndarray: Titrant volumes used in the fit, in L."""

    gran_function: np.ndarray
    r"""ndarray: The Gran function :math:`V_b\,10^{-pH}` at each volume in `V`."""

    slope: float
    """float: Least-squares slope of the Gran function vs. `V` (equal to :math:`-K_a`)."""

    intercept: float
    """float: Least-squares intercept of the Gran function vs. `V`."""

    @property
    def equivalence_volume(self) -> float:
        """float: The extrapolated x-intercept, ``-intercept / slope``: the equivalence volume."""
        return -self.intercept / self.slope

    @property
    def Ka(self) -> float:
        """float: The acid dissociation constant implied by the slope, ``-slope``."""
        return -self.slope


def gran_plot(V, pH) -> GranPlotResult:
    r"""Gran plot for a weak acid titrated by a strong base: locate the equivalence point by linear extrapolation.

    Before the equivalence volume :math:`V_e`, the buffer region obeys
    :math:`[H^+]=K_a\,n_{HA}/n_{A^-}=K_a(V_e-V_b)/V_b`, so the Gran function

    .. math::

        V_b\,10^{-pH} = K_a\,(V_e - V_b)

    is a straight line in :math:`V_b` whose x-intercept is :math:`V_e` and
    whose slope is :math:`-K_a` (G. Gran, *Analyst* 77, 661 (1952)). This
    turns the equivalence point into a linear extrapolation of data taken
    well *before* it, rather than a search for the steepest point of the
    sigmoidal curve (activity coefficients are neglected here).

    Parameters
    ----------
    V : array-like of float
        Titrant volumes in the buffer region (before the equivalence point), in L.
    pH : array-like of float
        Measured pH at each volume.

    Returns
    -------
    GranPlotResult

    Examples
    --------
    Ideal buffer-region data built from the Gran relation itself recover
    :math:`V_e` and :math:`K_a` exactly:

    >>> import numpy as np
    >>> Ka, Ve = 1.8e-5, 0.050
    >>> V = np.linspace(0.010, 0.045, 8)
    >>> pH = -np.log10(Ka * (Ve - V) / V)
    >>> result = gran_plot(V, pH)
    >>> round(result.equivalence_volume, 9), round(result.Ka / Ka, 9)
    (0.05, 1.0)
    """
    V = np.asarray(V, dtype=np.float64)
    pH = np.asarray(pH, dtype=np.float64)
    G = V * 10.0 ** (-pH)
    slope, intercept = np.polyfit(V, G, 1)
    return GranPlotResult(V=V, gran_function=G, slope=float(slope), intercept=float(intercept))
