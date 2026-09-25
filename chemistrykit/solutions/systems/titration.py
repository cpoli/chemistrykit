r"""Acid-base titration curves: strong/strong, weak-acid/strong-base, weak-base/strong-acid.

Every :class:`~chemistrykit.solutions.core.base_system.Titration`
subclass here computes pH as a function of titrant volume `Vb` from an
exact charge-balance equation -- see Harris, *Quantitative Chemical
Analysis*, 9th ed., Ch. 10-11, or Skoog, West, Holler & Crouch,
*Fundamentals of Analytical Chemistry*, 9th ed., Ch. 14, for the
derivations.

**Strong acid + strong base.** With :math:`\Delta C = (C_aV_a -
C_bV_b)/(V_a+V_b)` the net excess-acid concentration at a given titrant
volume, the charge balance :math:`[H^+] + C_bV_b/(V_a+V_b) = [OH^-] +
C_aV_a/(V_a+V_b)` rearranges (using :math:`[OH^-]=K_w/[H^+]`) to the
quadratic :math:`[H^+]^2 - \Delta C\,[H^+] - K_w = 0`, solved in closed
form.

**Weak acid + strong base.** Total acid concentration
:math:`C_{a,tot}=C_aV_a/(V_a+V_b)` is diluted as titrant is added, and
:math:`[A^-] = K_aC_{a,tot}/(K_a+[H^+])` from the equilibrium. The charge
balance :math:`C_bV_b/(V_a+V_b) + [H^+] = [A^-] + K_w/[H^+]` is then a
rational (not polynomial) equation in :math:`[H^+]`, solved numerically
via :func:`chemistrykit.solutions.utils.rootfinding.find_positive_root`
at each titrant volume. :class:`WeakBaseStrongAcidTitration` is the
mirror image.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.solutions.core.base_system import GranPlotResult, Titration
from chemistrykit.solutions.systems.acid_base import ph_from_h, poh_from_oh
from chemistrykit.solutions.utils.rootfinding import find_positive_root

__all__ = [
    "StrongAcidStrongBaseTitration",
    "WeakAcidStrongBaseTitration",
    "WeakBaseStrongAcidTitration",
    "gran_plot",
]


def _charge_balance_bracket(C_weak: float, C_strong: float, Kw: float) -> tuple:
    r"""Bracket ``[lo, hi]`` guaranteed to contain the root of a weak/strong titration charge balance.

    For :math:`r(x)=C_{strong}+x-K C_{weak}/(K+x)-K_w/x` (with `x` the
    :math:`[H^+]` of a weak-acid titration or the :math:`[OH^-]` of a
    weak-base one), `r` is increasing in `x`, and since the conjugate
    ion's concentration lies in :math:`[0, C_{weak}]`,
    :math:`r(hi)>0` at :math:`hi=2(C_{weak}+\sqrt{K_w})` and
    :math:`r(lo)<0` at :math:`lo=K_w/[2(C_{weak}+C_{strong}+\sqrt{K_w})]`,
    for any concentrations -- unlike a fixed pH 0-14 bracket, which
    concentrated solutions step outside of.
    """
    sqrt_Kw = np.sqrt(Kw)
    return Kw / (2.0 * (C_weak + C_strong + sqrt_Kw)), 2.0 * (C_weak + sqrt_Kw)


class StrongAcidStrongBaseTitration(Titration):
    r"""Titrating a strong acid (in the flask) with a strong base (the titrant).

    Parameters
    ----------
    Ca : float
        Initial strong-acid concentration, in mol/L.
    Va : float
        Initial acid volume, in L.
    Cb : float
        Titrant (strong base) concentration, in mol/L.
    Kw : float, default 1.0e-14
        Water autoionization constant.

    Examples
    --------
    The equivalence volume (where moles acid = moles base) is exactly
    ``Ca*Va/Cb``, and the numerically detected steepest-ascent point
    should match it closely:

    >>> import numpy as np
    >>> titration = StrongAcidStrongBaseTitration(Ca=0.100, Va=0.050, Cb=0.100)
    >>> Vb_equiv_exact = titration.equivalence_volume()
    >>> round(Vb_equiv_exact, 4)
    0.05
    >>> Vb = np.linspace(1e-6, 0.09, 20000)
    >>> Vb_equiv_numeric = titration.find_equivalence_point(Vb)
    >>> abs(Vb_equiv_numeric - Vb_equiv_exact) < 1e-3
    True
    """

    def __init__(self, Ca: float, Va: float, Cb: float, Kw: float = 1.0e-14):
        self.Ca = float(Ca)
        self.Va = float(Va)
        self.Cb = float(Cb)
        self.Kw = float(Kw)

    def equivalence_volume(self) -> float:
        """Return the exact (stoichiometric) equivalence volume, ``Ca*Va/Cb``.

        Returns
        -------
        float
        """
        return self.Ca * self.Va / self.Cb

    def pH_at(self, Vb):
        Vb = np.asarray(Vb, dtype=np.float64)
        V_tot = self.Va + Vb
        delta_C = (self.Ca * self.Va - self.Cb * Vb) / V_tot
        h = 0.5 * (delta_C + np.sqrt(delta_C**2 + 4.0 * self.Kw))
        return ph_from_h(h)


class WeakAcidStrongBaseTitration(Titration):
    r"""Titrating a weak acid HA (in the flask) with a strong base (the titrant).

    Parameters
    ----------
    Ca : float
        Initial weak-acid concentration, in mol/L.
    Va : float
        Initial acid volume, in L.
    Ka : float
        Acid dissociation constant of HA.
    Cb : float
        Titrant (strong base) concentration, in mol/L.
    Kw : float, default 1.0e-14
        Water autoionization constant.

    Examples
    --------
    At the half-equivalence volume, exactly half the acid has been
    converted to its conjugate base, so by Henderson-Hasselbalch
    ``pH = pKa`` there:

    >>> import numpy as np
    >>> titration = WeakAcidStrongBaseTitration(Ca=0.100, Va=0.050, Ka=1.8e-5, Cb=0.100)
    >>> Vb_half = titration.equivalence_volume() / 2.0
    >>> pKa = float(-np.log10(1.8e-5))
    >>> bool(abs(float(titration.pH_at(np.array([Vb_half]))[0]) - pKa) < 1e-3)
    True
    """

    def __init__(self, Ca: float, Va: float, Ka: float, Cb: float, Kw: float = 1.0e-14):
        self.Ca = float(Ca)
        self.Va = float(Va)
        self.Ka = float(Ka)
        self.Cb = float(Cb)
        self.Kw = float(Kw)

    def equivalence_volume(self) -> float:
        """Return the exact (stoichiometric) equivalence volume, ``Ca*Va/Cb``.

        Returns
        -------
        float
        """
        return self.Ca * self.Va / self.Cb

    def _residual(self, h: float, Ca_tot: float, Cb_added: float) -> float:
        A_minus = self.Ka * Ca_tot / (self.Ka + h)
        return Cb_added + h - A_minus - self.Kw / h

    def pH_at(self, Vb):
        Vb = np.atleast_1d(np.asarray(Vb, dtype=np.float64))
        h_values = np.empty_like(Vb)
        for i, vb in enumerate(Vb):
            V_tot = self.Va + vb
            Ca_tot = self.Ca * self.Va / V_tot
            Cb_added = self.Cb * vb / V_tot
            lo, hi = _charge_balance_bracket(Ca_tot, Cb_added, self.Kw)
            h_values[i] = find_positive_root(lambda h, Ca_tot=Ca_tot, Cb_added=Cb_added: self._residual(h, Ca_tot, Cb_added), lo=lo, hi=hi)
        return ph_from_h(h_values)


class WeakBaseStrongAcidTitration(Titration):
    r"""Titrating a weak base B (in the flask) with a strong acid (the titrant).

    The mirror image of :class:`WeakAcidStrongBaseTitration` (swap
    :math:`[H^+] \leftrightarrow [OH^-]`, :math:`K_a \leftrightarrow K_b`).

    Parameters
    ----------
    Cb : float
        Initial weak-base concentration, in mol/L.
    Vb0 : float
        Initial base volume, in L.
    Kb : float
        Base dissociation constant of B.
    Ca : float
        Titrant (strong acid) concentration, in mol/L.
    Kw : float, default 1.0e-14
        Water autoionization constant.

    Examples
    --------
    At the half-equivalence volume, ``pOH = pKb``:

    >>> import numpy as np
    >>> titration = WeakBaseStrongAcidTitration(Cb=0.100, Vb0=0.050, Kb=1.8e-5, Ca=0.100)
    >>> Va_half = titration.equivalence_volume() / 2.0
    >>> pKb = float(-np.log10(1.8e-5))
    >>> pH_half = float(titration.pH_at(np.array([Va_half]))[0])
    >>> pOH_half = 14.0 - pH_half  # Kw = 1e-14 -> pKw = 14
    >>> bool(abs(pOH_half - pKb) < 1e-3)
    True
    """

    def __init__(self, Cb: float, Vb0: float, Kb: float, Ca: float, Kw: float = 1.0e-14):
        self.Cb = float(Cb)
        self.Vb0 = float(Vb0)
        self.Kb = float(Kb)
        self.Ca = float(Ca)
        self.Kw = float(Kw)

    def equivalence_volume(self) -> float:
        """Return the exact (stoichiometric) equivalence volume, ``Cb*Vb0/Ca``.

        Returns
        -------
        float
        """
        return self.Cb * self.Vb0 / self.Ca

    def _residual(self, oh: float, Cb_tot: float, Ca_added: float) -> float:
        BH_plus = self.Kb * Cb_tot / (self.Kb + oh)
        return Ca_added + oh - BH_plus - self.Kw / oh

    def pH_at(self, Va):
        Va = np.atleast_1d(np.asarray(Va, dtype=np.float64))
        oh_values = np.empty_like(Va)
        for i, va in enumerate(Va):
            V_tot = self.Vb0 + va
            Cb_tot = self.Cb * self.Vb0 / V_tot
            Ca_added = self.Ca * va / V_tot
            lo, hi = _charge_balance_bracket(Cb_tot, Ca_added, self.Kw)
            oh_values[i] = find_positive_root(lambda oh, Cb_tot=Cb_tot, Ca_added=Ca_added: self._residual(oh, Cb_tot, Ca_added), lo=lo, hi=hi)
        pOH = poh_from_oh(oh_values)
        return -np.log10(self.Kw) - pOH


def gran_plot(Vb, pH, Va: float) -> GranPlotResult:
    r"""Gran's linearization of a strong-acid/strong-base titration curve.

    Before the equivalence point of a strong acid (initial volume
    :math:`V_a`) titrated with a strong base, the charge balance gives,
    once :math:`[OH^-]` is negligible,

    .. math::

        (V_a + V_b)\,10^{-\mathrm{pH}} = C_aV_a - C_bV_b,

    a straight line in :math:`V_b` whose x-intercept is exactly the
    equivalence volume :math:`V_e = C_aV_a/C_b` -- so :math:`V_e` can be
    found by extrapolating pre-equivalence data, without having to
    resolve the steep inflection itself (G. Gran, *Analyst* 77, 661
    (1952); Harris, *Quantitative Chemical Analysis*, 9th ed., Ch. 11-5).
    Neither :math:`C_a` nor :math:`C_b` needs to be known.

    Parameters
    ----------
    Vb : array-like of float
        Titrant volumes, in L, all *before* the equivalence point
        (conventionally the last 10-90% of the way to it).
    pH : array-like of float
        Measured pH at each volume in `Vb`.
    Va : float
        Initial volume of the acid being titrated, in L.

    Returns
    -------
    GranPlotResult

    Examples
    --------
    Recovering the 50.00 mL equivalence volume of 50.00 mL of 0.100 M HCl
    titrated with 0.100 M NaOH from pre-equivalence points only:

    >>> import numpy as np
    >>> titration = StrongAcidStrongBaseTitration(Ca=0.100, Va=0.050, Cb=0.100)
    >>> Vb = np.linspace(0.030, 0.045, 10)
    >>> result = gran_plot(Vb, titration.pH_at(Vb), Va=0.050)
    >>> round(result.equivalence_volume * 1000.0, 4)
    50.0
    >>> round(-result.slope, 6)  # slope is -Cb
    0.1
    """
    Vb = np.atleast_1d(np.asarray(Vb, dtype=np.float64))
    pH = np.atleast_1d(np.asarray(pH, dtype=np.float64))
    G = (Va + Vb) * 10.0**-pH
    slope, intercept = np.polyfit(Vb, G, 1)
    return GranPlotResult(Vb=Vb, G=G, slope=float(slope), intercept=float(intercept), equivalence_volume=float(-intercept / slope))
