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

from chemistrykit.solutions.core.base_system import Titration
from chemistrykit.solutions.systems.acid_base import ph_from_h, poh_from_oh
from chemistrykit.solutions.utils.rootfinding import find_positive_root

__all__ = [
    "StrongAcidStrongBaseTitration",
    "WeakAcidStrongBaseTitration",
    "WeakBaseStrongAcidTitration",
]


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
            h_values[i] = find_positive_root(lambda h, Ca_tot=Ca_tot, Cb_added=Cb_added: self._residual(h, Ca_tot, Cb_added))
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
            oh_values[i] = find_positive_root(lambda oh, Cb_tot=Cb_tot, Ca_added=Ca_added: self._residual(oh, Cb_tot, Ca_added))
        pOH = poh_from_oh(oh_values)
        return -np.log10(self.Kw) - pOH
