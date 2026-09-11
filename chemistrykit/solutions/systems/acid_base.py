r"""pH/pOH, weak acid/base equilibria, and Henderson-Hasselbalch buffers.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 6.6-6.8, or
Harris, *Quantitative Chemical Analysis*, 9th ed., Ch. 10-11, throughout.

**Exact weak-acid/base equilibrium.** For a weak acid HA at total
(analytical) concentration :math:`C_a`, the full system of equations --
mass balance :math:`[HA]+[A^-]=C_a`, charge balance :math:`[H^+] =
[A^-]+[OH^-]`, :math:`K_a=[H^+][A^-]/[HA]`, :math:`K_w=[H^+][OH^-]` --
reduces (eliminating :math:`[A^-]`, :math:`[HA]`, :math:`[OH^-]` in favor
of :math:`x=[H^+]`) to the cubic

.. math::

    x^3 + K_a x^2 - (K_a C_a + K_w) x - K_a K_w = 0

which has a unique positive real root (:class:`WeakAcid` solves it via
:func:`chemistrykit.solutions.utils.rootfinding.find_positive_real_root`,
rather than the more commonly taught simplification :math:`[H^+] \approx
\sqrt{K_a C_a}` that drops water autoionization -- valid when
:math:`K_a C_a \gg K_w`, but not otherwise). :class:`WeakBase` is the
exact mirror image with :math:`y=[OH^-]`, :math:`K_b`, :math:`C_b`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.solutions.core.base_system import WeakElectrolyte
from chemistrykit.solutions.utils.rootfinding import find_positive_real_root

__all__ = [
    "ph_from_h",
    "h_from_ph",
    "poh_from_oh",
    "oh_from_poh",
    "WeakAcid",
    "WeakBase",
    "Buffer",
    "henderson_hasselbalch_ph",
]


def ph_from_h(h_conc) -> float:
    r"""pH :math:`= -\log_{10}[H^+]`.

    Parameters
    ----------
    h_conc : float or array-like of float
        Hydrogen ion concentration, in mol/L.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> round(float(ph_from_h(1.0e-7)), 6)
    7.0
    """
    h_conc = np.asarray(h_conc, dtype=np.float64)
    return -np.log10(h_conc)


def h_from_ph(pH):
    r"""Inverse of :func:`ph_from_h`: :math:`[H^+] = 10^{-pH}`.

    Examples
    --------
    >>> round(float(h_from_ph(7.0)), 10)
    1e-07
    """
    pH = np.asarray(pH, dtype=np.float64)
    return 10.0**-pH


def poh_from_oh(oh_conc):
    r"""pOH :math:`= -\log_{10}[OH^-]`.

    Examples
    --------
    >>> round(float(poh_from_oh(1.0e-7)), 6)
    7.0
    """
    oh_conc = np.asarray(oh_conc, dtype=np.float64)
    return -np.log10(oh_conc)


def oh_from_poh(pOH):
    r"""Inverse of :func:`poh_from_oh`: :math:`[OH^-] = 10^{-pOH}`."""
    pOH = np.asarray(pOH, dtype=np.float64)
    return 10.0**-pOH


@dataclass
class WeakAcid(WeakElectrolyte):
    r"""A weak monoprotic acid HA in water, solved exactly (see module docstring).

    Parameters
    ----------
    Ca : float
        Total (analytical) acid concentration, in mol/L.
    Ka : float
        Acid dissociation constant.
    Kw : float, default 1.0e-14
        Water autoionization constant (default: 25 degC).

    Examples
    --------
    0.100 M acetic acid (Ka = 1.8e-5): the exact solution is very close to
    the textbook approximation :math:`[H^+] \approx \sqrt{K_a C_a}` here,
    since :math:`K_aC_a \gg K_w`:

    >>> acid = WeakAcid(Ca=0.100, Ka=1.8e-5)
    >>> round(acid.h_concentration(), 6)
    0.001333
    >>> round(float(np.sqrt(1.8e-5 * 0.100)), 6)
    0.001342
    >>> round(acid.pH(), 3)
    2.875
    """

    Ca: float
    Ka: float
    Kw: float = 1.0e-14

    @property
    def total_concentration(self) -> float:
        return self.Ca

    @property
    def equilibrium_constant(self) -> float:
        return self.Ka

    def _ion_concentration(self) -> float:
        coeffs = [1.0, self.Ka, -(self.Ka * self.Ca + self.Kw), -self.Ka * self.Kw]
        return find_positive_real_root(coeffs)

    def h_concentration(self) -> float:
        """Return the equilibrium :math:`[H^+]`, in mol/L.

        Returns
        -------
        float
        """
        return self._ion_concentration()

    def pH(self) -> float:
        """Return the equilibrium pH.

        Returns
        -------
        float
        """
        return float(ph_from_h(self.h_concentration()))


@dataclass
class WeakBase(WeakElectrolyte):
    r"""A weak base B in water, solved exactly (see module docstring).

    Parameters
    ----------
    Cb : float
        Total (analytical) base concentration, in mol/L.
    Kb : float
        Base dissociation constant.
    Kw : float, default 1.0e-14
        Water autoionization constant (default: 25 degC).

    Examples
    --------
    0.100 M ammonia (Kb = 1.8e-5), by the acid/base symmetry the same
    numbers as the :class:`WeakAcid` example, but for pOH:

    >>> base = WeakBase(Cb=0.100, Kb=1.8e-5)
    >>> round(base.oh_concentration(), 6)
    0.001333
    >>> round(base.pH(), 3)
    11.125
    """

    Cb: float
    Kb: float
    Kw: float = 1.0e-14

    @property
    def total_concentration(self) -> float:
        return self.Cb

    @property
    def equilibrium_constant(self) -> float:
        return self.Kb

    def _ion_concentration(self) -> float:
        coeffs = [1.0, self.Kb, -(self.Kb * self.Cb + self.Kw), -self.Kb * self.Kw]
        return find_positive_real_root(coeffs)

    def oh_concentration(self) -> float:
        """Return the equilibrium :math:`[OH^-]`, in mol/L.

        Returns
        -------
        float
        """
        return self._ion_concentration()

    def pOH(self) -> float:
        """Return the equilibrium pOH.

        Returns
        -------
        float
        """
        return float(poh_from_oh(self.oh_concentration()))

    def pH(self) -> float:
        r"""Return the equilibrium pH, via :math:`[H^+] = K_w/[OH^-]`.

        Returns
        -------
        float
        """
        h = self.Kw / self.oh_concentration()
        return float(ph_from_h(h))


def henderson_hasselbalch_ph(pKa: float, base_conc: float, acid_conc: float) -> float:
    r"""The Henderson-Hasselbalch equation: :math:`pH = pK_a + \log_{10}([A^-]/[HA])`.

    A buffer approximation valid when both `base_conc` and `acid_conc` are
    large enough that dissociation/hydrolysis negligibly perturbs their
    equilibrium ratio (L. J. Henderson, *Am. J. Physiol.* 21, 173 (1908);
    K. A. Hasselbalch, *Biochem. Z.* 78, 112 (1917); see Atkins & de
    Paula, *Physical Chemistry*, 11th ed., Ch. 6.7).

    Parameters
    ----------
    pKa : float
        :math:`-\log_{10}(K_a)` of the weak acid.
    base_conc : float
        Concentration (or moles) of the conjugate base, :math:`[A^-]`.
    acid_conc : float
        Concentration (or moles) of the weak acid, :math:`[HA]`.

    Returns
    -------
    float

    Examples
    --------
    An equimolar acid/conjugate-base buffer has pH = pKa exactly:

    >>> round(henderson_hasselbalch_ph(pKa=4.76, base_conc=0.10, acid_conc=0.10), 6)
    4.76
    """
    return float(pKa + np.log10(base_conc / acid_conc))


@dataclass
class Buffer:
    r"""An acid/conjugate-base buffer, via the Henderson-Hasselbalch approximation.

    Parameters
    ----------
    pKa : float
        :math:`-\log_{10}(K_a)` of the weak-acid component.
    acid_conc : float
        Concentration of the weak-acid form, :math:`[HA]`, in mol/L.
    base_conc : float
        Concentration of the conjugate-base form, :math:`[A^-]`, in mol/L.
    """

    pKa: float
    acid_conc: float
    base_conc: float

    def pH(self) -> float:
        """Return the buffer's pH via :func:`henderson_hasselbalch_ph`.

        Returns
        -------
        float
        """
        return henderson_hasselbalch_ph(self.pKa, self.base_conc, self.acid_conc)

    @classmethod
    def from_target_ph(cls, pKa: float, target_pH: float, total_conc: float) -> Buffer:
        r"""Build a :class:`Buffer` with a target pH from a fixed total buffer concentration.

        Inverting the Henderson-Hasselbalch equation for the
        base/acid ratio :math:`r=[A^-]/[HA] = 10^{pH-pKa}` and
        distributing the fixed total concentration
        :math:`C=[HA]+[A^-]` between the two forms in that ratio.

        Parameters
        ----------
        pKa : float
            :math:`-\log_{10}(K_a)` of the weak-acid component.
        target_pH : float
            Desired buffer pH.
        total_conc : float
            Total buffer concentration, :math:`[HA]+[A^-]`, in mol/L.

        Returns
        -------
        Buffer

        Examples
        --------
        Recovering the target pH exactly:

        >>> buf = Buffer.from_target_ph(pKa=4.76, target_pH=5.0, total_conc=0.20)
        >>> round(buf.pH(), 6)
        5.0
        >>> round(buf.acid_conc + buf.base_conc, 6)
        0.2
        """
        ratio = 10.0 ** (target_pH - pKa)
        acid_conc = total_conc / (1.0 + ratio)
        base_conc = total_conc - acid_conc
        return cls(pKa=pKa, acid_conc=acid_conc, base_conc=base_conc)
