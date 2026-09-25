"""Elementary integrated rate laws and half-lives (zero/first/second order).

These are the textbook closed-form solutions of ``d[A]/dt = -k[A]^n`` for
``n = 0, 1, 2`` -- see e.g. Atkins & de Paula, *Physical Chemistry*, 11th
ed., Ch. 20 ("The Rates of Chemical Reactions"), Table 20.3, or Levine,
*Physical Chemistry*, 6th ed., Ch. 17. No numerical integration is used:
each class implements the exact analytic solution directly.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.kinetics.core.base_system import RateLaw

__all__ = ["ZeroOrder", "FirstOrder", "SecondOrder"]


class ZeroOrder(RateLaw):
    r"""Zero-order kinetics: ``d[A]/dt = -k``.

    Integrated law: :math:`[A](t) = [A]_0 - kt`, valid only while
    :math:`[A](t) \geq 0` -- physically, the reaction stops (or crosses
    over to a different rate law) once the reactant is exhausted, at
    :math:`t = [A]_0 / k`. :meth:`concentration` and :meth:`rate` model
    that stop: past the exhaustion time the concentration stays at 0 and
    the rate drops to 0, rather than extrapolating the line to negative
    concentrations.

    Parameters
    ----------
    k : float
        Zero-order rate constant, in concentration/time.
    C0 : float
        Initial concentration :math:`[A]_0`.

    Examples
    --------
    >>> law = ZeroOrder(k=0.1, C0=1.0)
    >>> round(float(law.concentration(5.0)), 6)
    0.5
    >>> round(law.half_life(), 6)
    5.0
    >>> float(law.concentration(20.0))  # exhausted at t = C0/k = 10
    0.0
    """

    def __init__(self, k: float, C0: float):
        if k <= 0:
            raise ValueError("k must be positive")
        if C0 <= 0:
            raise ValueError("C0 must be positive")
        self.k = float(k)
        self.C0 = float(C0)

    def concentration(self, t):
        r"""Return :math:`\max([A]_0 - kt, 0)` at time(s) `t`.

        Parameters
        ----------
        t : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        t = np.asarray(t, dtype=np.float64)
        return np.maximum(self.C0 - self.k * t, 0.0)

    def rate(self, t=None):
        """Return the rate: ``k`` until the reactant is exhausted, 0 after.

        Parameters
        ----------
        t : float or array-like of float, optional
            Time at which to evaluate the rate; defaults to 0 (the
            initial rate ``k``). The rate does not depend on
            concentration, only on whether any reactant remains, i.e.
            whether ``t < C0 / k``.

        Returns
        -------
        float or ndarray
        """
        if t is None:
            return self.k
        t = np.asarray(t, dtype=np.float64)
        return np.where(t < self.C0 / self.k, self.k, 0.0)

    def half_life(self) -> float:
        r""":math:`t_{1/2} = [A]_0 / (2k)`, per Atkins & de Paula Table 20.3."""
        return self.C0 / (2.0 * self.k)


class FirstOrder(RateLaw):
    r"""First-order kinetics: ``d[A]/dt = -k[A]``.

    Integrated law: :math:`[A](t) = [A]_0 e^{-kt}`. The half-life
    :math:`t_{1/2} = \ln(2)/k` is independent of the initial
    concentration -- the defining signature of first-order kinetics
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 20, Table
    20.3), used e.g. for radioactive decay and many unimolecular
    reactions.

    Parameters
    ----------
    k : float
        First-order rate constant, in 1/time.
    C0 : float
        Initial concentration :math:`[A]_0`.

    Examples
    --------
    >>> import numpy as np
    >>> law = FirstOrder(k=np.log(2), C0=1.0)
    >>> round(float(law.concentration(1.0)), 6)
    0.5
    >>> round(law.half_life(), 6)
    1.0
    """

    def __init__(self, k: float, C0: float):
        if k <= 0:
            raise ValueError("k must be positive")
        if C0 <= 0:
            raise ValueError("C0 must be positive")
        self.k = float(k)
        self.C0 = float(C0)

    def concentration(self, t):
        t = np.asarray(t, dtype=np.float64)
        return self.C0 * np.exp(-self.k * t)

    def rate(self, t=None):
        """Return the instantaneous rate ``k * [A](t)``.

        Parameters
        ----------
        t : float or array-like of float, optional
            Time at which to evaluate the rate; defaults to 0 (the
            initial rate ``k * C0``).

        Returns
        -------
        float or ndarray
        """
        t = 0.0 if t is None else t
        return self.k * self.concentration(t)

    def half_life(self) -> float:
        r""":math:`t_{1/2} = \ln(2)/k`, independent of :math:`[A]_0`."""
        return float(np.log(2.0) / self.k)


class SecondOrder(RateLaw):
    r"""Second-order kinetics (single reactant): ``d[A]/dt = -k[A]^2``.

    Integrated law: :math:`1/[A](t) = 1/[A]_0 + kt`, i.e.
    :math:`[A](t) = [A]_0 / (1 + k[A]_0 t)`. Unlike first order, the
    half-life :math:`t_{1/2} = 1/(k[A]_0)` *does* depend on the initial
    concentration (Atkins & de Paula, *Physical Chemistry*, 11th ed.,
    Ch. 20, Table 20.3).

    Parameters
    ----------
    k : float
        Second-order rate constant, in 1/(concentration*time).
    C0 : float
        Initial concentration :math:`[A]_0`.

    Examples
    --------
    >>> law = SecondOrder(k=1.0, C0=1.0)
    >>> round(float(law.concentration(1.0)), 6)
    0.5
    >>> round(law.half_life(), 6)
    1.0
    """

    def __init__(self, k: float, C0: float):
        if k <= 0:
            raise ValueError("k must be positive")
        if C0 <= 0:
            raise ValueError("C0 must be positive")
        self.k = float(k)
        self.C0 = float(C0)

    def concentration(self, t):
        t = np.asarray(t, dtype=np.float64)
        return self.C0 / (1.0 + self.k * self.C0 * t)

    def rate(self, t=None):
        """Return the instantaneous rate ``k * [A](t)^2``.

        Parameters
        ----------
        t : float or array-like of float, optional
            Time at which to evaluate the rate; defaults to 0 (the
            initial rate ``k * C0^2``).

        Returns
        -------
        float or ndarray
        """
        t = 0.0 if t is None else t
        c = self.concentration(t)
        return self.k * c**2

    def half_life(self) -> float:
        r""":math:`t_{1/2} = 1/(k[A]_0)`, unlike first order, depends on :math:`[A]_0`."""
        return 1.0 / (self.k * self.C0)
