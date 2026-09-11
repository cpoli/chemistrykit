"""Shared nonlinear-equilibrium root-finding for solution chemistry.

Both a lone weak acid/base equilibrium
(:mod:`chemistrykit.solutions.systems.acid_base`) and a titration curve
(:mod:`chemistrykit.solutions.systems.titration`) reduce to the same
underlying numerical problem: given a charge-balance equation relating
:math:`[H^+]` (or :math:`[OH^-]`) to known total concentrations and
equilibrium constants, find the unique physically admissible
(positive, real) root -- this is that one shared numerical routine,
following the general house style of keeping numerics that support more
than one ``systems/`` model in ``utils/`` rather than duplicating it (cf.
``chemistrykit.kinetics.utils.linear_regression``,
``chemistrykit.thermo.utils.cubic_roots``).

Two flavors are provided:

* :func:`find_positive_real_root` -- for the *exact* charge-balance
  equation of a lone weak acid or base, which reduces to a cubic
  polynomial in the ion concentration (see
  :mod:`chemistrykit.solutions.systems.acid_base`'s module docstring for
  the derivation); solved via :func:`numpy.roots`.
* :func:`find_positive_root` -- for a general monotonic charge-balance
  *function* (not necessarily reducible to a low-degree polynomial by
  hand, e.g. mid-titration with an added strong-electrolyte term), solved
  by bracketing a sign change and refining with :func:`scipy.optimize.brentq`.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

__all__ = ["find_positive_real_root", "find_positive_root"]


def find_positive_real_root(coeffs, tol: float = 1e-9) -> float:
    r"""Return the unique positive real root of a polynomial.

    Parameters
    ----------
    coeffs : array-like of float
        Polynomial coefficients, highest degree first (as for
        :func:`numpy.roots`).
    tol : float, default 1e-9
        A root is treated as real if ``abs(imag part) <= tol *
        max(abs(real part), 1.0)``.

    Returns
    -------
    float

    Raises
    ------
    ValueError
        If no positive real root exists.

    Examples
    --------
    >>> round(find_positive_real_root([1.0, 1.0, -2.0]), 6)  # x^2+x-2=(x+2)(x-1)
    1.0
    """
    roots = np.roots(np.asarray(coeffs, dtype=np.float64))
    candidates = [r.real for r in roots if abs(r.imag) <= tol * max(abs(r.real), 1.0) and r.real > 0.0]
    if not candidates:
        raise ValueError("no positive real root found")
    # A charge-balance cubic for a physically sensible (positive Ka/Kb,
    # Ca/Cb) weak electrolyte has exactly one positive real root; if
    # floating-point noise ever produces more than one, the largest is
    # the physically dominant one (the others are spurious near-zero
    # artifacts of the water-autoionization term).
    return float(max(candidates))


def find_positive_root(residual, lo: float = 1e-14, hi: float = 1.0, n_scan: int = 400) -> float:
    r"""Find the unique positive root of a monotonic scalar function via bracketing + Brent's method.

    Scans `n_scan` log-spaced points in ``[lo, hi]`` for a sign change in
    `residual`, then refines the bracket with
    :func:`scipy.optimize.brentq`. Used where the charge-balance equation
    is not conveniently reducible to a polynomial by hand (e.g. a
    titration curve, where the total concentration itself depends on the
    volume of titrant added).

    Parameters
    ----------
    residual : callable
        ``residual(x) -> float``, monotonic in `x` on ``[lo, hi]``, with
        a single root.
    lo, hi : float, default 1e-14, 1.0
        Search bracket, in the same units as the equilibrium
        concentration being solved for (mol/L for :math:`[H^+]` or
        :math:`[OH^-]`, which is why the default bracket spans a physically
        sensible pH range of about 0-14).
    n_scan : int, default 400
        Number of log-spaced points used to search for a sign change.

    Returns
    -------
    float

    Raises
    ------
    ValueError
        If no sign change is found in ``[lo, hi]``.

    Examples
    --------
    >>> round(find_positive_root(lambda x: x**2 - 4.0, lo=1e-3, hi=10.0), 6)
    2.0
    """
    grid = np.logspace(np.log10(lo), np.log10(hi), n_scan)
    values = np.array([residual(x) for x in grid])
    sign_changes = np.where(np.diff(np.sign(values)) != 0)[0]
    if sign_changes.size == 0:
        raise ValueError("no sign change found in [lo, hi]; residual may not have a root in this bracket")
    i = sign_changes[0]
    return float(brentq(residual, grid[i], grid[i + 1]))
