"""A generic real-positive-root finder, shared by the cubic EOS models.

Both :class:`~chemistrykit.thermo.systems.equations_of_state.VanDerWaals`
and :class:`~chemistrykit.thermo.systems.equations_of_state.RedlichKwong`
need to solve a cubic polynomial in molar volume :math:`V_m` and then
discard the physically inadmissible roots (complex, or non-positive) --
this is that one shared numerical routine, following the general house
style of keeping numerics that support more than one ``systems/`` model
in ``utils/`` rather than duplicating it (cf.
``chemistrykit.kinetics.utils.linear_regression``).
"""

from __future__ import annotations

import numpy as np

__all__ = ["real_positive_roots"]


def real_positive_roots(coeffs, tol: float = 1e-7) -> np.ndarray:
    r"""Return the positive real roots of a polynomial, sorted ascending.

    A thin wrapper around :func:`numpy.roots` that filters out genuinely
    complex roots (imaginary part large relative to the real part) and
    non-positive roots. A cubic equation of state below its critical
    temperature can have up to three real positive roots (a metastable
    liquid branch, an unstable middle root, and a vapor branch); this
    function returns all of them so the caller can select a branch.

    Parameters
    ----------
    coeffs : array-like of float
        Polynomial coefficients, highest degree first (as for
        :func:`numpy.roots`).
    tol : float, default 1e-7
        A root is treated as real if ``abs(imag part) <= tol *
        max(abs(real part), 1.0)``.

    Returns
    -------
    ndarray
        Positive real roots, sorted ascending. Empty if none exist.

    Examples
    --------
    >>> import numpy as np
    >>> roots = real_positive_roots([1.0, -6.0, 11.0, -6.0])  # (x-1)(x-2)(x-3)
    >>> np.allclose(roots, [1.0, 2.0, 3.0])
    True

    A root at zero or a negative root is discarded:

    >>> real_positive_roots([1.0, 0.0, -1.0])  # (x-1)(x+1) -> roots +-1
    array([1.])
    """
    roots = np.roots(np.asarray(coeffs, dtype=np.float64))
    kept = [r.real for r in roots if abs(r.imag) <= tol * max(abs(r.real), 1.0) and r.real > 0.0]
    return np.sort(np.array(kept, dtype=np.float64))
