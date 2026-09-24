r"""Savitzky-Golay smoothing and differentiation of evenly sampled analytical signals.

Implemented from first principles as a moving local least-squares
polynomial fit (A. Savitzky and M. J. E. Golay, *Anal. Chem.* 36, 1627
(1964)): fitting a degree-`p` polynomial to each window of `2m+1` points
and evaluating it (or its derivative) at the window's center is
equivalent to convolving the data with a fixed set of weights, which
depend only on `m`, `p`, and the derivative order.
"""

from __future__ import annotations

from math import factorial

import numpy as np

__all__ = ["savitzky_golay_coefficients", "savitzky_golay"]


def savitzky_golay_coefficients(window_length: int, polyorder: int, deriv: int = 0, delta: float = 1.0) -> np.ndarray:
    r"""Convolution weights of a Savitzky-Golay filter.

    With window offsets :math:`z=-m,\ldots,m` and the design matrix
    :math:`J_{ij}=z_i^{\,j}` (:math:`j=0,\ldots,p`), the least-squares
    polynomial coefficients are :math:`(J^TJ)^{-1}J^T y`; row `d` of that
    matrix, times :math:`d!/\Delta^d`, gives the weights for the `d`-th
    derivative at the window center.

    Parameters
    ----------
    window_length : int
        Odd number of points in each window, :math:`2m+1`.
    polyorder : int
        Degree `p` of the local polynomial, less than `window_length`.
    deriv : int, default 0
        Derivative order (0 = smoothing).
    delta : float, default 1.0
        Sample spacing, used to scale derivatives.

    Returns
    -------
    ndarray, shape (window_length,)
        Weights :math:`c_i` such that the filtered value at point `k` is
        :math:`\sum_i c_i\,y_{k+i}`, `i` running from `-m` to `m`.

    Examples
    --------
    The classic 5-point quadratic smoothing weights from Savitzky and
    Golay's Table I, :math:`(-3, 12, 17, 12, -3)/35`:

    >>> import numpy as np
    >>> np.round(savitzky_golay_coefficients(5, 2) * 35, 10)
    array([-3., 12., 17., 12., -3.])
    """
    if window_length % 2 != 1 or window_length < 1:
        raise ValueError("window_length must be a positive odd integer")
    if not 0 <= polyorder < window_length:
        raise ValueError("polyorder must satisfy 0 <= polyorder < window_length")
    if not 0 <= deriv <= polyorder:
        raise ValueError("deriv must satisfy 0 <= deriv <= polyorder")
    m = window_length // 2
    z = np.arange(-m, m + 1, dtype=np.float64)
    J = np.vander(z, polyorder + 1, increasing=True)
    return np.linalg.pinv(J)[deriv] * factorial(deriv) / delta**deriv


def savitzky_golay(y, window_length: int, polyorder: int, deriv: int = 0, delta: float = 1.0) -> np.ndarray:
    r"""Smooth (or differentiate) a signal with a Savitzky-Golay filter.

    Interior points use the fixed convolution weights of
    :func:`savitzky_golay_coefficients`; the first and last `m` points,
    which have no full centered window, are taken from the polynomial
    fitted to the first and last full window (the same edge treatment as
    ``scipy.signal.savgol_filter(..., mode="interp")``).

    Parameters
    ----------
    y : array-like of float
        Evenly sampled signal, at least `window_length` points long.
    window_length, polyorder, deriv, delta
        As in :func:`savitzky_golay_coefficients`.

    Returns
    -------
    ndarray, same shape as `y`

    Examples
    --------
    A filter of polynomial order `p` passes any polynomial of degree <= `p`
    through unchanged, edges included:

    >>> import numpy as np
    >>> x = np.linspace(0.0, 1.0, 21)
    >>> y = 3.0 * x**2 - x + 0.5
    >>> bool(np.allclose(savitzky_golay(y, 7, 2), y))
    True
    """
    y = np.asarray(y, dtype=np.float64)
    n = len(y)
    if n < window_length:
        raise ValueError("signal must be at least window_length points long")
    m = window_length // 2
    c = savitzky_golay_coefficients(window_length, polyorder, deriv, delta)
    out = np.empty_like(y)
    out[m : n - m] = np.correlate(y, c, mode="valid")
    z = np.arange(window_length, dtype=np.float64)
    for sl, positions in ((slice(0, window_length), np.arange(m)), (slice(n - window_length, n), np.arange(window_length - m, window_length))):
        poly = np.polynomial.Polynomial.fit(z, y[sl], polyorder, domain=[0.0, 1.0], window=[0.0, 1.0]).deriv(deriv)
        vals = poly(positions.astype(np.float64)) / delta**deriv
        if sl.start == 0:
            out[:m] = vals
        else:
            out[n - m :] = vals
    return out
