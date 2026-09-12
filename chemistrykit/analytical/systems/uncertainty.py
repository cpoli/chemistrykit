r"""Propagation of uncertainty: closed-form rules for sums/products/powers, and a general numerical form.

First-order (linearized) error propagation: for a function :math:`y=f(x_1,
\ldots,x_n)` of independent (uncorrelated) measured quantities each with
random uncertainty :math:`\sigma_{x_i}`, the propagated uncertainty in
`y` is

.. math::

    \sigma_y=\sqrt{\sum_i\left(\frac{\partial f}{\partial x_i}\right)^2\sigma_{x_i}^2}

(Harris, *Quantitative Chemical Analysis*, 9th ed., Ch. 3, or Taylor, *An
Introduction to Error Analysis*, 2nd ed. (1997), Ch. 3). The three
closed-form special cases below (sum/difference, product/quotient,
power) are the textbook shortcuts of this general formula for the
specific operations chemists compute by hand most often;
:func:`propagate_uncertainty` implements the general formula directly
(via numerical partial derivatives) for an arbitrary function.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "propagate_sum",
    "propagate_product",
    "propagate_power",
    "propagate_uncertainty",
]


def propagate_sum(sigmas) -> float:
    r"""Absolute uncertainty of a sum or difference of independent quantities, :math:`\sigma=\sqrt{\sum_i\sigma_i^2}`.

    Applies identically whether the terms are added or subtracted (Harris,
    *Quantitative Chemical Analysis*, 9th ed., Ch. 3.2) -- only the
    *absolute* uncertainties of each term matter, not their signs.

    Parameters
    ----------
    sigmas : array-like of float
        Absolute uncertainty of each term.

    Returns
    -------
    float

    Examples
    --------
    Two equal uncertainties combine to :math:`\sigma\sqrt2`, not :math:`2\sigma`:

    >>> round(propagate_sum([0.02, 0.02]), 6)
    0.028284
    """
    sigmas = np.asarray(sigmas, dtype=np.float64)
    return float(np.sqrt(np.sum(sigmas**2)))


def propagate_product(values, sigmas) -> float:
    r"""Absolute uncertainty of a product or quotient of independent quantities.

    Relative uncertainties add in quadrature (Harris, *Quantitative
    Chemical Analysis*, 9th ed., Ch. 3.2):

    .. math::

        \frac{\sigma_y}{|y|}=\sqrt{\sum_i\left(\frac{\sigma_i}{x_i}\right)^2},
        \qquad y=\prod_i x_i^{\pm1}

    Parameters
    ----------
    values : array-like of float
        The value of each factor (as it appears in the product/quotient,
        e.g. use ``1/x`` for a quantity in the denominator -- see the
        Examples).
    sigmas : array-like of float
        Absolute uncertainty of each factor, same order as `values`.

    Returns
    -------
    float
        Absolute uncertainty of the product :math:`y=\prod_i` `values[i]`.

    Examples
    --------
    :math:`y=a\times b` with :math:`a=2.0\pm0.1`, :math:`b=3.0\pm0.2`:

    >>> y = 2.0 * 3.0
    >>> sigma_y = propagate_product([2.0, 3.0], [0.1, 0.2])
    >>> round(sigma_y, 6)
    0.5

    A quotient :math:`y=a/b` is a product with `b` replaced by :math:`1/b`
    (relative uncertainty of :math:`1/b` equals that of `b`):

    >>> a, sigma_a = 10.0, 0.5
    >>> b, sigma_b = 4.0, 0.2
    >>> y = a / b
    >>> sigma_y = propagate_product([a, 1.0 / b], [sigma_a, sigma_b / b**2])
    >>> round(sigma_y, 6)
    0.176777
    """
    values = np.asarray(values, dtype=np.float64)
    sigmas = np.asarray(sigmas, dtype=np.float64)
    y = np.prod(values)
    relative_variance = np.sum((sigmas / values) ** 2)
    return float(abs(y) * np.sqrt(relative_variance))


def propagate_power(value: float, sigma: float, exponent: float) -> float:
    r"""Absolute uncertainty of :math:`y=x^n`, :math:`\sigma_y=|n|\,|x|^{n-1}\sigma_x`.

    Equivalently :math:`\sigma_y/|y|=|n|\,\sigma_x/|x|` (Harris,
    *Quantitative Chemical Analysis*, 9th ed., Ch. 3.2).

    Parameters
    ----------
    value : float
        The base `x`.
    sigma : float
        Absolute uncertainty of `x`.
    exponent : float
        The exponent `n` (need not be an integer).

    Returns
    -------
    float

    Examples
    --------
    :math:`y=x^2` with :math:`x=3.0\pm0.1`: relative uncertainty doubles:

    >>> round(propagate_power(3.0, 0.1, 2.0), 6)
    0.6
    """
    return float(abs(exponent) * abs(value) ** (exponent - 1.0) * sigma)


def propagate_uncertainty(func, values, sigmas, h: float = 1e-6) -> float:
    r"""General first-order uncertainty propagation for an arbitrary differentiable function.

    Evaluates :math:`\sigma_y=\sqrt{\sum_i(\partial f/\partial
    x_i)^2\sigma_i^2}` (Harris, *Quantitative Chemical Analysis*, 9th ed.,
    Ch. 3.2, eq. 3.2) with each partial derivative estimated by a
    central finite difference -- useful for a function with no simple
    closed-form propagation rule (e.g. a pH from a Nernst equation, or a
    rate constant from an Arrhenius fit), at the cost of losing the exact
    closed-form results :func:`propagate_sum`/:func:`propagate_product`/
    :func:`propagate_power` give for their specific operations.

    Parameters
    ----------
    func : callable
        ``func(*values) -> float``.
    values : array-like of float
        The value of each independent variable, in the order `func` expects.
    sigmas : array-like of float
        Absolute uncertainty of each variable, same order as `values`.
    h : float, default 1e-6
        Relative step size for the central finite difference.

    Returns
    -------
    float

    Examples
    --------
    Reproduces :func:`propagate_product` for a simple product, since both
    implement the same general formula:

    >>> from chemistrykit.analytical.systems.uncertainty import propagate_product
    >>> sigma_general = propagate_uncertainty(lambda a, b: a * b, [2.0, 3.0], [0.1, 0.2])
    >>> sigma_closed_form = propagate_product([2.0, 3.0], [0.1, 0.2])
    >>> bool(abs(sigma_general - sigma_closed_form) < 1e-6)
    True
    """
    values = np.asarray(values, dtype=np.float64)
    sigmas = np.asarray(sigmas, dtype=np.float64)
    variance = 0.0
    for i in range(len(values)):
        step = h * max(abs(values[i]), 1.0)
        forward = values.copy()
        backward = values.copy()
        forward[i] += step
        backward[i] -= step
        partial = (func(*forward) - func(*backward)) / (2.0 * step)
        variance += (partial * sigmas[i]) ** 2
    return float(np.sqrt(variance))
