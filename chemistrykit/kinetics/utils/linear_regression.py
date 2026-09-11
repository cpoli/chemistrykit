"""Shared ordinary-least-squares linear regression.

Both :func:`chemistrykit.kinetics.systems.arrhenius.fit_arrhenius` and
:func:`chemistrykit.kinetics.systems.enzyme.fit_lineweaver_burk` reduce to
fitting a straight line to a linearized form of their respective rate law
(``ln k`` vs. ``1/T``, and ``1/v`` vs. ``1/[S]``, respectively) and
reporting a goodness-of-fit statistic -- this is that one shared
numerical routine, following the general house style of keeping numerics
that support more than one ``systems/`` model in ``utils/`` rather than
duplicating it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["LinearFit", "linear_fit"]


@dataclass
class LinearFit:
    """Result of an ordinary-least-squares straight-line fit ``y = slope*x + intercept``."""

    slope: float
    intercept: float
    r_squared: float


def linear_fit(x, y) -> LinearFit:
    """Fit ``y = slope*x + intercept`` by ordinary least squares.

    Parameters
    ----------
    x, y : array-like of float
        Data points (at least 2, with `x` not all identical).

    Returns
    -------
    LinearFit

    Examples
    --------
    >>> fit = linear_fit([0.0, 1.0, 2.0, 3.0], [1.0, 3.0, 5.0, 7.0])
    >>> round(fit.slope, 6), round(fit.intercept, 6), round(fit.r_squared, 6)
    (2.0, 1.0, 1.0)
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return LinearFit(slope=float(slope), intercept=float(intercept), r_squared=float(r_squared))
