"""Shared ordinary-least-squares linear regression, with the residual standard error a calibration curve needs.

Reimplemented here rather than imported cross-domain, following
:mod:`chemistrykit.surface.utils.regression`/:mod:`chemistrykit.electrochem.utils.regression`/
:mod:`chemistrykit.photochem.utils.regression`'s precedent of each domain
keeping its own small linear-fit utility rather than depending on a
sibling domain for it. Extended, relative to those, with the residual
standard error :math:`s_{y/x}` that
:mod:`chemistrykit.analytical.systems.calibration`'s LOD/LOQ calculation
needs.
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
    residual_std_error: float
    """float: The residual standard deviation about the regression line,
    :math:`s_{y/x}=\\sqrt{\\sum(y_i-\\hat y_i)^2/(n-2)}` (Harris,
    *Quantitative Chemical Analysis*, 9th ed., Ch. 4.5), the quantity a
    calibration curve's LOD/LOQ are built from."""


def linear_fit(x, y) -> LinearFit:
    """Fit ``y = slope*x + intercept`` by ordinary least squares.

    Parameters
    ----------
    x, y : array-like of float
        Data points (at least 3, so `residual_std_error` is defined, with
        `x` not all identical).

    Returns
    -------
    LinearFit

    Examples
    --------
    >>> fit = linear_fit([0.0, 1.0, 2.0, 3.0], [1.0, 3.0, 5.0, 7.0])
    >>> round(fit.slope, 6), round(fit.intercept, 6), round(fit.r_squared, 6)
    (2.0, 1.0, 1.0)
    >>> fit.residual_std_error < 1e-9
    True
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = len(x)
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    residuals = y - y_pred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    s_yx = np.sqrt(ss_res / (n - 2)) if n > 2 else 0.0
    return LinearFit(slope=float(slope), intercept=float(intercept), r_squared=float(r_squared), residual_std_error=float(s_yx))
