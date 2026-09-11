"""Shared ordinary-least-squares linear regression.

:func:`chemistrykit.thermo.systems.equilibrium.fit_van_t_hoff` fits a
straight line to the linearized van't Hoff equation (``ln K`` vs. ``1/T``)
and reports a goodness-of-fit statistic -- the same small numerical
routine as :mod:`chemistrykit.kinetics.utils.linear_regression` (used
there for the analogous Arrhenius-plot fit), reimplemented here rather
than imported cross-domain so that each chemistrykit domain subpackage
depends only on the shared :mod:`chemistrykit.constants` /
:mod:`chemistrykit.integrators` infrastructure, not on its sibling
domains.
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
