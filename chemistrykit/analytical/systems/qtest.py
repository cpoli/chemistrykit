r"""Dixon's Q-test for rejecting a single suspect outlier from a small data set.

See Harris, *Quantitative Chemical Analysis*, 9th ed., Ch. 4.6 and its
Table 4.5 (reproduced below), and D. B. Rorabacher, *Anal. Chem.* 63, 139
(1991) for the corrected/extended critical-value table this module uses
(Rorabacher's values supersede the original, less precise Dean & Dixon
(1951) table still seen in some older textbooks).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["Q_CRITICAL_TABLE", "QTestResult", "dixon_q_test"]

#: dict: Dixon's Q-test critical values, :math:`Q_{crit}(n,\text{confidence})`,
#: for sample size `n` (keys 3-10) at 90%, 95%, and 99% confidence
#: (Rorabacher, *Anal. Chem.* 63, 139 (1991), Table II).
Q_CRITICAL_TABLE: dict = {
    3: {0.90: 0.941, 0.95: 0.970, 0.99: 0.994},
    4: {0.90: 0.765, 0.95: 0.829, 0.99: 0.926},
    5: {0.90: 0.642, 0.95: 0.710, 0.99: 0.821},
    6: {0.90: 0.560, 0.95: 0.625, 0.99: 0.740},
    7: {0.90: 0.507, 0.95: 0.568, 0.99: 0.680},
    8: {0.90: 0.468, 0.95: 0.526, 0.99: 0.634},
    9: {0.90: 0.437, 0.95: 0.493, 0.99: 0.598},
    10: {0.90: 0.412, 0.95: 0.466, 0.99: 0.568},
}


@dataclass
class QTestResult:
    """Result of a :func:`dixon_q_test` call."""

    suspect_value: float
    """float: The data point tested (the minimum or maximum of the set)."""

    Q_statistic: float
    """float: The computed Q statistic, gap/range."""

    Q_critical: float
    """float: The critical value looked up from :data:`Q_CRITICAL_TABLE`."""

    reject: bool
    """bool: Whether the suspect value is rejected as an outlier (``Q_statistic > Q_critical``)."""


def dixon_q_test(data, confidence: float = 0.95, suspect: str = "auto") -> QTestResult:
    r"""Test whether the smallest or largest value in a small data set is an outlier.

    .. math::

        Q = \frac{\text{gap}}{\text{range}}

    where "gap" is the absolute difference between the suspect value and
    its nearest neighbor, and "range" is the difference between the
    largest and smallest values in the set (Harris, *Quantitative
    Chemical Analysis*, 9th ed., Ch. 4.6). The suspect value is rejected
    as an outlier if :math:`Q>Q_{crit}` for the sample size and
    confidence level chosen.

    Parameters
    ----------
    data : array-like of float
        The data set, 3 to 10 points.
    confidence : {0.90, 0.95, 0.99}, default 0.95
        Confidence level for the critical value.
    suspect : {"auto", "min", "max"}, default "auto"
        Which point to test. ``"auto"`` tests whichever of the minimum or
        maximum is farther (in absolute value) from the sample mean.

    Returns
    -------
    QTestResult

    Raises
    ------
    ValueError
        If `data` has fewer than 3 or more than 10 points, or
        `confidence` is not one of the tabulated levels.

    Examples
    --------
    A clear outlier in a small data set is correctly rejected at 95%
    confidence:

    >>> result = dixon_q_test([10.0, 10.1, 10.2, 10.3, 15.0], confidence=0.95)
    >>> result.suspect_value
    15.0
    >>> result.reject
    True

    A borderline value that is *not* extreme enough is retained:

    >>> result = dixon_q_test([10.0, 10.1, 10.2, 10.3, 10.9], confidence=0.99)
    >>> result.reject
    False
    """
    if not (3 <= len(data) <= 10):
        raise ValueError("Dixon's Q-test is tabulated here for 3-10 data points")
    if confidence not in Q_CRITICAL_TABLE[3]:
        raise ValueError(f"confidence must be one of {sorted(Q_CRITICAL_TABLE[3])}")

    data = np.sort(np.asarray(data, dtype=np.float64))
    n = len(data)
    data_range = data[-1] - data[0]

    gap_low = data[1] - data[0]
    gap_high = data[-1] - data[-2]

    if suspect == "min":
        suspect_value, gap = data[0], gap_low
    elif suspect == "max":
        suspect_value, gap = data[-1], gap_high
    elif suspect == "auto":
        mean = np.mean(data)
        if abs(data[0] - mean) >= abs(data[-1] - mean):
            suspect_value, gap = data[0], gap_low
        else:
            suspect_value, gap = data[-1], gap_high
    else:
        raise ValueError('suspect must be "auto", "min", or "max"')

    Q = gap / data_range if data_range > 0 else 0.0
    Q_crit = Q_CRITICAL_TABLE[n][confidence]
    return QTestResult(suspect_value=float(suspect_value), Q_statistic=float(Q), Q_critical=float(Q_crit), reject=bool(Q > Q_crit))
