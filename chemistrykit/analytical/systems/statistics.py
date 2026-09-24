r"""Small-sample statistics for replicate measurements: Student's t confidence interval, Grubbs' outlier test, and the Horwitz precision function.

See Harris, *Quantitative Chemical Analysis*, 9th ed., Ch. 4, or Miller &
Miller, *Statistics and Chemometrics for Analytical Chemistry*, 6th ed.,
Ch. 2-3, throughout; each function cites its original source below.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

__all__ = [
    "ConfidenceIntervalResult",
    "t_confidence_interval",
    "GrubbsTestResult",
    "grubbs_critical_value",
    "grubbs_test",
    "horwitz_rsd",
    "horrat",
]


@dataclass
class ConfidenceIntervalResult:
    """Result of a :func:`t_confidence_interval` call."""

    mean: float
    """float: Sample mean of the replicates."""

    std: float
    """float: Sample standard deviation `s` (``ddof=1``)."""

    n: int
    """int: Number of replicates."""

    t_critical: float
    """float: Two-sided Student's `t` critical value for ``n - 1`` degrees of freedom."""

    half_width: float
    r"""float: :math:`t\,s/\sqrt n`, the interval's half-width."""

    @property
    def lower(self) -> float:
        """float: Lower end of the interval."""
        return self.mean - self.half_width

    @property
    def upper(self) -> float:
        """float: Upper end of the interval."""
        return self.mean + self.half_width


def t_confidence_interval(data, confidence: float = 0.95) -> ConfidenceIntervalResult:
    r"""Student's `t` confidence interval for the true mean of a few replicate measurements.

    .. math::

        \mu = \bar x \pm \frac{t\,s}{\sqrt n}

    with `s` the sample standard deviation and `t` the two-sided critical
    value of Student's distribution with :math:`n-1` degrees of freedom
    ("Student" [W. S. Gosset], *Biometrika* 6, 1 (1908)). For small `n`,
    `t` is much larger than the normal-distribution value 1.96, widening
    the interval to account for `s` itself being an uncertain estimate of
    the true standard deviation.

    Parameters
    ----------
    data : array-like of float
        Replicate measurements (at least 2).
    confidence : float, default 0.95
        Two-sided confidence level.

    Returns
    -------
    ConfidenceIntervalResult

    Examples
    --------
    For five replicates, the 95% two-sided `t` value is the tabulated 2.776:

    >>> ci = t_confidence_interval([10.1, 10.3, 9.9, 10.2, 10.0])
    >>> round(ci.t_critical, 3), round(ci.mean, 6)
    (2.776, 10.1)
    >>> round(ci.half_width, 4)
    0.1963
    """
    x = np.asarray(data, dtype=np.float64)
    n = len(x)
    if n < 2:
        raise ValueError("a t confidence interval needs at least 2 measurements")
    mean = float(np.mean(x))
    s = float(np.std(x, ddof=1))
    t_crit = float(stats.t.ppf(0.5 + confidence / 2.0, df=n - 1))
    return ConfidenceIntervalResult(mean=mean, std=s, n=n, t_critical=t_crit, half_width=float(t_crit * s / np.sqrt(n)))


@dataclass
class GrubbsTestResult:
    """Result of a :func:`grubbs_test` call."""

    suspect_value: float
    """float: The data point tested (the one farthest from the mean)."""

    G_statistic: float
    r"""float: :math:`G=|x_{suspect}-\bar x|/s`."""

    G_critical: float
    """float: The critical value from :func:`grubbs_critical_value`."""

    reject: bool
    """bool: Whether the suspect value is rejected as an outlier (``G_statistic > G_critical``)."""


def grubbs_critical_value(n: int, alpha: float = 0.05) -> float:
    r"""Two-sided critical value of Grubbs' statistic for a sample of size `n`.

    .. math::

        G_{crit} = \frac{n-1}{\sqrt n}\sqrt{\frac{t^2}{n-2+t^2}}

    with `t` the upper :math:`\alpha/(2n)` critical value of Student's
    distribution with :math:`n-2` degrees of freedom (F. E. Grubbs,
    *Technometrics* 11, 1 (1969)).

    Parameters
    ----------
    n : int
        Sample size (at least 3).
    alpha : float, default 0.05
        Significance level.

    Returns
    -------
    float

    Examples
    --------
    Matches the tabulated two-sided 5% value for ten observations:

    >>> round(grubbs_critical_value(10, alpha=0.05), 3)
    2.29
    """
    if n < 3:
        raise ValueError("Grubbs' test needs at least 3 observations")
    t = stats.t.ppf(1.0 - alpha / (2.0 * n), df=n - 2)
    return float((n - 1) / np.sqrt(n) * np.sqrt(t**2 / (n - 2 + t**2)))


def grubbs_test(data, alpha: float = 0.05) -> GrubbsTestResult:
    r"""Grubbs' test for a single outlier: the most extreme value, in units of the sample standard deviation.

    .. math::

        G = \frac{\max_i|x_i-\bar x|}{s}

    compared against :func:`grubbs_critical_value` (F. E. Grubbs, *Ann.
    Math. Statist.* 21, 27 (1950)). Unlike Dixon's gap-over-range
    :func:`~chemistrykit.analytical.dixon_q_test`, `G` uses every point
    through the mean and standard deviation, and its critical value comes
    from a formula rather than a finite table, so any `n` >= 3 works.

    Parameters
    ----------
    data : array-like of float
        The data set, at least 3 points.
    alpha : float, default 0.05
        Significance level of the two-sided test.

    Returns
    -------
    GrubbsTestResult

    Examples
    --------
    >>> result = grubbs_test([10.0, 10.1, 10.2, 10.3, 15.0])
    >>> result.suspect_value, result.reject
    (15.0, True)
    >>> grubbs_test([10.0, 10.1, 10.2, 10.3, 10.4]).reject
    False
    """
    x = np.asarray(data, dtype=np.float64)
    n = len(x)
    mean = np.mean(x)
    s = np.std(x, ddof=1)
    idx = int(np.argmax(np.abs(x - mean)))
    G = float(abs(x[idx] - mean) / s) if s > 0 else 0.0
    G_crit = grubbs_critical_value(n, alpha)
    return GrubbsTestResult(suspect_value=float(x[idx]), G_statistic=G, G_critical=G_crit, reject=bool(G > G_crit))


def horwitz_rsd(mass_fraction):
    r"""Horwitz's predicted between-laboratory relative standard deviation, in percent.

    .. math::

        \text{RSD}_R(\%) = 2^{\,1-0.5\log_{10}C} = 2\,C^{-0.1505}

    with `C` the analyte concentration as a dimensionless mass fraction
    (1 for a pure substance, :math:`10^{-6}` for 1 ppm) -- the empirical
    "Horwitz trumpet" fitted to thousands of collaborative-study results
    (W. Horwitz, L. R. Kamps, K. W. Boyer, *J. Assoc. Off. Anal. Chem.*
    63, 1344 (1980)): each 100-fold drop in concentration doubles the
    expected RSD.

    Parameters
    ----------
    mass_fraction : float or array-like of float
        Analyte concentration as a mass fraction, in (0, 1].

    Returns
    -------
    float or ndarray

    Examples
    --------
    2% for a pure substance, 16% at 1 ppm:

    >>> round(float(horwitz_rsd(1.0)), 6), round(float(horwitz_rsd(1e-6)), 6)
    (2.0, 16.0)
    """
    C = np.asarray(mass_fraction, dtype=np.float64)
    result = 2.0 ** (1.0 - 0.5 * np.log10(C))
    return float(result) if result.ndim == 0 else result


def horrat(observed_rsd, mass_fraction):
    r"""The HorRat ratio: an observed between-laboratory RSD divided by :func:`horwitz_rsd`.

    Values roughly between 0.5 and 2 are conventionally taken to indicate
    a method with normal, acceptable interlaboratory precision.

    Parameters
    ----------
    observed_rsd : float or array-like of float
        Observed reproducibility RSD, in percent.
    mass_fraction : float or array-like of float
        Analyte concentration as a mass fraction.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> round(float(horrat(8.0, 1e-6)), 6)
    0.5
    """
    result = np.asarray(observed_rsd, dtype=np.float64) / np.asarray(horwitz_rsd(mass_fraction), dtype=np.float64)
    return float(result) if result.ndim == 0 else result
