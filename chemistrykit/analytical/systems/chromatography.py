r"""Chromatographic plate theory: theoretical plates, the van Deemter equation, resolution, and selectivity.

See Harris, *Quantitative Chemical Analysis*, 9th ed., Ch. 23 ("An
Introduction to Chromatographic Separations"), or Skoog, West, Holler &
Crouch, *Fundamentals of Analytical Chemistry*, 9th ed., Ch. 26,
throughout.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.spectro.utils.lineshapes import gaussian

__all__ = [
    "theoretical_plates",
    "plate_height",
    "van_deemter_H",
    "optimum_flow_velocity",
    "minimum_plate_height",
    "retention_factor",
    "selectivity_factor",
    "resolution",
    "simulate_chromatogram",
    "kovats_retention_index",
    "purnell_resolution",
]


def theoretical_plates(retention_time: float, peak_width: float, width_type: str = "base"):
    r"""Number of theoretical plates `N` from a peak's retention time and width.

    For a Gaussian elution peak (Harris, *Quantitative Chemical
    Analysis*, 9th ed., eq. 23.19-23.20):

    .. math::

        N = 16\left(\frac{t_R}{w_{base}}\right)^2 = 5.545\left(\frac{t_R}{w_{1/2}}\right)^2

    where :math:`w_{base}` is the width at the peak base (tangents to the
    inflection points) and :math:`w_{1/2}` is the full width at half
    maximum -- the numeric prefactor differs only because a Gaussian's
    base width (4 standard deviations) and FWHM
    (:math:`2\sqrt{2\ln2}\,\sigma`) are different multiples of `sigma`.

    Parameters
    ----------
    retention_time : float or array-like of float
        Peak retention time :math:`t_R`.
    peak_width : float or array-like of float
        Peak width, in the same units as `retention_time`, of the type
        selected by `width_type`.
    width_type : {"base", "half_height"}, default "base"

    Returns
    -------
    float or ndarray

    Examples
    --------
    A peak eluting at 10.0 min with a 0.5 min base width:

    >>> round(float(theoretical_plates(10.0, 0.5, width_type="base")), 1)
    6400.0

    The two width conventions must agree for a Gaussian peak of
    consistent shape (base width = 4 sigma, FWHM =
    :math:`2\sqrt{2\ln2}\,\sigma\approx2.3548\sigma`):

    >>> sigma = 0.125
    >>> N_base = theoretical_plates(10.0, 4.0 * sigma, width_type="base")
    >>> N_half = theoretical_plates(10.0, 2.3548 * sigma, width_type="half_height")
    >>> bool(abs(N_base - N_half) / N_base < 1e-3)
    True
    """
    tR = np.asarray(retention_time, dtype=np.float64)
    w = np.asarray(peak_width, dtype=np.float64)
    if width_type == "base":
        result = 16.0 * (tR / w) ** 2
    elif width_type == "half_height":
        result = 5.545 * (tR / w) ** 2
    else:
        raise ValueError('width_type must be "base" or "half_height"')
    return float(result) if result.ndim == 0 else result


def plate_height(column_length: float, N):
    r"""Plate height (HETP) :math:`H=L/N`, the column length per theoretical plate.

    Parameters
    ----------
    column_length : float
        Column length `L`, e.g. in cm.
    N : float or array-like of float
        Number of theoretical plates.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> round(float(plate_height(column_length=25.0, N=6400.0)), 6)
    0.003906
    """
    N = np.asarray(N, dtype=np.float64)
    result = column_length / N
    return float(result) if result.ndim == 0 else result


def van_deemter_H(u, A: float, B: float, C: float):
    r"""The van Deemter equation: plate height `H` as a function of mobile-phase linear velocity `u`.

    .. math::

        H = A + \frac{B}{u} + Cu

    with `A` the eddy-diffusion term (velocity-independent band
    broadening from unequal flow paths), `B` the longitudinal-molecular-
    diffusion term (dominant at low `u`), and `C` the mass-transfer-
    resistance term (dominant at high `u`) (J. J. van Deemter, F. J.
    Zuiderweg, A. Klinkenberg, *Chem. Eng. Sci.* 5, 271 (1956); Harris,
    *Quantitative Chemical Analysis*, 9th ed., Ch. 23.4).

    Parameters
    ----------
    u : float or array-like of float
        Mobile-phase linear velocity.
    A, B, C : float
        Van Deemter coefficients (all non-negative for a physical column).

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> round(float(van_deemter_H(u=2.0, A=1.0, B=2.0, C=0.05)), 4)
    2.1
    """
    u = np.asarray(u, dtype=np.float64)
    result = A + B / u + C * u
    return float(result) if result.ndim == 0 else result


def optimum_flow_velocity(B: float, C: float) -> float:
    r"""The flow velocity :math:`u_{opt}=\sqrt{B/C}` minimizing the van Deemter equation.

    Found by :math:`dH/du=-B/u^2+C=0` (Harris, *Quantitative Chemical
    Analysis*, 9th ed., Ch. 23.4).

    Parameters
    ----------
    B, C : float
        Van Deemter longitudinal-diffusion and mass-transfer coefficients.

    Returns
    -------
    float

    Examples
    --------
    >>> round(optimum_flow_velocity(B=2.0, C=0.05), 6)
    6.324555
    """
    return float(np.sqrt(B / C))


def minimum_plate_height(A: float, B: float, C: float) -> float:
    r"""The minimum plate height :math:`H_{min}=A+2\sqrt{BC}`, at the optimum flow velocity.

    Substituting :func:`optimum_flow_velocity` into
    :func:`van_deemter_H` (Harris, *Quantitative Chemical Analysis*, 9th
    ed., Ch. 23.4).

    Parameters
    ----------
    A, B, C : float

    Returns
    -------
    float

    Examples
    --------
    Matches direct evaluation of :func:`van_deemter_H` at the optimum
    velocity from :func:`optimum_flow_velocity`:

    >>> A, B, C = 1.0, 2.0, 0.05
    >>> u_opt = optimum_flow_velocity(B, C)
    >>> bool(round(minimum_plate_height(A, B, C), 6) == round(float(van_deemter_H(u_opt, A, B, C)), 6))
    True
    """
    return A + 2.0 * np.sqrt(B * C)


def retention_factor(retention_time: float, dead_time: float):
    r"""The retention (capacity) factor :math:`k=(t_R-t_0)/t_0`.

    Parameters
    ----------
    retention_time : float or array-like of float
        Analyte retention time :math:`t_R`.
    dead_time : float
        Column dead time (retention time of an unretained species) :math:`t_0`.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> round(float(retention_factor(retention_time=12.0, dead_time=2.0)), 6)
    5.0
    """
    tR = np.asarray(retention_time, dtype=np.float64)
    result = (tR - dead_time) / dead_time
    return float(result) if result.ndim == 0 else result


def selectivity_factor(k1: float, k2: float) -> float:
    r"""The selectivity (relative retention) factor :math:`\alpha=k_2/k_1` (:math:`k_2\ge k_1` by convention).

    Parameters
    ----------
    k1, k2 : float
        Retention factors of the earlier- and later-eluting peaks.

    Returns
    -------
    float

    Examples
    --------
    >>> round(selectivity_factor(k1=2.0, k2=5.0), 6)
    2.5
    """
    return k2 / k1


def resolution(tR1: float, tR2: float, w1: float, w2: float) -> float:
    r"""Chromatographic resolution :math:`R_s=2(t_{R,2}-t_{R,1})/(w_1+w_2)` between two adjacent peaks.

    :math:`R_s\ge1.5` is the conventional criterion for baseline
    separation (Harris, *Quantitative Chemical Analysis*, 9th ed., eq.
    23.21).

    Parameters
    ----------
    tR1, tR2 : float
        Retention times of the earlier- and later-eluting peaks.
    w1, w2 : float
        Base widths of the two peaks, same time units.

    Returns
    -------
    float

    Examples
    --------
    >>> round(resolution(tR1=9.0, tR2=10.0, w1=0.5, w2=0.5), 6)
    2.0
    """
    return 2.0 * (tR2 - tR1) / (w1 + w2)


def simulate_chromatogram(t, centers, retention_times=None, N: float = 10000.0, amplitudes=None):
    r"""Simulate a chromatogram as a sum of Gaussian elution peaks of plate-count-consistent width.

    Reuses :func:`chemistrykit.spectro.utils.lineshapes.gaussian` (each
    elution peak is, to a good approximation, Gaussian -- the same
    plate-theory result that underlies :func:`theoretical_plates`), with
    each peak's FWHM set by :math:`w_{1/2}=t_R\sqrt{8\ln2/N}` (inverting
    :func:`theoretical_plates`'s half-height form,
    :math:`N=5.545(t_R/w_{1/2})^2` with :math:`5.545\approx8\ln2`).

    Parameters
    ----------
    t : array-like of float
        Time grid to evaluate the chromatogram on.
    centers : array-like of float
        Retention time of each peak (alias for `retention_times`, kept
        for a natural call signature with :func:`chemistrykit.spectro.utils.lineshapes.broaden_stick_spectrum`).
    retention_times : array-like of float, optional
        If given, overrides `centers`.
    N : float, default 10000.0
        Number of theoretical plates (assumed equal for every peak, the
        common simplifying assumption for a single column/method).
    amplitudes : array-like of float, optional
        Relative peak heights; defaults to 1.0 for every peak.

    Returns
    -------
    ndarray, shape matching `t`

    Examples
    --------
    A single simulated peak's apparent plate count (recovered from its
    numerically-measured FWHM) matches the `N` it was built from:

    >>> import numpy as np
    >>> t = np.linspace(8.0, 12.0, 200001)
    >>> N_true = 10000.0
    >>> chrom = simulate_chromatogram(t, centers=[10.0], N=N_true)
    >>> half_max = chrom.max() / 2.0
    >>> above = t[chrom >= half_max]
    >>> fwhm = above.max() - above.min()
    >>> N_recovered = theoretical_plates(10.0, fwhm, width_type="half_height")
    >>> bool(abs(N_recovered - N_true) / N_true < 0.01)
    True
    """
    t = np.asarray(t, dtype=np.float64)
    centers = np.asarray(retention_times if retention_times is not None else centers, dtype=np.float64)
    amps = np.ones_like(centers) if amplitudes is None else np.asarray(amplitudes, dtype=np.float64)
    chrom = np.zeros_like(t)
    for tR, amp in zip(centers, amps, strict=True):
        fwhm = tR * np.sqrt(8.0 * np.log(2.0) / N)
        chrom += amp * gaussian(t, tR, fwhm)
    return chrom


def kovats_retention_index(t_x, t_n, t_N, n: int, N: int | None = None, dead_time: float = 0.0):
    r"""Kováts retention index of an analyte, bracketed by two n-alkane standards.

    Under isothermal conditions the logarithm of an n-alkane's adjusted
    retention time :math:`t'=t_R-t_0` grows linearly with its carbon
    number, so an analyte eluting between the n-alkanes with `n` and `N`
    carbons is assigned the logarithmically interpolated index (E. Kováts,
    *Helv. Chim. Acta* 41, 1915 (1958)):

    .. math::

        I = 100\left[n + (N-n)\,\frac{\log t'_x-\log t'_n}{\log t'_N-\log t'_n}\right]

    With the usual choice :math:`N=n+1` this is Kováts' original
    definition; every n-alkane itself has, by construction, :math:`I=100n`.

    Parameters
    ----------
    t_x : float or array-like of float
        Retention time(s) of the analyte.
    t_n, t_N : float
        Retention times of the smaller (`n` carbons) and larger (`N`
        carbons) bracketing n-alkanes.
    n : int
        Carbon number of the earlier-eluting n-alkane.
    N : int, optional
        Carbon number of the later-eluting n-alkane; defaults to ``n + 1``.
    dead_time : float, default 0.0
        Column dead time :math:`t_0`, subtracted from every retention time
        to form the adjusted retention times :math:`t'`.

    Returns
    -------
    float or ndarray

    Examples
    --------
    An analyte eluting at the geometric mean of the adjusted retention
    times of n-octane and n-nonane sits exactly halfway, at `I` = 850:

    >>> import numpy as np
    >>> t0, t8, t9 = 1.0, 5.0, 9.0
    >>> t_x = t0 + np.sqrt((t8 - t0) * (t9 - t0))
    >>> round(float(kovats_retention_index(t_x, t8, t9, n=8, dead_time=t0)), 6)
    850.0
    """
    if N is None:
        N = n + 1
    tx = np.asarray(t_x, dtype=np.float64) - dead_time
    tn = t_n - dead_time
    tN = t_N - dead_time
    result = 100.0 * (n + (N - n) * (np.log10(tx) - np.log10(tn)) / (np.log10(tN) - np.log10(tn)))
    return float(result) if result.ndim == 0 else result


def purnell_resolution(N, alpha: float, k2: float):
    r"""Purnell's resolution equation: resolution in terms of efficiency, selectivity, and retention.

    .. math::

        R_s = \frac{\sqrt N}{4}\,\frac{\alpha-1}{\alpha}\,\frac{k_2}{1+k_2}

    (J. H. Purnell, *J. Chem. Soc.* 1960, 1268), with `N` the plate count,
    :math:`\alpha=k_2/k_1` the selectivity factor, and :math:`k_2` the
    retention factor of the later-eluting peak. It is exact when both
    peaks have the base width :math:`4t_{R,2}/\sqrt N` of the second peak,
    and separates the three independent levers a chromatographer can pull:
    efficiency (:math:`\sqrt N`), selectivity, and retention.

    Parameters
    ----------
    N : float or array-like of float
        Number of theoretical plates.
    alpha : float
        Selectivity factor :math:`\alpha\ge1`.
    k2 : float
        Retention factor of the later-eluting peak.

    Returns
    -------
    float or ndarray

    Examples
    --------
    Agrees with the direct :func:`resolution` formula for two peaks of
    equal base width :math:`4t_{R,2}/\sqrt N`:

    >>> import numpy as np
    >>> N, t0, k1, k2 = 10000.0, 1.0, 4.0, 4.4
    >>> tR1, tR2 = t0 * (1 + k1), t0 * (1 + k2)
    >>> w = 4.0 * tR2 / np.sqrt(N)
    >>> bool(abs(purnell_resolution(N, alpha=k2 / k1, k2=k2) - resolution(tR1, tR2, w, w)) < 1e-12)
    True

    Resolution grows only as :math:`\sqrt N`: quadrupling the plate count
    doubles it.

    >>> round(float(purnell_resolution(40000.0, 1.1, 4.4) / purnell_resolution(10000.0, 1.1, 4.4)), 6)
    2.0
    """
    N = np.asarray(N, dtype=np.float64)
    result = np.sqrt(N) / 4.0 * (alpha - 1.0) / alpha * k2 / (1.0 + k2)
    return float(result) if result.ndim == 0 else result
