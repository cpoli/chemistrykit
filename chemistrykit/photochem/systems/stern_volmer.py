r"""Stern-Volmer fluorescence quenching, and distinguishing static from dynamic quenching.

See Lakowicz, *Principles of Fluorescence Spectroscopy*, 3rd ed., Ch. 8,
throughout. The original: O. Stern & M. Volmer, *Phys. Z.* 20, 183
(1919).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.photochem.utils.regression import LinearFit, linear_fit

__all__ = [
    "stern_volmer_ratio",
    "dynamic_quenching_constant",
    "SternVolmerFit",
    "fit_stern_volmer",
    "classify_quenching_mechanism",
]


def stern_volmer_ratio(Ksv: float, Q):
    r"""The Stern-Volmer equation: :math:`I_0/I = 1 + K_{sv}[Q]`.

    A plot of the fluorescence-intensity ratio :math:`I_0/I` (unquenched
    over quenched intensity) against quencher concentration :math:`[Q]`
    is linear with slope :math:`K_{sv}`, the Stern-Volmer constant
    (Lakowicz, *Principles of Fluorescence Spectroscopy*, 3rd ed., Ch. 8,
    eq. 8.4). For *purely dynamic* (collisional) quenching,
    :math:`K_{sv}=k_q\tau_0` (see :func:`dynamic_quenching_constant`) and
    this same linear relationship also holds for the fluorescence
    *lifetime* ratio :math:`\tau_0/\tau`; for *purely static* quenching
    (ground-state complex formation) the intensity ratio is still linear
    in :math:`[Q]` (with :math:`K_{sv}` reinterpreted as the complex's
    formation constant) but the lifetime ratio stays exactly 1 -- see
    :func:`classify_quenching_mechanism`.

    Parameters
    ----------
    Ksv : float
        Stern-Volmer constant, in units inverse to `Q`'s (e.g. 1/M for
        `Q` in mol/L).
    Q : float or array-like of float
        Quencher concentration.

    Returns
    -------
    float or ndarray
        :math:`I_0/I` (or, for a dynamic mechanism, equally
        :math:`\tau_0/\tau`).

    Examples
    --------
    At zero quencher, the ratio is exactly 1 (no quenching):

    >>> round(float(stern_volmer_ratio(Ksv=5.0, Q=0.0)), 6)
    1.0

    The ratio grows linearly with quencher concentration:

    >>> Q = np.array([0.0, 0.01, 0.02, 0.03])
    >>> ratio = stern_volmer_ratio(Ksv=20.0, Q=Q)
    >>> bool(np.allclose(np.diff(ratio), 20.0 * np.diff(Q)))
    True
    """
    result = 1.0 + Ksv * np.asarray(Q, dtype=np.float64)
    return float(result) if result.ndim == 0 else result


def dynamic_quenching_constant(kq: float, tau0: float) -> float:
    r"""The Stern-Volmer constant for purely dynamic (collisional) quenching: :math:`K_{sv} = k_q\tau_0`.

    :math:`k_q` is the bimolecular quenching rate constant (an
    encounter/collision frequency, at most diffusion-controlled) and
    :math:`\tau_0` the unquenched excited-state lifetime -- a longer-
    lived excited state has more opportunity to be collisionally
    quenched before it decays on its own, hence the product form
    (Lakowicz, *Principles of Fluorescence Spectroscopy*, 3rd ed., Ch.
    8, eq. 8.5).

    Parameters
    ----------
    kq : float
        Bimolecular quenching rate constant, in 1/(M*s).
    tau0 : float
        Unquenched excited-state lifetime, in s.

    Returns
    -------
    float
        :math:`K_{sv}`, in 1/M.

    Examples
    --------
    >>> round(dynamic_quenching_constant(kq=1.0e10, tau0=5.0e-9), 4)
    50.0
    """
    return kq * tau0


@dataclass
class SternVolmerFit:
    """Result of fitting intensity-ratio-vs-quencher-concentration data to the Stern-Volmer equation."""

    Ksv: float
    """float: Fitted Stern-Volmer constant."""

    r_squared: float
    """float: Coefficient of determination of the linear fit."""

    def predict(self, Q):
        """Evaluate the fitted Stern-Volmer line at quencher concentration(s) `Q`.

        Parameters
        ----------
        Q : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        return stern_volmer_ratio(self.Ksv, Q)


def fit_stern_volmer(Q, intensity_ratio) -> SternVolmerFit:
    r"""Fit :math:`I_0/I` (or :math:`\tau_0/\tau`) vs. :math:`[Q]` data to the Stern-Volmer equation.

    Fits :math:`I_0/I = 1 + K_{sv}[Q]` by ordinary least squares,
    constraining the intercept to 1 exactly (the physically required
    value at :math:`[Q]=0`) by fitting the *shift* :math:`I_0/I - 1`
    through the origin against `Q` -- structurally similar to
    :func:`chemistrykit.kinetics.systems.arrhenius.fit_arrhenius`'s
    linearized fit, but with a fixed rather than free intercept.

    Parameters
    ----------
    Q : array-like of float
        Quencher concentrations (at least 2 distinct values, including a
        nonzero one).
    intensity_ratio : array-like of float
        Corresponding measured :math:`I_0/I` (or :math:`\tau_0/\tau`).

    Returns
    -------
    SternVolmerFit

    Examples
    --------
    Generate exact data from a known :math:`K_{sv}` and recover it:

    >>> Q = np.array([0.0, 0.005, 0.010, 0.020, 0.040])
    >>> ratio = stern_volmer_ratio(Ksv=35.0, Q=Q)
    >>> fit = fit_stern_volmer(Q, ratio)
    >>> round(fit.Ksv, 6)
    35.0
    >>> round(fit.r_squared, 6)
    1.0
    """
    Q = np.asarray(Q, dtype=np.float64)
    intensity_ratio = np.asarray(intensity_ratio, dtype=np.float64)
    fit: LinearFit = linear_fit(Q, intensity_ratio)
    return SternVolmerFit(Ksv=fit.slope, r_squared=fit.r_squared)


def classify_quenching_mechanism(intensity_ratio_slope: float, lifetime_ratio_slope: float, rtol: float = 0.1) -> str:
    r"""Classify a quenching mechanism as dynamic, static, or mixed from the two Stern-Volmer slopes.

    Measuring *both* the steady-state intensity ratio's Stern-Volmer
    slope and the excited-state lifetime ratio's is the standard
    diagnostic (Lakowicz, *Principles of Fluorescence Spectroscopy*, 3rd
    ed., Ch. 8.2): pure dynamic (collisional) quenching reduces the
    lifetime and the intensity by the same factor, so the two slopes
    agree (:math:`(I_0/I)` and :math:`(\tau_0/\tau)` are superimposable);
    pure static quenching (ground-state complexation) removes a fraction
    of fluorophores from the *population* entirely without shortening
    the lifetime of the ones that remain uncomplexed, so
    :math:`\tau_0/\tau=1` (slope 0) while :math:`I_0/I` still rises
    linearly. A slope ratio between these extremes indicates both
    mechanisms are contributing.

    Parameters
    ----------
    intensity_ratio_slope : float
        Fitted Stern-Volmer slope of :math:`I_0/I` vs. :math:`[Q]`.
    lifetime_ratio_slope : float
        Fitted Stern-Volmer slope of :math:`\tau_0/\tau` vs. :math:`[Q]`.
    rtol : float, default 0.1
        Relative tolerance for judging the two slopes "equal" (dynamic)
        vs. the lifetime slope "zero" (static).

    Returns
    -------
    str
        One of ``"dynamic"``, ``"static"``, or ``"mixed"``.

    Examples
    --------
    Equal slopes indicate purely dynamic quenching:

    >>> classify_quenching_mechanism(intensity_ratio_slope=40.0, lifetime_ratio_slope=40.0)
    'dynamic'

    A vanishing lifetime slope with a nonzero intensity slope indicates
    purely static quenching:

    >>> classify_quenching_mechanism(intensity_ratio_slope=40.0, lifetime_ratio_slope=0.0)
    'static'

    Something in between indicates a mixed mechanism:

    >>> classify_quenching_mechanism(intensity_ratio_slope=40.0, lifetime_ratio_slope=20.0)
    'mixed'
    """
    if abs(lifetime_ratio_slope) <= rtol * abs(intensity_ratio_slope):
        return "static"
    if abs(intensity_ratio_slope - lifetime_ratio_slope) <= rtol * abs(intensity_ratio_slope):
        return "dynamic"
    return "mixed"
