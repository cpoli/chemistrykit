r"""Butler-Volmer electrode kinetics, exchange current density, and Tafel-plot linearization.

See Bard & Faulkner, *Electrochemical Methods: Fundamentals and
Applications*, 2nd ed., Ch. 3.3-3.4, throughout.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.constants import FARADAY, STANDARD_TEMPERATURE, R
from chemistrykit.electrochem.utils.regression import LinearFit, linear_fit

__all__ = [
    "butler_volmer_current_density",
    "exchange_current_density",
    "tafel_slope",
    "tafel_overpotential",
    "TafelFit",
    "fit_tafel_plot",
]


def butler_volmer_current_density(i0: float, eta, alpha: float = 0.5, n: int = 1, T: float = STANDARD_TEMPERATURE, R_gas: float = R, F: float = FARADAY):
    r"""The Butler-Volmer equation for net electrode current density.

    .. math::

        i = i_0\left[\exp\!\left(\frac{\alpha n F\eta}{RT}\right)
                     -\exp\!\left(-\frac{(1-\alpha)nF\eta}{RT}\right)\right]

    where :math:`\eta = E - E_{eq}` is the overpotential, :math:`i_0` the
    exchange current density, and :math:`\alpha` the *anodic* charge
    transfer coefficient, so :math:`1-\alpha` is the cathodic one (Bard &
    Faulkner, *Electrochemical Methods*, 2nd ed., Ch. 3.4, eq. 3.4.11,
    written there with anodic current positive and cathodic coefficient
    :math:`\alpha_c=1-\alpha`). The same anodic convention is used by
    every function in this module. The first term is the anodic
    (oxidation) partial current, the second the cathodic (reduction)
    partial current; at :math:`\eta=0` they exactly cancel, so no net
    current flows at equilibrium even though both partial reactions are
    still occurring (a genuinely dynamic equilibrium, unlike a
    thermodynamic "nothing is happening").

    Parameters
    ----------
    i0 : float
        Exchange current density, in A/m^2 (or any consistent current-
        density unit; the same unit is returned).
    eta : float or array-like of float
        Overpotential, in V.
    alpha : float, default 0.5
        Anodic charge transfer (symmetry) coefficient, in :math:`(0, 1)`;
        0.5 is the common approximation for a simple, symmetric
        one-electron step.
    n : int, default 1
        Number of electrons transferred in the rate-determining step.
    T : float, default :data:`chemistrykit.constants.STANDARD_TEMPERATURE`
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.
    F : float, default :data:`chemistrykit.constants.FARADAY`
        Faraday constant, in C/mol.

    Returns
    -------
    float or ndarray

    Examples
    --------
    At zero overpotential the net current is exactly zero:

    >>> round(float(butler_volmer_current_density(i0=1e-6, eta=0.0)), 12)
    0.0

    A positive overpotential drives net anodic (positive) current:

    >>> bool(butler_volmer_current_density(1e-6, eta=0.1) > 0)
    True
    """
    eta = np.asarray(eta, dtype=np.float64)
    f = n * F / (R_gas * T)
    result = i0 * (np.exp(alpha * f * eta) - np.exp(-(1.0 - alpha) * f * eta))
    return float(result) if result.ndim == 0 else result


def exchange_current_density(k0: float, C_ox: float, C_red: float, n: int = 1, alpha: float = 0.5, F: float = FARADAY) -> float:
    r"""Exchange current density from a standard heterogeneous rate constant and bulk concentrations.

    .. math::

        i_0 = nFk^0 C_{ox}^{\alpha}C_{red}^{1-\alpha}

    (Bard & Faulkner, *Electrochemical Methods*, 2nd ed., Ch. 3.4, eq.
    3.4.6, :math:`i_0=nFk^0C_{ox}^{1-\alpha_c}C_{red}^{\alpha_c}`,
    rewritten with this module's anodic :math:`\alpha=1-\alpha_c`) -- the current density that flows equally in both directions
    at equilibrium, before any net overpotential is applied.

    Parameters
    ----------
    k0 : float
        Standard heterogeneous electron-transfer rate constant, in m/s
        (or a consistent length/time unit).
    C_ox, C_red : float
        Bulk concentrations of the oxidized and reduced forms, in
        mol/m^3 (consistent with `k0`'s length unit).
    n : int, default 1
        Electrons transferred.
    alpha : float, default 0.5
        Anodic charge transfer coefficient, as in
        :func:`butler_volmer_current_density`.
    F : float, default :data:`chemistrykit.constants.FARADAY`
        Faraday constant, in C/mol.

    Returns
    -------
    float
        Exchange current density, in A/m^2.

    Examples
    --------
    Symmetric concentrations and alpha=0.5 give a simple, order-1-free-of-composition result:

    >>> i0 = exchange_current_density(k0=1e-5, C_ox=1.0, C_red=1.0, n=1)
    >>> round(i0, 6)
    0.964853
    """
    return n * F * k0 * C_ox**alpha * C_red ** (1.0 - alpha)


def tafel_slope(alpha: float, n: int = 1, T: float = STANDARD_TEMPERATURE, R_gas: float = R, F: float = FARADAY, branch: str = "anodic") -> float:
    r"""The Tafel slope :math:`b`, in V per decade of current.

    From linearizing the dominant (far-from-equilibrium) exponential
    term of the Butler-Volmer equation and rewriting in terms of
    :math:`\log_{10}i` (Bard & Faulkner, *Electrochemical Methods*, 2nd
    ed., Ch. 3.4, eq. 3.4.11-3.4.13):

    .. math::

        b_{anodic} = \frac{2.303RT}{\alpha n F}, \qquad
        b_{cathodic} = \frac{2.303RT}{(1-\alpha) n F}

    Parameters
    ----------
    alpha : float
        Anodic charge transfer coefficient, as in
        :func:`butler_volmer_current_density`.
    n : int, default 1
        Electrons transferred.
    T : float, default :data:`chemistrykit.constants.STANDARD_TEMPERATURE`
    R_gas : float, default :data:`chemistrykit.constants.R`
    F : float, default :data:`chemistrykit.constants.FARADAY`
    branch : {"anodic", "cathodic"}
        Which high-overpotential branch's slope to return.

    Returns
    -------
    float
        Tafel slope, in V/decade (positive for the anodic branch,
        negative for the cathodic one, matching the sign convention
        that :math:`\eta` increases with :math:`\log_{10}i` on the
        anodic branch and decreases with it on the cathodic one).

    Examples
    --------
    With `alpha=0.5` and a one-electron step at 25 degC, the anodic and
    cathodic slopes have equal magnitude (the symmetric-barrier case):

    >>> b_a = tafel_slope(alpha=0.5, n=1, branch="anodic")
    >>> b_c = tafel_slope(alpha=0.5, n=1, branch="cathodic")
    >>> round(b_a, 4), round(b_c, 4)
    (0.1183, -0.1183)
    """
    if branch == "anodic":
        return 2.303 * R_gas * T / (alpha * n * F)
    elif branch == "cathodic":
        return -2.303 * R_gas * T / ((1.0 - alpha) * n * F)
    raise ValueError("branch must be 'anodic' or 'cathodic'")


def tafel_overpotential(
    i, i0: float, alpha: float = 0.5, n: int = 1, T: float = STANDARD_TEMPERATURE, R_gas: float = R, F: float = FARADAY, branch: str = "anodic"
):
    r"""The linearized Tafel-equation overpotential at high :math:`|\eta|`.

    :math:`\eta = b\log_{10}(i/i_0)`, dropping the Butler-Volmer
    equation's back-reaction exponential term (valid once
    :math:`|\eta|\gtrsim100/n` mV, where that term is negligible
    compared to the dominant one -- Bard & Faulkner, *Electrochemical
    Methods*, 2nd ed., Ch. 3.4). **Approximation flagged**: this is the
    high-overpotential limit of :func:`butler_volmer_current_density`,
    not the exact equation -- it is inaccurate near :math:`\eta=0` where
    both exponential terms matter (see the module examples gallery for a
    numerical comparison across overpotential).

    Parameters
    ----------
    i : float or array-like of float
        Current density (same sign/branch as `branch`), in the same
        units as `i0`.
    i0 : float
        Exchange current density.
    alpha : float, default 0.5
    n : int, default 1
    T : float, default :data:`chemistrykit.constants.STANDARD_TEMPERATURE`
    R_gas : float, default :data:`chemistrykit.constants.R`
    F : float, default :data:`chemistrykit.constants.FARADAY`
    branch : {"anodic", "cathodic"}

    Returns
    -------
    float or ndarray
        Overpotential, in V.

    Examples
    --------
    At high overpotential the Tafel approximation agrees closely with
    the full Butler-Volmer equation (checked here to better than 0.5%
    at 300 mV):

    >>> i0, alpha, n = 1e-6, 0.5, 1
    >>> eta_true = 0.30
    >>> i_full = butler_volmer_current_density(i0, eta_true, alpha=alpha, n=n)
    >>> eta_tafel = tafel_overpotential(i_full, i0, alpha=alpha, n=n, branch="anodic")
    >>> bool(abs(eta_tafel - eta_true) / eta_true < 0.005)
    True
    """
    b = tafel_slope(alpha=alpha, n=n, T=T, R_gas=R_gas, F=F, branch=branch)
    i = np.asarray(i, dtype=np.float64)
    result = b * np.log10(np.abs(i) / i0)
    return float(result) if result.ndim == 0 else result


@dataclass
class TafelFit:
    r"""Result of fitting :math:`\log_{10}|i|` vs. :math:`\eta` data to a Tafel line.

    Since :math:`\eta = b\log_{10}i - b\log_{10}i_0`, a plot of `eta`
    against :math:`\log_{10}|i|` is a straight line of slope `b` (the
    Tafel slope) and intercept :math:`-b\log_{10}i_0`, from which the
    exchange current density :math:`i_0=10^{-\text{intercept}/b}`
    follows.
    """

    tafel_slope: float
    """float: Fitted Tafel slope `b`, in V/decade."""

    exchange_current_density: float
    """float: Fitted exchange current density :math:`i_0`."""

    r_squared: float
    """float: Coefficient of determination of the linear fit."""

    def predict(self, i):
        """Evaluate the fitted Tafel line's overpotential at current density(-ies) `i`.

        Parameters
        ----------
        i : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        i = np.asarray(i, dtype=np.float64)
        return self.tafel_slope * np.log10(np.abs(i) / self.exchange_current_density)


def fit_tafel_plot(eta, i) -> TafelFit:
    r"""Fit high-overpotential (`eta`, `i`) data to a Tafel line and recover `b` and :math:`i_0`.

    Linearizes :math:`\eta = b\log_{10}i - b\log_{10}i_0` and fits by
    ordinary least squares against :math:`x=\log_{10}|i|` -- the
    experimental Tafel-plot analysis (Bard & Faulkner, *Electrochemical
    Methods*, 2nd ed., Ch. 3.4), structurally identical to
    :func:`chemistrykit.kinetics.systems.arrhenius.fit_arrhenius`'s
    Arrhenius plot.

    Parameters
    ----------
    eta : array-like of float
        Overpotentials, in V (should be restricted to the high-`|eta|`
        Tafel regime for the linear approximation to hold).
    i : array-like of float
        Corresponding current densities.

    Returns
    -------
    TafelFit

    Examples
    --------
    Generate exact high-overpotential data from a known :math:`(b, i_0)`
    and recover both:

    >>> import numpy as np
    >>> i0_true, alpha, n = 2e-6, 0.5, 1
    >>> eta = np.linspace(0.2, 0.4, 10)
    >>> i = butler_volmer_current_density(i0_true, eta, alpha=alpha, n=n)
    >>> fit = fit_tafel_plot(eta, i)
    >>> round(fit.exchange_current_density / i0_true, 2)
    1.0
    """
    eta = np.asarray(eta, dtype=np.float64)
    i = np.asarray(i, dtype=np.float64)
    fit: LinearFit = linear_fit(np.log10(np.abs(i)), eta)
    b = fit.slope
    i0 = 10.0 ** (-fit.intercept / b)
    return TafelFit(tafel_slope=float(b), exchange_current_density=float(i0), r_squared=fit.r_squared)
