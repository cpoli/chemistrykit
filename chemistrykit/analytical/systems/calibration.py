r"""Linear-regression calibration curves, and IUPAC-convention limits of detection/quantitation.

See Harris, *Quantitative Chemical Analysis*, 9th ed., Ch. 4.5 ("Method
Validation") and Ch. 5, or the IUPAC recommendation (G. L. Long & J. D.
Winefordner, *Anal. Chem.* 55, 712A (1983)), for the standard :math:`3.3
\sigma/m` (LOD) and :math:`10\sigma/m` (LOQ) convention, with `sigma` the
calibration curve's residual standard error and `m` its slope.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.analytical.utils.regression import linear_fit

__all__ = ["LinearCalibration", "fit_calibration"]

#: float: IUPAC LOD multiplier (Long & Winefordner, 1983).
_LOD_MULTIPLIER = 3.3

#: float: IUPAC LOQ multiplier (Long & Winefordner, 1983).
_LOQ_MULTIPLIER = 10.0


@dataclass
class LinearCalibration:
    """A fitted instrument-response-vs-concentration calibration curve, ``signal = slope*conc + intercept``."""

    slope: float
    """float: Sensitivity, signal units per concentration unit."""

    intercept: float
    """float: Signal at zero concentration (ideally the blank signal)."""

    r_squared: float
    """float: Coefficient of determination of the fit."""

    residual_std_error: float
    """float: :math:`s_{y/x}`, the residual standard error about the fit
    (Harris, *Quantitative Chemical Analysis*, 9th ed., Ch. 4.5)."""

    def predict_signal(self, concentration):
        """Predict the instrument signal at given concentration(s).

        Parameters
        ----------
        concentration : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        c = np.asarray(concentration, dtype=np.float64)
        result = self.slope * c + self.intercept
        return float(result) if result.ndim == 0 else result

    def predict_concentration(self, signal):
        r"""Invert the calibration to estimate concentration from a measured signal.

        .. math::

            \hat c = (y_{measured} - b)/m

        Parameters
        ----------
        signal : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        y = np.asarray(signal, dtype=np.float64)
        result = (y - self.intercept) / self.slope
        return float(result) if result.ndim == 0 else result

    def lod(self) -> float:
        r"""Limit of detection, :math:`\text{LOD}=3.3\,s_{y/x}/|m|` (IUPAC convention).

        Returns
        -------
        float
        """
        return _LOD_MULTIPLIER * self.residual_std_error / abs(self.slope)

    def loq(self) -> float:
        r"""Limit of quantitation, :math:`\text{LOQ}=10\,s_{y/x}/|m|` (IUPAC convention).

        Returns
        -------
        float
        """
        return _LOQ_MULTIPLIER * self.residual_std_error / abs(self.slope)


def fit_calibration(concentration, signal) -> LinearCalibration:
    r"""Fit a linear calibration curve (signal vs. concentration) by ordinary least squares.

    Parameters
    ----------
    concentration : array-like of float
        Known standard concentrations (at least 3 distinct values, so the
        residual standard error is defined).
    signal : array-like of float
        Corresponding measured instrument signals.

    Returns
    -------
    LinearCalibration

    Examples
    --------
    LOQ is always exactly :math:`10/3.3` times LOD, by definition,
    regardless of the data:

    >>> conc = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    >>> signal = [0.02, 1.05, 1.98, 3.10, 3.95, 5.08]
    >>> cal = fit_calibration(conc, signal)
    >>> round(cal.loq() / cal.lod(), 6) == round(10.0 / 3.3, 6)
    True

    A perfect (noiseless) calibration has zero residual error and hence
    zero LOD/LOQ:

    >>> cal_perfect = fit_calibration([0.0, 1.0, 2.0, 3.0], [1.0, 3.0, 5.0, 7.0])
    >>> cal_perfect.lod() < 1e-9
    True
    """
    fit = linear_fit(concentration, signal)
    return LinearCalibration(slope=fit.slope, intercept=fit.intercept, r_squared=fit.r_squared, residual_std_error=fit.residual_std_error)
