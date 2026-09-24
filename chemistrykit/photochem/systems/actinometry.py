r"""Chemical actinometry with potassium ferrioxalate (Hatchard & Parker, 1956).

Irradiating acidic :math:`\mathrm{K_3[Fe(C_2O_4)_3]}` photoreduces
Fe(III) to Fe(II) with an accurately known quantum yield (about 1.2 in
the near UV; C. G. Hatchard and C. A. Parker, *Proc. R. Soc. Lond. A*
235, 518 (1956)). The Fe(II) formed is determined spectrophotometrically
as its red 1,10-phenanthroline complex (:math:`\varepsilon_{510}\approx
1.11\times10^4` M\ :sup:`-1` cm\ :sup:`-1`), and the photon flux follows
from

.. math::

    q_p = \frac{n_{\mathrm{Fe^{2+}}}}{\Phi_{\mathrm{Fe^{2+}}}\,t\,(1-10^{-A})},

where :math:`A` is the actinometer's absorbance at the irradiation
wavelength (Kuhn, Braslavsky & Schmidt, "Chemical Actinometry," *Pure
Appl. Chem.* 76, 2105 (2004)).
"""

from __future__ import annotations

import numpy as np

from chemistrykit.spectro.systems.beer_lambert import transmittance

__all__ = ["ferrioxalate_fe2_moles", "ferrioxalate_photon_flux"]


def ferrioxalate_fe2_moles(absorbance_510, volume_L: float, path_length_cm: float = 1.0, epsilon_510: float = 1.11e4):
    r"""Moles of Fe(II) from the absorbance of its phenanthroline complex at 510 nm.

    Beer-Lambert: :math:`n = A_{510} V/(\varepsilon_{510}\ell)`.

    Parameters
    ----------
    absorbance_510 : float or array-like of float
        Absorbance change of the developed solution at 510 nm.
    volume_L : float
        Volume of the developed solution, in L.
    path_length_cm : float, default 1.0
        Cuvette path length, in cm.
    epsilon_510 : float, default 1.11e4
        Molar absorptivity of Fe(phen)3(2+), in M^-1 cm^-1.

    Returns
    -------
    float or ndarray
        Moles of Fe(II).

    Examples
    --------
    >>> round(ferrioxalate_fe2_moles(1.11, volume_L=0.01) * 1e6, 6)
    1.0
    """
    result = np.asarray(absorbance_510, dtype=np.float64) * volume_L / (epsilon_510 * path_length_cm)
    return float(result) if result.ndim == 0 else result


def ferrioxalate_photon_flux(moles_fe2, time_s, absorbance: float = np.inf, quantum_yield: float = 1.21):
    r"""Incident photon flux :math:`q_p = n_{Fe^{2+}}/(\Phi\,t\,(1-10^{-A}))` from ferrioxalate actinometry.

    Parameters
    ----------
    moles_fe2 : float or array-like of float
        Moles of Fe(II) formed (e.g. from :func:`ferrioxalate_fe2_moles`).
    time_s : float or array-like of float
        Irradiation time, in s.
    absorbance : float, default inf
        Actinometer absorbance at the irradiation wavelength (``inf``
        means total absorption, the usual operating regime).
    quantum_yield : float, default 1.21
        Fe(II) quantum yield at the irradiation wavelength (1.21 near
        365 nm per Hatchard & Parker).

    Returns
    -------
    float or ndarray
        Incident photon flux, in einstein/s (mol photons/s).

    Examples
    --------
    1.21 micromol Fe(II) in 100 s under total absorption means
    :math:`10^{-8}` einstein/s:

    >>> round(ferrioxalate_photon_flux(1.21e-6, 100.0) * 1e8, 6)
    1.0
    """
    frac_abs = 1.0 if np.isinf(absorbance) else 1.0 - transmittance(absorbance)
    result = np.asarray(moles_fe2, dtype=np.float64) / (quantum_yield * np.asarray(time_s, dtype=np.float64) * frac_abs)
    return float(result) if result.ndim == 0 else result
