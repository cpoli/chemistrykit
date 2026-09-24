r"""Excitation energy transfer: Förster (dipole-dipole) and Dexter (exchange) mechanisms.

* Förster resonance energy transfer (FRET; Th. Förster, *Ann. Phys.*
  437, 55 (1948)) -- long-range dipole-dipole coupling with rate
  :math:`k_T=(1/\tau_D)(R_0/r)^6`. See Lakowicz, *Principles of
  Fluorescence Spectroscopy*, 3rd ed., Ch. 13.
* Dexter exchange transfer (D. L. Dexter, *J. Chem. Phys.* 21, 836
  (1953)) -- short-range, orbital-overlap-mediated transfer with rate
  :math:`k_{ET}=KJ\exp(-2r/L)`. See Turro, Ramamurthy & Scaiano, *Modern
  Molecular Photochemistry of Organic Molecules*, Ch. 7.
"""

from __future__ import annotations

import numpy as np

__all__ = ["forster_radius", "forster_rate", "forster_efficiency", "dexter_rate"]


def forster_radius(kappa2: float, n: float, quantum_yield_donor: float, overlap_J: float) -> float:
    r"""Förster distance :math:`R_0 = 0.211\,(\kappa^2 n^{-4} Q_D J)^{1/6}` Å.

    The donor-acceptor distance at which transfer is 50% efficient
    (Lakowicz, 3rd ed., eq. 13.12), with the spectral overlap integral
    :math:`J` in M\ :sup:`-1` cm\ :sup:`-1` nm\ :sup:`4`.

    Parameters
    ----------
    kappa2 : float
        Orientation factor :math:`\kappa^2` (2/3 for isotropic averaging).
    n : float
        Refractive index of the medium.
    quantum_yield_donor : float
        Donor fluorescence quantum yield without acceptor.
    overlap_J : float
        Spectral overlap integral, in M^-1 cm^-1 nm^4.

    Returns
    -------
    float
        :math:`R_0`, in Å.

    Examples
    --------
    A typical dye pair (:math:`J=10^{15}`, :math:`Q_D=0.5`, :math:`n=1.4`)
    gives :math:`R_0` of about 4.4 nm:

    >>> round(forster_radius(2 / 3, 1.4, 0.5, 1.0e15), 1)
    44.4
    """
    return 0.211 * (kappa2 * n**-4 * quantum_yield_donor * overlap_J) ** (1.0 / 6.0)


def forster_rate(r, R0: float, tau_D: float):
    r"""Förster transfer rate :math:`k_T = \tau_D^{-1}(R_0/r)^6`.

    Parameters
    ----------
    r : float or array-like of float
        Donor-acceptor distance (same unit as `R0`).
    R0 : float
        Förster distance.
    tau_D : float
        Donor lifetime without acceptor.

    Returns
    -------
    float or ndarray
        Transfer rate constant, in 1/(unit of `tau_D`).

    Examples
    --------
    At :math:`r=R_0` transfer competes equally with donor decay:

    >>> forster_rate(5.0, R0=5.0, tau_D=2.0)
    0.5
    """
    result = (1.0 / tau_D) * (R0 / np.asarray(r, dtype=np.float64)) ** 6
    return float(result) if result.ndim == 0 else result


def forster_efficiency(r, R0: float):
    r"""FRET efficiency :math:`E = 1/(1+(r/R_0)^6)`.

    Parameters
    ----------
    r : float or array-like of float
        Donor-acceptor distance.
    R0 : float
        Förster distance (same unit as `r`).

    Returns
    -------
    float or ndarray
        Efficiency in :math:`[0,1]`.

    Examples
    --------
    >>> forster_efficiency(5.0, R0=5.0)
    0.5
    >>> round(forster_efficiency(10.0, R0=5.0), 5)
    0.01538
    """
    result = 1.0 / (1.0 + (np.asarray(r, dtype=np.float64) / R0) ** 6)
    return float(result) if result.ndim == 0 else result


def dexter_rate(r, K: float, J: float, L: float):
    r"""Dexter exchange energy-transfer rate :math:`k_{ET} = KJ\exp(-2r/L)`.

    Exchange transfer requires the donor's and acceptor's electron
    clouds to overlap, so the rate decays exponentially with separation
    rather than as the :math:`r^{-6}` power law of Förster transfer, and
    it allows spin-forbidden triplet-triplet transfer that dipole
    coupling does not (Dexter 1953; Turro et al., Ch. 7).

    Parameters
    ----------
    r : float or array-like of float
        Donor-acceptor edge-to-edge distance.
    K : float
        Pre-exponential factor (orbital interaction strength).
    J : float
        Normalized spectral overlap integral.
    L : float
        Effective average van der Waals radius (same unit as `r`).

    Returns
    -------
    float or ndarray

    Examples
    --------
    Each increase of :math:`r` by :math:`L/2` reduces the rate by a
    factor of :math:`e`:

    >>> k1, k2 = dexter_rate(1.0, K=1.0, J=1.0, L=2.0), dexter_rate(2.0, K=1.0, J=1.0, L=2.0)
    >>> round(k1 / k2, 6) == round(float(np.e), 6)
    True
    """
    result = K * J * np.exp(-2.0 * np.asarray(r, dtype=np.float64) / L)
    return float(result) if result.ndim == 0 else result
