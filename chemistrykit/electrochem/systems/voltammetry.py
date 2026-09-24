r"""Diffusion-limited currents in electroanalysis: chronoamperometry, polarography, and voltammetry.

See Bard & Faulkner, *Electrochemical Methods: Fundamentals and
Applications*, 2nd ed., Ch. 5.2 (the Cottrell equation), Ch. 7.1 (the
dropping mercury electrode and the Ilkovič equation), Ch. 5.4 (the
reversible wave shape), and Ch. 6.2 (the Randles-Ševčík equation).

Unless stated otherwise, quantities are SI: area in m^2, concentration in
mol/m^3 (numerically equal to mmol/L), diffusion coefficient in m^2/s,
current in A.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import FARADAY, STANDARD_TEMPERATURE, R

__all__ = [
    "cottrell_current",
    "ilkovic_diffusion_current",
    "polarographic_wave_current",
    "randles_sevcik_peak_current",
]

#: Dimensionless peak value :math:`\sqrt{\pi}\,\chi(\sigma t)_{max}` of the
#: reversible linear-sweep voltammogram (Bard & Faulkner, 2nd ed., eq. 6.2.18).
_RANDLES_SEVCIK_COEFFICIENT = 0.4463


def cottrell_current(t, n: int, A: float, C: float, D: float, F: float = FARADAY):
    r"""The Cottrell equation: :math:`i(t) = nFAC\sqrt{D/(\pi t)}`.

    Diffusion-limited current at a planar electrode after a potential step
    that drives the surface concentration of the electroactive species to
    zero (Bard & Faulkner, *Electrochemical Methods*, 2nd ed., eq. 5.2.11).

    Parameters
    ----------
    t : float or array-like of float
        Time since the potential step, in s (must be positive).
    n : int
        Electrons transferred per molecule.
    A : float
        Electrode area, in m^2.
    C : float
        Bulk concentration, in mol/m^3.
    D : float
        Diffusion coefficient, in m^2/s.
    F : float, default :data:`chemistrykit.constants.FARADAY`

    Returns
    -------
    float or ndarray
        Current, in A.

    Examples
    --------
    The Cottrell signature -- :math:`i\sqrt{t}` is constant in time:

    >>> a = cottrell_current(1.0, n=1, A=1e-6, C=1.0, D=1e-9)
    >>> b = cottrell_current(4.0, n=1, A=1e-6, C=1.0, D=1e-9)
    >>> round(a / b, 12)
    2.0
    """
    t = np.asarray(t, dtype=np.float64)
    result = n * F * A * C * np.sqrt(D / (np.pi * t))
    return float(result) if result.ndim == 0 else result


def ilkovic_diffusion_current(n: int, D: float, m: float, t_drop: float, C: float, average: bool = False) -> float:
    r"""The Ilkovič equation for the polarographic diffusion current at a dropping mercury electrode.

    :math:`i_d = k\,n D^{1/2} m^{2/3} t^{1/6} C` with :math:`k=708` for the
    maximum current at the end of each drop's life and :math:`k=607` for
    the drop-averaged current (Bard & Faulkner, *Electrochemical Methods*,
    2nd ed., eqs. 7.1.6-7.1.7). **Traditional polarographic units** are
    used here, not SI, because the numerical constants are defined in them.

    Parameters
    ----------
    n : int
        Electrons transferred per molecule.
    D : float
        Diffusion coefficient, in cm^2/s.
    m : float
        Mercury flow rate, in mg/s.
    t_drop : float
        Drop time, in s.
    C : float
        Bulk concentration, in mmol/L.
    average : bool, default False
        Return the drop-averaged (607) rather than the maximum (708) current.

    Returns
    -------
    float
        Diffusion current, in microamperes.

    Examples
    --------
    >>> round(ilkovic_diffusion_current(n=2, D=1e-5, m=2.0, t_drop=4.0, C=1.0), 2)
    8.96

    The diffusion current is strictly proportional to concentration -- the
    basis of quantitative polarography:

    >>> round(ilkovic_diffusion_current(2, 1e-5, 2.0, 4.0, C=3.0) / ilkovic_diffusion_current(2, 1e-5, 2.0, 4.0, C=1.0), 12)
    3.0
    """
    k = 607.0 if average else 708.0
    return float(k * n * np.sqrt(D) * m ** (2.0 / 3.0) * t_drop ** (1.0 / 6.0) * C)


def polarographic_wave_current(E, E_half: float, i_d: float, n: int, T: float = STANDARD_TEMPERATURE, R_gas: float = R, F: float = FARADAY):
    r"""Heyrovský-Ilkovič reversible wave: :math:`E = E_{1/2} + \frac{RT}{nF}\ln\frac{i_d-i}{i}`.

    Solved for the (cathodic) current, :math:`i = i_d/[1+\exp(nF(E-E_{1/2})/RT)]`:
    a sigmoidal wave rising from 0 to the diffusion plateau :math:`i_d` as
    the potential is swept negative through the half-wave potential
    (Bard & Faulkner, *Electrochemical Methods*, 2nd ed., eq. 5.4.22).

    Parameters
    ----------
    E : float or array-like of float
        Electrode potential, in V.
    E_half : float
        Half-wave potential, in V (characteristic of the species).
    i_d : float
        Limiting diffusion current (any unit; the result has the same unit).
    n : int
        Electrons transferred.
    T : float, default :data:`chemistrykit.constants.STANDARD_TEMPERATURE`
    R_gas : float, default :data:`chemistrykit.constants.R`
    F : float, default :data:`chemistrykit.constants.FARADAY`

    Returns
    -------
    float or ndarray
        Current, in the unit of `i_d`.

    Examples
    --------
    At the half-wave potential the current is exactly half the plateau:

    >>> polarographic_wave_current(-0.40, E_half=-0.40, i_d=8.0, n=2)
    4.0
    """
    E = np.asarray(E, dtype=np.float64)
    result = i_d / (1.0 + np.exp(n * F * (E - E_half) / (R_gas * T)))
    return float(result) if result.ndim == 0 else result


def randles_sevcik_peak_current(v, n: int, A: float, C: float, D: float, T: float = STANDARD_TEMPERATURE, R_gas: float = R, F: float = FARADAY):
    r"""The Randles-Ševčík equation: :math:`i_p = 0.4463\,nFAC\sqrt{nFvD/(RT)}`.

    Peak current of a reversible, diffusion-controlled linear-sweep or
    cyclic voltammogram at a planar electrode (Bard & Faulkner,
    *Electrochemical Methods*, 2nd ed., eq. 6.2.19). At 25 degC this is
    the familiar :math:`i_p = (2.69\times10^5)\,n^{3/2}AD^{1/2}Cv^{1/2}`
    in cm/mol/s units.

    Parameters
    ----------
    v : float or array-like of float
        Scan rate, in V/s.
    n : int
        Electrons transferred.
    A : float
        Electrode area, in m^2.
    C : float
        Bulk concentration, in mol/m^3.
    D : float
        Diffusion coefficient, in m^2/s.
    T : float, default :data:`chemistrykit.constants.STANDARD_TEMPERATURE`
    R_gas : float, default :data:`chemistrykit.constants.R`
    F : float, default :data:`chemistrykit.constants.FARADAY`

    Returns
    -------
    float or ndarray
        Peak current, in A.

    Examples
    --------
    1 mM ferrocene-like couple (n=1, D=1e-9 m^2/s) at a 3 mm diameter disc,
    100 mV/s:

    >>> import math
    >>> A = math.pi * (1.5e-3) ** 2
    >>> round(randles_sevcik_peak_current(0.1, n=1, A=A, C=1.0, D=1e-9) * 1e6, 1)  # microamps
    19.0

    Peak current grows as the square root of scan rate:

    >>> round(randles_sevcik_peak_current(0.4, 1, A, 1.0, 1e-9) / randles_sevcik_peak_current(0.1, 1, A, 1.0, 1e-9), 12)
    2.0
    """
    v = np.asarray(v, dtype=np.float64)
    result = _RANDLES_SEVCIK_COEFFICIENT * n * F * A * C * np.sqrt(n * F * v * D / (R_gas * T))
    return float(result) if result.ndim == 0 else result
