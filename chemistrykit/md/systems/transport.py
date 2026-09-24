r"""Transport coefficients from MD trajectories: Einstein and Green-Kubo diffusion.

Two exactly equivalent routes from a simulated trajectory to the
self-diffusion coefficient :math:`D` of a :math:`d`-dimensional fluid:

- Einstein (1905): the mean-squared displacement grows linearly at long
  times, :math:`\langle|\vec r(t)-\vec r(0)|^2\rangle \to 2dDt`
  (A. Einstein, *Ann. Phys.* 17, 549 (1905)).
- Green-Kubo: :math:`D` is the time integral of the velocity
  autocorrelation function,
  :math:`D = \frac{1}{d}\int_0^\infty \langle\vec v(0)\cdot\vec v(t)\rangle\,dt`
  (M. S. Green, *J. Chem. Phys.* 22, 398 (1954); R. Kubo, *J. Phys. Soc.
  Jpn.* 12, 570 (1957)).

See Allen & Tildesley, *Computer Simulation of Liquids*, 2nd ed., Ch. 2.7
and 8.2. All averages below run over particles and over every available
time origin.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

__all__ = [
    "unwrap_trajectory",
    "mean_squared_displacement",
    "velocity_autocorrelation",
    "einstein_diffusion_coefficient",
    "green_kubo_diffusion_coefficient",
]


def unwrap_trajectory(positions, box_length: Optional[float]):
    """Undo periodic wrapping of a sampled trajectory.

    Accumulates the minimum-image displacement between consecutive frames,
    so it is exact provided no particle moves more than half a box length
    between two samples.

    Parameters
    ----------
    positions : array-like, shape (n_frames, n_particles, n_dim)
    box_length : float or None
        Cubic periodic box side length; ``None`` or ``<= 0`` returns the
        positions unchanged.

    Returns
    -------
    ndarray, shape (n_frames, n_particles, n_dim)

    Examples
    --------
    >>> import numpy as np
    >>> wrapped = np.array([[[9.5]], [[0.3]], [[1.1]]])  # crosses the box edge at 10
    >>> unwrap_trajectory(wrapped, box_length=10.0).ravel()
    array([ 9.5, 10.3, 11.1])
    """
    positions = np.asarray(positions, dtype=np.float64)
    if box_length is None or box_length <= 0:
        return positions.copy()
    steps = np.diff(positions, axis=0)
    steps -= box_length * np.round(steps / box_length)
    unwrapped = np.empty_like(positions)
    unwrapped[0] = positions[0]
    unwrapped[1:] = positions[0] + np.cumsum(steps, axis=0)
    return unwrapped


def mean_squared_displacement(positions, box_length: Optional[float] = None, max_lag: Optional[int] = None):
    r"""Time-origin-averaged mean-squared displacement.

    .. math::

        \mathrm{MSD}(\ell) = \left\langle |\vec r_i(t_0+\ell) - \vec r_i(t_0)|^2 \right\rangle_{i,\,t_0}

    Parameters
    ----------
    positions : array-like, shape (n_frames, n_particles, n_dim)
        Sampled positions (e.g. :attr:`chemistrykit.md.MDResult.positions`).
        Wrapped periodic positions are unwrapped first via
        :func:`unwrap_trajectory` when `box_length` is given.
    box_length : float or None, optional
    max_lag : int, optional
        Largest lag, in frames; defaults to ``n_frames // 2``.

    Returns
    -------
    ndarray, shape (max_lag + 1,)
        MSD at lags ``0, 1, ..., max_lag`` frames.

    Examples
    --------
    Particles moving ballistically at unit speed have MSD equal to the lag squared:

    >>> import numpy as np
    >>> t = np.arange(5.0)
    >>> positions = np.stack([t, np.zeros(5), np.zeros(5)], axis=-1)[:, None, :]
    >>> mean_squared_displacement(positions, max_lag=3)
    array([0., 1., 4., 9.])
    """
    positions = unwrap_trajectory(positions, box_length)
    n_frames = positions.shape[0]
    max_lag = n_frames // 2 if max_lag is None else int(max_lag)
    if not 0 <= max_lag < n_frames:
        raise ValueError("max_lag must be in [0, n_frames)")
    msd = np.empty(max_lag + 1)
    for lag in range(max_lag + 1):
        disp = positions[lag:] - positions[: n_frames - lag]
        msd[lag] = np.mean(np.sum(disp * disp, axis=-1))
    return msd


def velocity_autocorrelation(velocities, max_lag: Optional[int] = None, normalize: bool = False):
    r"""Time-origin-averaged velocity autocorrelation function.

    .. math::

        C(\ell) = \left\langle \vec v_i(t_0)\cdot\vec v_i(t_0+\ell) \right\rangle_{i,\,t_0}

    Parameters
    ----------
    velocities : array-like, shape (n_frames, n_particles, n_dim)
    max_lag : int, optional
        Largest lag, in frames; defaults to ``n_frames // 2``.
    normalize : bool, default False
        Divide by :math:`C(0)`.

    Returns
    -------
    ndarray, shape (max_lag + 1,)

    Examples
    --------
    >>> import numpy as np
    >>> v = np.array([[[1.0, 0.0]], [[0.0, 1.0]], [[-1.0, 0.0]]])  # a velocity rotating by 90 degrees per frame
    >>> velocity_autocorrelation(v, max_lag=2)
    array([ 1.,  0., -1.])
    """
    velocities = np.asarray(velocities, dtype=np.float64)
    n_frames = velocities.shape[0]
    max_lag = n_frames // 2 if max_lag is None else int(max_lag)
    if not 0 <= max_lag < n_frames:
        raise ValueError("max_lag must be in [0, n_frames)")
    c = np.empty(max_lag + 1)
    for lag in range(max_lag + 1):
        c[lag] = np.mean(np.sum(velocities[lag:] * velocities[: n_frames - lag], axis=-1))
    if normalize:
        c = c / c[0]
    return c


def einstein_diffusion_coefficient(t, msd, n_dim: int = 3, fit_from: float = 0.0) -> float:
    r"""Self-diffusion coefficient from the Einstein relation :math:`\mathrm{MSD}=2dDt`.

    Fits a straight line to ``msd`` against ``t`` for ``t >= fit_from``
    (skipping the short-time ballistic regime) and returns slope/(2d).

    Parameters
    ----------
    t : array-like
        Lag times.
    msd : array-like
        Mean-squared displacement at each lag (see :func:`mean_squared_displacement`).
    n_dim : int, default 3
    fit_from : float, default 0.0

    Returns
    -------
    float

    Examples
    --------
    >>> import numpy as np
    >>> t = np.linspace(0.0, 10.0, 11)
    >>> round(einstein_diffusion_coefficient(t, 6.0 * 0.25 * t + 0.1), 12)
    0.25
    """
    t = np.asarray(t, dtype=np.float64)
    msd = np.asarray(msd, dtype=np.float64)
    mask = t >= fit_from
    if mask.sum() < 2:
        raise ValueError("need at least two points with t >= fit_from")
    slope = np.polyfit(t[mask], msd[mask], 1)[0]
    return float(slope / (2.0 * n_dim))


def green_kubo_diffusion_coefficient(t, vacf, n_dim: int = 3) -> float:
    r"""Self-diffusion coefficient from the Green-Kubo integral of the VACF.

    .. math::

        D = \frac{1}{d}\int_0^{t_{\max}} \langle\vec v(0)\cdot\vec v(t)\rangle\,dt

    evaluated by the trapezoidal rule.

    Parameters
    ----------
    t : array-like
        Lag times, starting at 0.
    vacf : array-like
        Unnormalized velocity autocorrelation (see :func:`velocity_autocorrelation`).
    n_dim : int, default 3

    Returns
    -------
    float

    Examples
    --------
    An exponentially decaying VACF :math:`C(t)=3(k_BT/m)e^{-\gamma t}` (Langevin
    dynamics) integrates to the Stokes-Einstein-like value :math:`k_BT/(m\gamma)`:

    >>> import numpy as np
    >>> t = np.linspace(0.0, 40.0, 40001)
    >>> round(green_kubo_diffusion_coefficient(t, 3.0 * 1.5 * np.exp(-2.0 * t)), 6)
    0.75
    """
    t = np.asarray(t, dtype=np.float64)
    vacf = np.asarray(vacf, dtype=np.float64)
    integral = float(np.sum(0.5 * (vacf[1:] + vacf[:-1]) * np.diff(t)))
    return integral / n_dim
