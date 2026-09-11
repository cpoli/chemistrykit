r"""The Maxwell-Boltzmann speed distribution.

.. math::

    f(v) = 4\pi\left(\frac{m}{2\pi k_BT}\right)^{3/2}v^2\,e^{-mv^2/(2k_BT)}

the probability density of molecular *speeds* (as opposed to velocity
*components*, which are each independently Gaussian) in a classical
ideal gas at equilibrium -- J. C. Maxwell, *Philos. Mag.* 19, 19 (1860);
L. Boltzmann, *Wiener Berichte* 66, 275 (1872); McQuarrie, *Statistical
Mechanics*, Ch. 27; Atkins & de Paula, *Physical Chemistry*, 11th ed.,
Ch. 1.2b.

A molecular-dynamics trajectory's velocity ensemble is a direct
cross-check for this distribution: at equilibrium, an
:class:`chemistrykit.md.systems.lj_fluid.LJFluid` run's particle speeds
should follow :class:`MaxwellBoltzmannSpeedDistribution` at the
simulation's instantaneous temperature (see
``examples/md/lj_fluid/plot_02_maxwell_boltzmann_check.py``).
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf

from chemistrykit.constants import K_B

__all__ = [
    "MaxwellBoltzmannSpeedDistribution",
    "most_probable_speed",
    "mean_speed",
    "rms_speed",
]


def most_probable_speed(mass, temperature, k_b=K_B):
    r"""Most probable speed, :math:`v_p=\sqrt{2k_BT/m}` -- where :math:`f(v)` peaks.

    Parameters
    ----------
    mass : float
    temperature : float
    k_b : float, default :data:`chemistrykit.constants.K_B`

    Returns
    -------
    float
    """
    return np.sqrt(2.0 * k_b * temperature / mass)


def mean_speed(mass, temperature, k_b=K_B):
    r"""Mean speed, :math:`\bar v=\sqrt{8k_BT/(\pi m)}`.

    Parameters
    ----------
    mass : float
    temperature : float
    k_b : float, default :data:`chemistrykit.constants.K_B`

    Returns
    -------
    float
    """
    return np.sqrt(8.0 * k_b * temperature / (np.pi * mass))


def rms_speed(mass, temperature, k_b=K_B):
    r"""Root-mean-square speed, :math:`v_{\text{rms}}=\sqrt{3k_BT/m}`.

    Parameters
    ----------
    mass : float
    temperature : float
    k_b : float, default :data:`chemistrykit.constants.K_B`

    Returns
    -------
    float
    """
    return np.sqrt(3.0 * k_b * temperature / mass)


class MaxwellBoltzmannSpeedDistribution:
    r"""The Maxwell-Boltzmann speed distribution for a gas of mass `mass` at `temperature`.

    Parameters
    ----------
    mass : float
        Molecular mass (in kg for SI/real units; in reduced units, e.g.
        those of :class:`chemistrykit.md.systems.lj_fluid.LJFluid`, pass
        the reduced mass and ``k_b=1.0``).
    temperature : float
        Absolute temperature.
    k_b : float, default :data:`chemistrykit.constants.K_B`
        Boltzmann constant, in the same unit system as `mass`/`temperature`.

    Examples
    --------
    The three characteristic speeds stand in the fixed ratio
    :math:`v_p:\bar v:v_{\text{rms}} = \sqrt2:\sqrt{8/\pi}:\sqrt3`,
    independent of mass or temperature -- a purely geometric property of
    the distribution's shape (Atkins & de Paula, *Physical Chemistry*,
    11th ed., Ch. 1.2b):

    >>> import numpy as np
    >>> dist = MaxwellBoltzmannSpeedDistribution(mass=6.63e-26, temperature=298.15)
    >>> vp, vbar, vrms = dist.most_probable_speed(), dist.mean_speed(), dist.rms_speed()
    >>> np.allclose([vbar / vp, vrms / vp], [np.sqrt(8.0 / np.pi) / np.sqrt(2.0), np.sqrt(3.0) / np.sqrt(2.0)])
    True

    The pdf integrates to 1 over all speeds:

    >>> from scipy.integrate import quad
    >>> total, _ = quad(dist.pdf, 0.0, np.inf)
    >>> round(total, 6)
    1.0
    """

    def __init__(self, mass: float, temperature: float, k_b: float = K_B):
        if mass <= 0 or temperature <= 0:
            raise ValueError("mass and temperature must be positive")
        self.mass = float(mass)
        self.temperature = float(temperature)
        self.k_b = float(k_b)

    def pdf(self, v):
        """Probability density f(v) at speed(s) `v`.

        Parameters
        ----------
        v : float or array-like of float
            Must be non-negative.

        Returns
        -------
        float or ndarray
        """
        v = np.asarray(v, dtype=np.float64)
        prefactor = 4.0 * np.pi * (self.mass / (2.0 * np.pi * self.k_b * self.temperature)) ** 1.5
        return prefactor * v**2 * np.exp(-self.mass * v**2 / (2.0 * self.k_b * self.temperature))

    def cdf(self, v):
        r"""Cumulative probability :math:`P(V\leq v)`.

        .. math::

            F(v) = \mathrm{erf}(av) - \frac{2}{\sqrt\pi}av\,e^{-(av)^2},
            \qquad a=\sqrt{\frac{m}{2k_BT}}

        Parameters
        ----------
        v : float or array-like of float

        Returns
        -------
        float or ndarray

        Examples
        --------
        >>> dist = MaxwellBoltzmannSpeedDistribution(mass=6.63e-26, temperature=298.15)
        >>> float(dist.cdf(0.0))
        0.0
        """
        v = np.asarray(v, dtype=np.float64)
        a = np.sqrt(self.mass / (2.0 * self.k_b * self.temperature))
        return erf(a * v) - (2.0 / np.sqrt(np.pi)) * a * v * np.exp(-((a * v) ** 2))

    def most_probable_speed(self) -> float:
        """See :func:`most_probable_speed`."""
        return most_probable_speed(self.mass, self.temperature, self.k_b)

    def mean_speed(self) -> float:
        """See :func:`mean_speed`."""
        return mean_speed(self.mass, self.temperature, self.k_b)

    def rms_speed(self) -> float:
        """See :func:`rms_speed`."""
        return rms_speed(self.mass, self.temperature, self.k_b)

    def sample(self, n: int, rng=None):
        r"""Draw `n` random speeds from this distribution.

        Each of the 3 Cartesian velocity components of a classical ideal
        gas is independently Gaussian-distributed with variance
        :math:`k_BT/m`; the *speed* (the norm of the velocity vector) of
        such a vector is exactly Maxwell-Boltzmann distributed, which is
        the simplest way to sample it (no inverse-CDF / rejection-sampling
        machinery needed).

        Parameters
        ----------
        n : int
        rng : int, numpy.random.Generator, or None

        Returns
        -------
        ndarray, shape (n,)

        Examples
        --------
        >>> dist = MaxwellBoltzmannSpeedDistribution(mass=6.63e-26, temperature=298.15)
        >>> speeds = dist.sample(200_000, rng=0)
        >>> bool(abs(speeds.mean() - dist.mean_speed()) / dist.mean_speed() < 0.01)
        True
        """
        rng = np.random.default_rng(rng)
        sigma = np.sqrt(self.k_b * self.temperature / self.mass)
        components = rng.normal(0.0, sigma, size=(n, 3))
        return np.linalg.norm(components, axis=1)
