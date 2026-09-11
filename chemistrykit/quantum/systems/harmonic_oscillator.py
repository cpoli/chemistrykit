r"""The quantum harmonic oscillator, and its anharmonic (Morse) generalization.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 7.4 for the
harmonic-oscillator energy levels and Ch. 8.1 (real molecular vibrations)
for the Morse potential's exact vibrational eigenvalues, or Levine,
*Quantum Chemistry*, 7th ed., Ch. 2.5 and 4.3.
"""

from __future__ import annotations

import numpy as np
from scipy.special import factorial, hermite

from chemistrykit.constants import HBAR
from chemistrykit.quantum.core.base_system import QuantumSystem

__all__ = ["QuantumHarmonicOscillator", "MorseOscillator", "compare_harmonic_vs_morse"]


class QuantumHarmonicOscillator(QuantumSystem):
    r"""A 1D quantum harmonic oscillator, potential :math:`V(x)=\frac12kx^2`.

    .. math::

        \omega = \sqrt{k/m}, \qquad E_v = \hbar\omega\left(v+\frac12\right), \qquad v=0,1,2,\dots

    (Atkins & de Paula, *Physical Chemistry*, 11th ed., eq. 7.16-7.17.)
    Energy levels are exactly evenly spaced by :math:`\hbar\omega` -- the
    defining feature that :class:`MorseOscillator` below (a real bond's
    anharmonic potential) violates.

    Parameters
    ----------
    mass : float
        Oscillator mass (the reduced mass, for a two-body vibration), in kg.
    force_constant : float
        Force constant `k`, in N/m.

    Examples
    --------
    Consecutive levels are always separated by exactly :math:`\hbar\omega`:

    >>> ho = QuantumHarmonicOscillator(mass=1.6e-27, force_constant=500.0)
    >>> spacing1 = ho.energy(1) - ho.energy(0)
    >>> spacing2 = ho.energy(5) - ho.energy(4)
    >>> bool(round(float(spacing1), 30) == round(float(spacing2), 30))
    True
    >>> round(float(spacing1 / (HBAR * ho.angular_frequency)), 9)
    1.0
    """

    def __init__(self, mass: float, force_constant: float):
        if mass <= 0:
            raise ValueError("mass must be positive")
        if force_constant <= 0:
            raise ValueError("force_constant must be positive")
        self.mass = float(mass)
        self.force_constant = float(force_constant)

    @property
    def angular_frequency(self) -> float:
        r"""float: The classical angular frequency :math:`\omega=\sqrt{k/m}`, in rad/s."""
        return np.sqrt(self.force_constant / self.mass)

    @property
    def zero_point_energy(self) -> float:
        r"""float: The zero-point energy :math:`E_0=\frac12\hbar\omega`, in J."""
        return 0.5 * HBAR * self.angular_frequency

    def energy(self, v):
        """Return :math:`E_v=\\hbar\\omega(v+1/2)`.

        Parameters
        ----------
        v : int or array-like of int
            Vibrational quantum number(s), each >= 0.

        Returns
        -------
        float or ndarray
            Energy, in J.
        """
        v = np.asarray(v, dtype=np.float64)
        if np.any(v < 0):
            raise ValueError("v must be >= 0")
        return HBAR * self.angular_frequency * (v + 0.5)

    def wavefunction(self, x, v: int):
        r"""Return the harmonic-oscillator eigenfunction :math:`\psi_v(x)`.

        .. math::

            \psi_v(x) = \left(\frac{m\omega}{\pi\hbar}\right)^{1/4}
                        \frac{1}{\sqrt{2^vv!}}H_v(\xi)e^{-\xi^2/2},
                        \qquad \xi=x\sqrt{m\omega/\hbar}

        where :math:`H_v` is the physicists' Hermite polynomial (Atkins &
        de Paula, *Physical Chemistry*, 11th ed., eq. 7.18; here evaluated
        via :func:`scipy.special.hermite`).

        Parameters
        ----------
        x : float or array-like of float
            Displacement from equilibrium, in m.
        v : int
            Vibrational quantum number, >= 0.

        Returns
        -------
        float or ndarray
            Amplitude, in m^-1/2.

        Examples
        --------
        The ground state is a normalized Gaussian, maximal at ``x=0``:

        >>> ho = QuantumHarmonicOscillator(mass=1.6e-27, force_constant=500.0)
        >>> psi0_at_0 = ho.wavefunction(0.0, v=0)
        >>> psi0_away = ho.wavefunction(1.0e-11, v=0)
        >>> bool(psi0_at_0 > psi0_away > 0)
        True
        """
        if v < 0:
            raise ValueError("v must be >= 0")
        x = np.asarray(x, dtype=np.float64)
        beta = np.sqrt(self.mass * self.angular_frequency / HBAR)
        xi = x * beta
        prefactor = (self.mass * self.angular_frequency / (np.pi * HBAR)) ** 0.25 / np.sqrt(2.0**v * factorial(v))
        return prefactor * hermite(v)(xi) * np.exp(-0.5 * xi**2)


class MorseOscillator(QuantumSystem):
    r"""A diatomic-like vibration in the Morse potential, the standard anharmonic model of a real chemical bond.

    .. math::

        V(x) = D_e\left(1-e^{-ax}\right)^2, \qquad
        E_v = \hbar\omega\left(v+\frac12\right)
              - \frac{(\hbar\omega)^2}{4D_e}\left(v+\frac12\right)^2

    (P. M. Morse, *Phys. Rev.* 34, 57 (1929); Atkins & de Paula, *Physical
    Chemistry*, 11th ed., eq. 8.9, giving the *exact* vibrational
    eigenvalues of the Morse potential -- an important contrast with
    :class:`QuantumHarmonicOscillator`, whose energies are exact only for
    the harmonic (quadratic) approximation to the true potential, not for
    any real bond.) The bond dissociates once the levels stop being bound
    -- :attr:`v_max` gives the last bound vibrational level, beyond which
    this closed form no longer applies (Atkins & de Paula, *Physical
    Chemistry*, 11th ed., Ch. 8.1).

    `omega` is fixed by the curvature at the bottom of the well, exactly
    as for the harmonic oscillator, so the two are directly comparable at
    a shared force constant (see :func:`compare_harmonic_vs_morse`).

    Parameters
    ----------
    mass : float
        Reduced mass, in kg.
    force_constant : float
        Force constant at the potential minimum, `k`, in N/m (fixes
        :math:`\omega=\sqrt{k/m}`, identical to the harmonic oscillator's).
    dissociation_energy : float
        Well depth :math:`D_e` (measured from the potential minimum, *not*
        from the zero-point level), in J.

    Examples
    --------
    Level spacing decreases monotonically with `v` (anharmonicity), unlike
    the harmonic oscillator's perfectly even spacing:

    >>> morse = MorseOscillator(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19)
    >>> spacing = [morse.energy(v + 1) - morse.energy(v) for v in range(4)]
    >>> bool(all(spacing[i] > spacing[i + 1] for i in range(3)))
    True
    """

    def __init__(self, mass: float, force_constant: float, dissociation_energy: float):
        if mass <= 0:
            raise ValueError("mass must be positive")
        if force_constant <= 0:
            raise ValueError("force_constant must be positive")
        if dissociation_energy <= 0:
            raise ValueError("dissociation_energy must be positive")
        self.mass = float(mass)
        self.force_constant = float(force_constant)
        self.dissociation_energy = float(dissociation_energy)

    @property
    def angular_frequency(self) -> float:
        r"""float: :math:`\omega=\sqrt{k/m}`, in rad/s -- same formula as the harmonic oscillator."""
        return np.sqrt(self.force_constant / self.mass)

    @property
    def anharmonicity_constant(self) -> float:
        r"""float: The dimensionless anharmonicity constant :math:`x_e=\hbar\omega/(4D_e)`."""
        return HBAR * self.angular_frequency / (4.0 * self.dissociation_energy)

    @property
    def v_max(self) -> int:
        r"""int: The highest bound vibrational quantum number.

        Found from :math:`dE_v/dv=0` at :math:`v=v_{max}+\frac12`, i.e.
        the last integer `v` before the Morse levels turn over and start
        decreasing (an artifact of the quadratic-in-v formula beyond the
        true dissociation limit; Atkins & de Paula, *Physical Chemistry*,
        11th ed., Ch. 8.1).
        """
        return int(np.floor(1.0 / (2.0 * self.anharmonicity_constant) - 0.5))

    def energy(self, v):
        r"""Return the exact Morse vibrational energy :math:`E_v`.

        Parameters
        ----------
        v : int or array-like of int
            Vibrational quantum number(s), each >= 0 (and, physically,
            ``<= v_max``, though this is not enforced).

        Returns
        -------
        float or ndarray
            Energy, in J, measured from the bottom of the well.
        """
        v = np.asarray(v, dtype=np.float64)
        if np.any(v < 0):
            raise ValueError("v must be >= 0")
        hw = HBAR * self.angular_frequency
        return hw * (v + 0.5) - hw**2 / (4.0 * self.dissociation_energy) * (v + 0.5) ** 2


def compare_harmonic_vs_morse(mass: float, force_constant: float, dissociation_energy: float, v_max: int):
    r"""Tabulate harmonic-oscillator vs. Morse-potential vibrational energies at a shared force constant.

    Both models share the same curvature at the well minimum (same
    `mass`, `force_constant`, hence the same :math:`\omega`), so they
    agree closely for the lowest levels (:math:`v\approx0`, where the
    Morse potential is well approximated by its harmonic term) and
    diverge increasingly as `v` grows: the Morse level spacing shrinks
    toward zero as :math:`v\to v_{max}` (approaching dissociation, where
    the vibrational levels become a near-continuum), while the harmonic
    spacing stays exactly :math:`\hbar\omega` at every `v` -- the harmonic
    approximation's defining failure mode for a real bond (Atkins & de
    Paula, *Physical Chemistry*, 11th ed., Ch. 8.1).

    Parameters
    ----------
    mass : float
        Reduced mass, in kg.
    force_constant : float
        Shared force constant `k`, in N/m.
    dissociation_energy : float
        Morse well depth :math:`D_e`, in J.
    v_max : int
        Highest vibrational quantum number to tabulate (may safely exceed
        the Morse oscillator's own bound-level cutoff for comparison
        purposes -- the *harmonic* levels stay well-defined at any `v`).

    Returns
    -------
    v : ndarray, shape (v_max + 1,)
        Vibrational quantum numbers ``0, 1, ..., v_max``.
    harmonic_energies : ndarray, shape (v_max + 1,)
        Harmonic-oscillator energies, in J.
    morse_energies : ndarray, shape (v_max + 1,)
        Morse-potential energies, in J.

    Examples
    --------
    The two models agree closely at `v=0` (small relative difference)
    and disagree much more by `v=5` (anharmonicity accumulates):

    >>> v, E_harmonic, E_morse = compare_harmonic_vs_morse(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19, v_max=5)
    >>> relative_diff_0 = abs(E_harmonic[0] - E_morse[0]) / E_harmonic[0]
    >>> relative_diff_5 = abs(E_harmonic[5] - E_morse[5]) / E_harmonic[5]
    >>> bool(relative_diff_5 > relative_diff_0)
    True
    """
    harmonic = QuantumHarmonicOscillator(mass, force_constant)
    morse = MorseOscillator(mass, force_constant, dissociation_energy)
    v = np.arange(0, v_max + 1)
    return v, harmonic.energy(v), morse.energy(v)
