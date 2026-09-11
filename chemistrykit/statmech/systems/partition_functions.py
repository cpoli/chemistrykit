r"""Translational, rotational, and vibrational partition functions.

See McQuarrie, *Statistical Mechanics* (1976), Ch. 6, or Atkins & de
Paula, *Physical Chemistry*, 11th ed., Ch. 13, for the derivation of each
partition function below and the thermodynamic-function formulas
(:math:`U`, :math:`S`, :math:`C_V`) built from them.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from chemistrykit.constants import K_B, NA, C, H, R
from chemistrykit.statmech.core.base_system import PartitionFunction, ThermodynamicFunctions

__all__ = [
    "TranslationalPartitionFunction",
    "RotationalPartitionFunctionLinear",
    "VibrationalPartitionFunctionHarmonic",
    "IdealGasMolecule",
    "sackur_tetrode_entropy",
]


class TranslationalPartitionFunction(PartitionFunction):
    r"""The translational partition function of a particle in a box (volume) `V`.

    .. math::

        q_{\text{trans}} = \left(\frac{2\pi mk_BT}{h^2}\right)^{3/2}V

    the high-temperature (particle-in-a-box energy levels essentially
    continuous) limit of the exact 3D particle-in-a-box partition
    function, valid whenever :math:`V^{1/3}\gg\Lambda` (the thermal de
    Broglie wavelength, :func:`chemistrykit.statmech.utils.thermal_wavelength.thermal_de_broglie_wavelength`)
    -- true for any real gas well above its condensation point (McQuarrie,
    *Statistical Mechanics*, Ch. 6.1). Because translational motion
    permutes identical molecules through the same physical states (unlike
    the rotational/vibrational modes below, which are internal to each
    molecule), the `N`-molecule thermodynamic functions here include the
    :math:`1/N!` indistinguishability correction, giving the
    Sackur-Tetrode entropy (see :meth:`entropy`, :func:`sackur_tetrode_entropy`).

    Parameters
    ----------
    mass : float
        Molecular mass, in kg.
    volume : float
        Container volume, in m^3.

    Examples
    --------
    Internal energy is exactly the equipartition value :math:`\frac32Nk_BT`
    (3 translational quadratic degrees of freedom) at every temperature --
    the classical limit is essentially exact here:

    >>> q = TranslationalPartitionFunction(mass=6.63e-26, volume=1.0e-3)
    >>> round(q.internal_energy(298.15, N=1.0) / (1.5 * 1.380649e-23 * 298.15), 6)
    1.0
    """

    def __init__(self, mass: float, volume: float):
        if mass <= 0 or volume <= 0:
            raise ValueError("mass and volume must be positive")
        self.mass = float(mass)
        self.volume = float(volume)

    def value(self, T):
        return (2.0 * np.pi * self.mass * K_B * T / H**2) ** 1.5 * self.volume

    def internal_energy(self, T, N: float = NA):
        return 1.5 * N * K_B * T

    def entropy(self, T, N: float = NA):
        r"""Return the Sackur-Tetrode translational entropy.

        .. math::

            S = Nk_B\left[\ln\left(\frac{q_{\text{trans}}}{N}\right) + \frac{5}{2}\right]

        (O. Sackur, *Ann. Phys.* 36, 958 (1911); H. Tetrode, *Ann.
        Phys.* 38, 434 (1912); Atkins & de Paula, *Physical Chemistry*,
        11th ed., eq. 13.24.) See also the standalone convenience
        function :func:`sackur_tetrode_entropy` for computing this
        directly from `P` and `T` without building a container volume by
        hand.
        """
        q = self.value(T)
        return N * K_B * (np.log(q / N) + 2.5)

    def heat_capacity_v(self, T, N: float = NA):
        return 1.5 * N * K_B


def sackur_tetrode_entropy(mass: float, T: float, P: float, R_gas: float = R):
    r"""Molar translational entropy of an ideal gas via the Sackur-Tetrode equation.

    A convenience wrapper around
    :meth:`TranslationalPartitionFunction.entropy` that takes the
    pressure directly, using the ideal-gas volume per molecule
    :math:`V/N=k_BT/P`.

    Parameters
    ----------
    mass : float
        Molecular mass, in kg.
    T : float
        Absolute temperature, in K.
    P : float
        Pressure, in Pa.
    R_gas : float, default :data:`chemistrykit.constants.R`

    Returns
    -------
    float
        Molar entropy, in J/(mol K).

    Examples
    --------
    The textbook value for argon gas at 298.15 K and 1 bar is
    :math:`S_m^\circ\approx154.8` J/(mol K) (Atkins & de Paula, *Physical
    Chemistry*, 11th ed., Table 13.1):

    >>> import scipy.constants as sc
    >>> S = sackur_tetrode_entropy(mass=39.948 * sc.atomic_mass, T=298.15, P=1.0e5)
    >>> round(float(S), 1)
    154.8
    """
    volume_per_molecule = K_B * T / P
    q_per_molecule = (2.0 * np.pi * mass * K_B * T / H**2) ** 1.5 * volume_per_molecule
    return R_gas * (np.log(q_per_molecule) + 2.5)


class RotationalPartitionFunctionLinear(PartitionFunction):
    r"""The rotational partition function of a linear rigid rotor, high-temperature (classical) limit.

    .. math::

        q_{\text{rot}} = \frac{T}{\sigma\Theta_{\text{rot}}}, \qquad
        \Theta_{\text{rot}} = \frac{h^2}{8\pi^2Ik_B}

    where `I` is the moment of inertia and :math:`\sigma` the rotational
    symmetry number (2 for a homonuclear diatomic, 1 for heteronuclear).
    This is the classical (high-temperature) approximation to the exact
    sum over quantized rotational levels, valid whenever
    :math:`T\gg\Theta_{\text{rot}}` -- true at room temperature for
    essentially every molecule except H2 and its isotopologues (McQuarrie,
    *Statistical Mechanics*, Ch. 6.3; Atkins & de Paula, *Physical
    Chemistry*, 11th ed., Ch. 13.2c). Because this classical limit is
    used, :meth:`internal_energy` and :meth:`heat_capacity_v` are the
    exact equipartition values (:math:`Nk_BT` and :math:`Nk_B`
    respectively, from 2 rotational quadratic degrees of freedom) rather
    than a sum over levels -- this is the approximation being flagged.

    Parameters
    ----------
    moment_of_inertia : float
        Moment of inertia, in kg m^2.
    symmetry_number : int, default 1
        Rotational symmetry number :math:`\sigma`.

    Examples
    --------
    >>> q = RotationalPartitionFunctionLinear(moment_of_inertia=1.45e-46, symmetry_number=1)  # ~HCl-like
    >>> round(q.heat_capacity_v(298.15, N=1.0) / 1.380649e-23, 6)
    1.0
    """

    def __init__(self, moment_of_inertia: float, symmetry_number: int = 1):
        if moment_of_inertia <= 0:
            raise ValueError("moment_of_inertia must be positive")
        if symmetry_number < 1:
            raise ValueError("symmetry_number must be a positive integer")
        self.moment_of_inertia = float(moment_of_inertia)
        self.symmetry_number = int(symmetry_number)

    @property
    def rotational_temperature(self) -> float:
        r"""The rotational temperature :math:`\Theta_{\text{rot}}=h^2/(8\pi^2Ik_B)`, in K."""
        return H**2 / (8.0 * np.pi**2 * self.moment_of_inertia * K_B)

    def value(self, T):
        return T / (self.symmetry_number * self.rotational_temperature)

    def internal_energy(self, T, N: float = NA):
        return N * K_B * T

    def entropy(self, T, N: float = NA):
        return N * K_B * (np.log(self.value(T)) + 1.0)

    def heat_capacity_v(self, T, N: float = NA):
        return N * K_B


class VibrationalPartitionFunctionHarmonic(PartitionFunction):
    r"""The vibrational partition function of a single harmonic-oscillator mode.

    .. math::

        q_{\text{vib}} = \frac{1}{1-e^{-\Theta_{\text{vib}}/T}}, \qquad
        \Theta_{\text{vib}} = \frac{h\nu}{k_B}

    with energies measured from the :math:`v=0` ground vibrational level
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 13.2d).
    Unlike :class:`RotationalPartitionFunctionLinear`'s classical limit,
    this is an *exact* closed form within the harmonic approximation,
    valid at every temperature -- it is the harmonic-potential assumption
    itself (rather than a classical/high-`T` truncation) that is
    approximate, since real bonds are anharmonic (see
    :class:`chemistrykit.md.systems.pair_potentials.Morse`).
    :meth:`heat_capacity_v` is the Einstein-solid heat-capacity formula
    (A. Einstein, *Ann. Phys.* 22, 180 (1907)), reused here for a single
    vibrational mode rather than a 3D solid's lattice vibrations.

    Parameters
    ----------
    frequency : float
        Vibrational frequency, in Hz.

    Examples
    --------
    The Einstein/vibrational heat capacity interpolates between 0 (T=0,
    the mode is frozen out) and the classical equipartition value
    :math:`Nk_B` (:math:`T\to\infty`, one quadratic kinetic + one
    quadratic potential term):

    >>> q = VibrationalPartitionFunctionHarmonic(frequency=8.7e13)  # ~HCl-like, ~2900 cm^-1
    >>> round(float(q.heat_capacity_v(1.0, N=1.0)), 30)
    0.0
    >>> ratio = q.heat_capacity_v(1.0e6, N=1.0) / 1.380649e-23
    >>> bool(0.99 < ratio < 1.0)
    True
    """

    def __init__(self, frequency: float):
        if frequency <= 0:
            raise ValueError("frequency must be positive")
        self.frequency = float(frequency)

    @classmethod
    def from_wavenumber(cls, wavenumber_cm_inv: float) -> VibrationalPartitionFunctionHarmonic:
        r"""Build from a vibrational wavenumber in cm^-1 via :math:`\nu=c\tilde\nu`.

        Parameters
        ----------
        wavenumber_cm_inv : float
            Wavenumber :math:`\tilde\nu`, in cm^-1 (the usual IR-spectroscopy unit).

        Returns
        -------
        VibrationalPartitionFunctionHarmonic

        Examples
        --------
        >>> q = VibrationalPartitionFunctionHarmonic.from_wavenumber(2886.0)  # HCl fundamental
        >>> round(q.frequency / 1.0e13, 3)
        8.652
        """
        frequency = C * wavenumber_cm_inv * 100.0
        return cls(frequency=frequency)

    @property
    def vibrational_temperature(self) -> float:
        r"""The vibrational temperature :math:`\Theta_{\text{vib}}=h\nu/k_B`, in K."""
        return H * self.frequency / K_B

    def value(self, T):
        x = self.vibrational_temperature / T
        return 1.0 / (1.0 - np.exp(-x))

    def internal_energy(self, T, N: float = NA):
        x = self.vibrational_temperature / T
        return N * K_B * self.vibrational_temperature * np.exp(-x) / (1.0 - np.exp(-x))

    def entropy(self, T, N: float = NA):
        x = self.vibrational_temperature / T
        return N * K_B * (x / (np.exp(x) - 1.0) - np.log(1.0 - np.exp(-x)))

    def heat_capacity_v(self, T, N: float = NA):
        x = self.vibrational_temperature / T
        return N * K_B * x**2 * np.exp(-x) / (1.0 - np.exp(-x)) ** 2


class IdealGasMolecule:
    r"""Combined translational + rotational + vibrational thermodynamic functions for one ideal-gas species.

    Sums the independent contributions of a
    :class:`TranslationalPartitionFunction`, an optional
    :class:`RotationalPartitionFunctionLinear`, and zero or more
    :class:`VibrationalPartitionFunctionHarmonic` modes -- valid because
    the total molecular partition function factorizes,
    :math:`q=q_{\text{trans}}q_{\text{rot}}q_{\text{vib},1}q_{\text{vib},2}\cdots`,
    whenever these degrees of freedom are independent (the standard
    approximation for a rigid-rotor/harmonic-oscillator molecule; Atkins
    & de Paula, *Physical Chemistry*, 11th ed., Ch. 13.1).

    Parameters
    ----------
    mass : float
        Molecular mass, in kg.
    volume : float
        Container volume, in m^3.
    moment_of_inertia : float, optional
        Moment of inertia for a linear-rotor rotational contribution; omit
        for an atom (no rotational degrees of freedom).
    symmetry_number : int, default 1
    vibrational_frequencies : sequence of float, optional
        Vibrational mode frequencies, in Hz.

    Examples
    --------
    At very high temperature every mode approaches its classical
    equipartition heat capacity: :math:`\frac32R` (translation) +
    :math:`R` (linear rotation) + :math:`R` per vibrational mode:

    >>> molecule = IdealGasMolecule(mass=6.63e-26, volume=1.0e-3, moment_of_inertia=1.45e-46, vibrational_frequencies=[8.7e13])
    >>> Cv_over_R = molecule.heat_capacity_v(1.0e7) / 8.31446261815324
    >>> round(float(Cv_over_R), 2)
    3.5
    """

    def __init__(
        self,
        mass: float,
        volume: float,
        moment_of_inertia: float = None,
        symmetry_number: int = 1,
        vibrational_frequencies: Sequence[float] = (),
    ):
        self.translational = TranslationalPartitionFunction(mass, volume)
        self.rotational = RotationalPartitionFunctionLinear(moment_of_inertia, symmetry_number) if moment_of_inertia is not None else None
        self.vibrational_modes = [VibrationalPartitionFunctionHarmonic(f) for f in vibrational_frequencies]

    def _modes(self):
        modes = [self.translational]
        if self.rotational is not None:
            modes.append(self.rotational)
        modes.extend(self.vibrational_modes)
        return modes

    def internal_energy(self, T, N: float = NA):
        """Total internal energy: sum of each mode's contribution. See :meth:`chemistrykit.statmech.core.base_system.PartitionFunction.internal_energy`."""
        return sum(mode.internal_energy(T, N) for mode in self._modes())

    def entropy(self, T, N: float = NA):
        """Total entropy: sum of each mode's contribution."""
        return sum(mode.entropy(T, N) for mode in self._modes())

    def heat_capacity_v(self, T, N: float = NA):
        """Total heat capacity: sum of each mode's contribution."""
        return sum(mode.heat_capacity_v(T, N) for mode in self._modes())

    def thermodynamic_functions(self, T, N: float = NA) -> ThermodynamicFunctions:
        """Bundle q, U, S, Cv, and A at temperature `T` into a :class:`~chemistrykit.statmech.core.base_system.ThermodynamicFunctions`.

        Parameters
        ----------
        T : float
        N : float, default :data:`chemistrykit.constants.NA`

        Returns
        -------
        chemistrykit.statmech.core.base_system.ThermodynamicFunctions
        """
        q = float(np.prod([mode.value(T) for mode in self._modes()]))
        U = self.internal_energy(T, N)
        S = self.entropy(T, N)
        Cv = self.heat_capacity_v(T, N)
        return ThermodynamicFunctions(T=T, q=q, U=U, S=S, Cv=Cv, A=U - T * S)
