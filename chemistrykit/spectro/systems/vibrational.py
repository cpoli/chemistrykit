r"""Vibrational (IR) band positions: harmonic vs. Morse (anharmonic) fundamentals/overtones, and triatomic normal modes.

Builds on :class:`chemistrykit.quantum.systems.harmonic_oscillator.QuantumHarmonicOscillator`
and :class:`~chemistrykit.quantum.systems.harmonic_oscillator.MorseOscillator`
for the underlying vibrational energy levels -- this module only converts
level *differences* to the IR-observable wavenumbers and adds the
diatomic-to-triatomic generalization (normal-mode analysis via the
Wilson GF-matrix method). See Atkins & de Paula, *Physical Chemistry*,
11th ed., Ch. 12.4-12.5 (vibrational spectroscopy) and Wilson, Decius &
Cross, *Molecular Vibrations: The Theory of Infrared and Raman
Vibrational Spectra* (1955), Ch. 4 (internal coordinates and the GF
method) throughout.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla

from chemistrykit.constants import C
from chemistrykit.quantum.systems.harmonic_oscillator import MorseOscillator, QuantumHarmonicOscillator
from chemistrykit.quantum.utils.secular_equation import solve_secular_equation
from chemistrykit.spectro.systems.rotational import energy_to_wavenumber

__all__ = [
    "harmonic_fundamental_wavenumber",
    "morse_transition_wavenumbers",
    "anharmonicity_from_overtones",
    "TriatomicNormalModes",
    "NormalModeResult",
]


def harmonic_fundamental_wavenumber(force_constant: float, reduced_mass: float) -> float:
    r"""The harmonic-oscillator fundamental (:math:`v=0\to1`) IR band position.

    .. math::

        \tilde\nu = \frac{1}{2\pi c}\sqrt{\frac{k}{\mu}}

    Every harmonic-oscillator transition sits at exactly this same
    wavenumber (evenly spaced levels), unlike the Morse oscillator's
    overtones (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch.
    12.4).

    Parameters
    ----------
    force_constant : float
        Bond force constant `k`, in N/m.
    reduced_mass : float
        Vibrational reduced mass, in kg.

    Returns
    -------
    float
        Wavenumber, in cm^-1.

    Examples
    --------
    >>> import scipy.constants as sc
    >>> mu = (1.008 * 34.97) / (1.008 + 34.97) * sc.atomic_mass  # HCl-like reduced mass
    >>> nu = harmonic_fundamental_wavenumber(force_constant=480.0, reduced_mass=mu)
    >>> bool(2500.0 < nu < 3200.0)  # HCl's real fundamental is 2886 cm^-1
    True
    """
    ho = QuantumHarmonicOscillator(mass=reduced_mass, force_constant=force_constant)
    delta_e = ho.energy(1) - ho.energy(0)
    return float(energy_to_wavenumber(delta_e))


def morse_transition_wavenumbers(morse: MorseOscillator, v_max: int) -> np.ndarray:
    r"""Wavenumbers of the :math:`0\to v` Morse-oscillator "overtone" transitions, :math:`v=1,\dots,v_{max}`.

    Unlike the harmonic oscillator, successive Morse transitions
    (fundamental :math:`0\to1`, first overtone :math:`0\to2`, ...) are
    *not* simple multiples of each other -- the level spacing shrinks
    with `v` (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch.
    12.4).

    Parameters
    ----------
    morse : chemistrykit.quantum.systems.harmonic_oscillator.MorseOscillator
    v_max : int
        Highest final vibrational level to include.

    Returns
    -------
    ndarray, shape (v_max,)
        Wavenumbers of the :math:`0\to1,0\to2,\dots,0\to v_{max}`
        transitions, in cm^-1.

    Examples
    --------
    The first overtone (:math:`0\to2`) falls short of exactly twice the
    fundamental (:math:`0\to1`) -- the textbook anharmonicity signature:

    >>> morse = MorseOscillator(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19)
    >>> wavenumbers = morse_transition_wavenumbers(morse, v_max=2)
    >>> bool(wavenumbers[1] < 2.0 * wavenumbers[0])
    True
    """
    v_final = np.arange(1, v_max + 1)
    delta_e = morse.energy(v_final) - morse.energy(0)
    return energy_to_wavenumber(delta_e)


def anharmonicity_from_overtones(wavenumber_01: float, wavenumber_02: float) -> tuple:
    r"""Recover the Morse parameters :math:`\omega_e`, :math:`\omega_ex_e` from two measured band positions.

    The Morse-oscillator term-value formula gives the fundamental and
    first-overtone wavenumbers as

    .. math::

        \tilde\nu_{0\to1} = \omega_e - 2\omega_ex_e, \qquad
        \tilde\nu_{0\to2} = 2\omega_e - 6\omega_ex_e

    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 12.4;
    Herzberg, *Molecular Spectra and Molecular Structure I*, 2nd ed., Ch.
    III), a linear system solved here exactly (no fitting) for the two
    spectroscopic constants :math:`\omega_e` (the harmonic wavenumber)
    and :math:`\omega_ex_e` (the anharmonicity constant) from two
    observed transitions -- the standard way real vibrational
    anharmonicity constants are extracted from an observed overtone
    spectrum.

    Parameters
    ----------
    wavenumber_01 : float
        Observed :math:`0\to1` (fundamental) wavenumber, in cm^-1.
    wavenumber_02 : float
        Observed :math:`0\to2` (first overtone) wavenumber, in cm^-1.

    Returns
    -------
    omega_e : float
        Harmonic wavenumber, in cm^-1.
    omega_e_xe : float
        Anharmonicity constant, in cm^-1.

    Examples
    --------
    Round-trip against a :class:`~chemistrykit.quantum.systems.harmonic_oscillator.MorseOscillator`
    built from known `omega_e`/`omega_e_xe`-equivalent parameters
    recovers them from its own predicted band positions:

    >>> morse = MorseOscillator(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19)
    >>> nu01, nu02 = morse_transition_wavenumbers(morse, v_max=2)
    >>> omega_e, omega_e_xe = anharmonicity_from_overtones(nu01, nu02)
    >>> bool(omega_e > nu01 > 0)  # the harmonic wavenumber exceeds the (red-shifted) observed fundamental
    True
    >>> bool(omega_e_xe > 0)
    True
    """
    # [[1, -2], [2, -6]] @ [omega_e, omega_e_xe] = [nu01, nu02]
    A = np.array([[1.0, -2.0], [2.0, -6.0]])
    b = np.array([wavenumber_01, wavenumber_02])
    omega_e, omega_e_xe = np.linalg.solve(A, b)
    return float(omega_e), float(omega_e_xe)


@dataclass
class NormalModeResult:
    """The output of :meth:`TriatomicNormalModes.solve`."""

    wavenumbers: np.ndarray
    """ndarray, shape (3,): Normal-mode wavenumbers, in cm^-1, ascending."""

    internal_coordinate_vectors: np.ndarray
    """ndarray, shape (3, 3): Column `k` gives mode `k`'s displacement
    pattern in the (r1, r2, bend) internal-coordinate basis (see
    :class:`TriatomicNormalModes`), useful for classifying a mode as
    (anti)symmetric stretch or bend by inspection of the relative signs
    and magnitudes of its `r1`/`r2` components."""

    is_linear: bool
    """bool: Whether the equilibrium geometry is linear (in which case
    the computed bend mode is one member of a doubly degenerate pair --
    see :class:`TriatomicNormalModes`'s docstring)."""


class TriatomicNormalModes:
    r"""Normal-mode vibrational frequencies of an A-B-A triatomic, via the Wilson GF-matrix method.

    A genuine (not tabulated) normal-mode calculation: an internal-
    coordinate force field is combined with the molecule's actual
    equilibrium geometry and atomic masses through Wilson's `G` (kinetic
    energy) and `F` (force constant) matrices, and the resulting
    generalized eigenvalue problem is solved for the true normal-mode
    frequencies (Wilson, Decius & Cross, *Molecular Vibrations*, 1955,
    Ch. 4). Three internal coordinates are used: the two bond stretches
    :math:`r_1`, :math:`r_2`, and the bond angle (or, for a linear
    equilibrium geometry, a linearized bend coordinate -- see below).

    **The `G` matrix** (mass-weighted, purely geometric) is built here
    *numerically*: the Wilson B-matrix,
    :math:`B_{ki}=\partial S_k/\partial x_i` (internal coordinate `k`
    with respect to Cartesian displacement `i`), is obtained by central
    finite differences of the internal-coordinate functions, then
    :math:`G=BM^{-1}B^T` with `M` the diagonal atomic-mass matrix -- an
    exact numerical evaluation of the standard analytic B-matrix
    formulas, not an independent approximation.

    **The `F` matrix** is a simple (diagonal) valence force field:
    :math:`F=\mathrm{diag}(k_{r_1},k_{r_2},k_\theta)`, with *no*
    stretch-stretch or stretch-bend interaction constants. **Flagged as
    an approximation**: a full experimental force field for a real
    molecule (e.g. CO2, H2O) generally does include small off-diagonal
    interaction terms, fit to reproduce observed frequencies to high
    precision; the diagonal approximation used here still gives the
    correct qualitative mode pattern (symmetric stretch, antisymmetric
    stretch, bend) and the right order of magnitude for reasonable
    literature-typical `k_r`/`k_theta`, but not spectroscopic-accuracy
    absolute wavenumbers.

    **Linear equilibrium geometries** (`is_linear=True`) need a bend
    coordinate that stays well-defined exactly at 180 degrees, where the
    ordinary :math:`\theta=\arccos(\dots)` bond-angle formula's
    derivative diverges. The standard fix (Wilson, Decius & Cross, Ch.
    4) is used: the bend coordinate is instead the linearized transverse
    (in-plane) displacement combination
    :math:`S_{\text{bend}}=(x_0-2x_1+x_2)/r_{eq}`, well-defined at
    exactly linear. A linear triatomic's bend is genuinely doubly
    degenerate (in-plane and out-of-plane bending are equivalent by the
    molecule's cylindrical symmetry); this module solves only the
    in-plane 3-internal-coordinate problem, so the reported bend
    wavenumber should be understood as representing *both* degenerate
    components (Herzberg, *Molecular Spectra and Molecular Structure
    II*, 1945, Ch. I.3).

    The resulting symmetric generalized eigenvalue problem, `GFL=L\\Lambda`,
    is solved here as the equivalent symmetric form
    :math:`Fc=\lambda G^{-1}c` via
    :func:`chemistrykit.quantum.utils.secular_equation.solve_secular_equation`
    (with `F` playing the role of the Hamiltonian and :math:`G^{-1}` the
    role of the overlap matrix) -- reusing the same shared-eigenvalue-
    solve machinery :mod:`chemistrykit.quantum`'s variational models use,
    since both are literally the same kind of generalized symmetric
    eigenproblem, :math:`\lambda=\omega^2` in SI units, converted to a
    wavenumber via :math:`\tilde\nu=\omega/(2\pi c)`.

    Parameters
    ----------
    masses : array-like of float, shape (3,)
        Atomic masses (atom0, central atom1, atom2), in kg.
    equilibrium_coordinates : array-like of float, shape (3, 3)
        Cartesian equilibrium coordinates, in meters.
    k_r1, k_r2 : float
        Stretch force constants for bonds (atom0-atom1) and
        (atom2-atom1), in N/m.
    k_theta : float
        Bend force constant, in J/rad^2 (equivalently N m/rad^2).

    Examples
    --------
    See :meth:`linear` and :meth:`bent` for the usual way to construct
    one; :meth:`solve` returns three positive wavenumbers for a stable
    equilibrium geometry.
    """

    def __init__(self, masses, equilibrium_coordinates, k_r1: float, k_r2: float, k_theta: float):
        self.masses = np.asarray(masses, dtype=np.float64)
        if self.masses.shape != (3,):
            raise ValueError("masses must have shape (3,)")
        self.equilibrium_coordinates = np.asarray(equilibrium_coordinates, dtype=np.float64)
        if self.equilibrium_coordinates.shape != (3, 3):
            raise ValueError("equilibrium_coordinates must have shape (3, 3)")
        if k_r1 <= 0 or k_r2 <= 0 or k_theta <= 0:
            raise ValueError("force constants must be positive")
        self.k_r1 = float(k_r1)
        self.k_r2 = float(k_r2)
        self.k_theta = float(k_theta)
        self.is_linear = self._check_linear()
        self._r_eq = np.linalg.norm(self.equilibrium_coordinates[0] - self.equilibrium_coordinates[1])

    def _check_linear(self, tol: float = 1e-6) -> bool:
        v1 = self.equilibrium_coordinates[0] - self.equilibrium_coordinates[1]
        v2 = self.equilibrium_coordinates[2] - self.equilibrium_coordinates[1]
        cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return bool(cos_theta < -1.0 + tol)

    @classmethod
    def linear(cls, mass_terminal: float, mass_central: float, bond_length: float, k_r: float, k_theta: float) -> TriatomicNormalModes:
        """Build a symmetric linear A-B-A triatomic (e.g. CO2-like) along the `z` axis.

        Parameters
        ----------
        mass_terminal : float
            Mass of each terminal atom (A), in kg.
        mass_central : float
            Mass of the central atom (B), in kg.
        bond_length : float
            Equilibrium A-B bond length, in meters.
        k_r : float
            Shared stretch force constant for both bonds, in N/m.
        k_theta : float
            Bend force constant, in J/rad^2.

        Returns
        -------
        TriatomicNormalModes
        """
        coordinates = np.array([[0.0, 0.0, -bond_length], [0.0, 0.0, 0.0], [0.0, 0.0, bond_length]])
        masses = np.array([mass_terminal, mass_central, mass_terminal])
        return cls(masses, coordinates, k_r, k_r, k_theta)

    @classmethod
    def bent(cls, mass_terminal: float, mass_central: float, bond_length: float, bond_angle_degrees: float, k_r: float, k_theta: float) -> TriatomicNormalModes:
        """Build a symmetric bent A-B-A triatomic (e.g. H2O-like) in the `xz` plane.

        Parameters
        ----------
        mass_terminal : float
            Mass of each terminal atom (A), in kg.
        mass_central : float
            Mass of the central atom (B), in kg.
        bond_length : float
            Equilibrium A-B bond length, in meters.
        bond_angle_degrees : float
            Equilibrium A-B-A bond angle, in degrees (< 180).
        k_r : float
            Shared stretch force constant for both bonds, in N/m.
        k_theta : float
            Bend force constant, in J/rad^2.

        Returns
        -------
        TriatomicNormalModes
        """
        half_angle = np.radians(bond_angle_degrees) / 2.0
        coordinates = np.array(
            [
                [bond_length * np.sin(half_angle), 0.0, bond_length * np.cos(half_angle)],
                [0.0, 0.0, 0.0],
                [-bond_length * np.sin(half_angle), 0.0, bond_length * np.cos(half_angle)],
            ]
        )
        masses = np.array([mass_terminal, mass_central, mass_terminal])
        return cls(masses, coordinates, k_r, k_r, k_theta)

    def _internal_coordinates(self, cartesian: np.ndarray) -> np.ndarray:
        r1 = np.linalg.norm(cartesian[0] - cartesian[1])
        r2 = np.linalg.norm(cartesian[2] - cartesian[1])
        if self.is_linear:
            bend = (cartesian[0, 0] - 2.0 * cartesian[1, 0] + cartesian[2, 0]) / self._r_eq
        else:
            v1 = (cartesian[0] - cartesian[1]) / r1
            v2 = (cartesian[2] - cartesian[1]) / r2
            cos_theta = np.clip(np.dot(v1, v2), -1.0, 1.0)
            bend = np.arccos(cos_theta)
        return np.array([r1, r2, bend])

    def _b_matrix(self) -> np.ndarray:
        flat = self.equilibrium_coordinates.flatten()
        n = flat.size
        h = 1e-6 * self._r_eq
        B = np.zeros((3, n))
        for j in range(n):
            step = np.zeros(n)
            step[j] = h
            plus = self._internal_coordinates((flat + step).reshape(3, 3))
            minus = self._internal_coordinates((flat - step).reshape(3, 3))
            B[:, j] = (plus - minus) / (2.0 * h)
        return B

    def _g_matrix(self) -> np.ndarray:
        B = self._b_matrix()
        inverse_masses = np.repeat(1.0 / self.masses, 3)
        return B @ np.diag(inverse_masses) @ B.T

    def _f_matrix(self) -> np.ndarray:
        return np.diag([self.k_r1, self.k_r2, self.k_theta])

    def solve(self) -> NormalModeResult:
        r"""Diagonalize the GF-matrix eigenproblem and return the normal-mode wavenumbers.

        Returns
        -------
        NormalModeResult

        Examples
        --------
        A bent triatomic (H2O-like) has 3 real, positive vibrational
        wavenumbers, with the bend lower in energy than either stretch --
        the correct qualitative ordering for a real bent AX2 molecule:

        >>> import scipy.constants as sc
        >>> water = TriatomicNormalModes.bent(
        ...     mass_terminal=1.008 * sc.atomic_mass,
        ...     mass_central=15.999 * sc.atomic_mass,
        ...     bond_length=95.8e-12,
        ...     bond_angle_degrees=104.5,
        ...     k_r=770.0,
        ...     k_theta=0.7e-18,
        ... )
        >>> result = water.solve()
        >>> bool(np.all(result.wavenumbers > 0))
        True
        >>> bool(result.wavenumbers[0] < result.wavenumbers[1] < result.wavenumbers[2])
        True
        """
        F = self._f_matrix()
        G = self._g_matrix()
        G_inv = np.linalg.inv(G)
        try:
            eigenvalues, eigenvectors = solve_secular_equation(F, G_inv)
        except ValueError:
            # G_inv from a numerically imperfect linear geometry can pick
            # up a tiny asymmetry that solve_secular_equation's strict
            # symmetry check rejects; scipy's eigh tolerates it directly.
            eigenvalues, eigenvectors = sla.eigh(F, G_inv)
        order = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]
        angular_frequencies = np.sqrt(np.clip(eigenvalues, 0.0, None))
        wavenumbers = angular_frequencies / (2.0 * np.pi * C) / 100.0
        return NormalModeResult(wavenumbers=wavenumbers, internal_coordinate_vectors=eigenvectors, is_linear=self.is_linear)
