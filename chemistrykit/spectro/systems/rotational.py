r"""Rigid-rotor rotational (microwave) spectra: transition wavenumbers, relative intensities, and isotope shifts.

Builds on :class:`chemistrykit.quantum.systems.rigid_rotor.RigidRotor`
for the underlying energy levels and :math:`\Delta J=\pm1` selection-rule
transition energies -- this module does not re-derive the rigid-rotor
eigenvalues, only converts them to the spectroscopic quantities
(wavenumber positions, a relative-intensity model, isotope shifts) an
actual rotational spectrum is reported in. See Atkins & de Paula,
*Physical Chemistry*, 11th ed., Ch. 12.2 (rotational spectroscopy), and
Ch. 12.3 (isotope effects) throughout.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import K_B, C, H
from chemistrykit.quantum.systems.rigid_rotor import RigidRotor
from chemistrykit.spectro.core.base_system import Spectrum

__all__ = ["energy_to_wavenumber", "rotational_line_wavenumbers", "rotational_spectrum", "isotope_shift_ratio"]


def energy_to_wavenumber(energy_j) -> float:
    r"""Convert an energy (in J) to a spectroscopic wavenumber :math:`\tilde\nu=E/(hc)`, in cm^-1.

    Parameters
    ----------
    energy_j : float or array-like of float

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> round(float(energy_to_wavenumber(H * C * 100.0)), 6)
    1.0
    """
    energy_j = np.asarray(energy_j, dtype=np.float64)
    return energy_j / (H * C) / 100.0


def rotational_line_wavenumbers(rotor: RigidRotor, J_max: int) -> np.ndarray:
    r"""Wavenumbers of the :math:`J\to J+1` rotational absorption lines, :math:`J=0,\dots,J_{max}`.

    Each line sits at :math:`\tilde\nu_J=2B(J+1)` (in wavenumber units),
    evenly spaced by :math:`2B` -- the defining signature of a rigid-
    rotor rotational spectrum (Atkins & de Paula, *Physical Chemistry*,
    11th ed., Ch. 12.2).

    Parameters
    ----------
    rotor : chemistrykit.quantum.systems.rigid_rotor.RigidRotor
    J_max : int
        Highest initial `J` to include.

    Returns
    -------
    ndarray, shape (J_max + 1,)
        Wavenumbers, in cm^-1.

    Examples
    --------
    >>> import scipy.constants as sc
    >>> rotor = RigidRotor.from_diatomic(mass1=1.008 * sc.atomic_mass, mass2=34.97 * sc.atomic_mass, bond_length=127.5e-12)
    >>> lines = rotational_line_wavenumbers(rotor, J_max=4)
    >>> spacings = np.diff(lines)
    >>> bool(np.allclose(spacings, spacings[0]))
    True
    """
    J_values = np.arange(0, J_max + 1)
    transition_energies = np.array([rotor.transition_energy(int(J)) for J in J_values])
    return energy_to_wavenumber(transition_energies)


def rotational_spectrum(rotor: RigidRotor, J_max: int, temperature: float) -> Spectrum:
    r"""Build a rotational-spectrum :class:`~chemistrykit.spectro.core.base_system.Spectrum`, with relative intensities.

    **Approximation flagged explicitly**: the relative intensity of the
    :math:`J\to J+1` line is modeled here as the thermal (Boltzmann)
    population of the initial level `J`, :math:`(2J+1)e^{-E_J/k_BT}`,
    times the degeneracy-averaged squared transition dipole
    :math:`|\mu_{J+1,J}|^2\propto(J+1)/(2J+1)` of a linear rotor (Atkins &
    de Paula, *Physical Chemistry*, 11th ed., Ch. 12.2(b); Bernath,
    *Spectra of Atoms and Molecules*, 2nd ed., Ch. 6.5). The
    :math:`(2J+1)` factors cancel, leaving

    .. math::

        I_J \propto (J+1)\exp\!\left[-\frac{E_J}{k_BT}\right]

    The frequency factor and stimulated-emission correction of a full
    absorption-coefficient calculation are omitted.

    This reproduces the qualitative textbook feature of a rotational
    spectrum -- intensity rising from `J=0`, peaking at some intermediate
    `J`, then falling off -- but is not a first-principles transition-
    dipole calculation for any particular real molecule.

    Parameters
    ----------
    rotor : chemistrykit.quantum.systems.rigid_rotor.RigidRotor
    J_max : int
        Highest initial `J` to include.
    temperature : float
        Temperature, in K, for the Boltzmann population weighting.

    Returns
    -------
    Spectrum
        `positions` in cm^-1; `intensities` relative (not normalized to
        any absolute scale).

    Examples
    --------
    At a low enough temperature, only the lowest few `J` are
    significantly populated, so the peak intensity line has a small `J`:

    >>> import scipy.constants as sc
    >>> rotor = RigidRotor.from_diatomic(mass1=1.008 * sc.atomic_mass, mass2=34.97 * sc.atomic_mass, bond_length=127.5e-12)
    >>> spectrum = rotational_spectrum(rotor, J_max=15, temperature=300.0)
    >>> peak_J = int(np.argmax(spectrum.intensities))
    >>> bool(0 < peak_J < 15)
    True
    """
    J_values = np.arange(0, J_max + 1)
    positions = rotational_line_wavenumbers(rotor, J_max)
    energies = np.array([rotor.energy(int(J)) for J in J_values])
    intensities = (J_values + 1.0) * np.exp(-energies / (K_B * temperature))
    labels = [f"J={J}->{J + 1}" for J in J_values]
    return Spectrum(positions=positions, intensities=intensities, labels=labels)


def isotope_shift_ratio(mass1a: float, mass2: float, mass1b: float) -> float:
    r"""Ratio of rotational constants (and transition wavenumbers) between two isotopologues.

    Substituting atom 1 (e.g. H) for a heavier isotope (e.g. D) changes
    only the reduced mass :math:`\mu=m_1m_2/(m_1+m_2)`, and the
    rotational constant :math:`B\propto1/I=1/(\mu r^2)` scales inversely
    with it (bond length `r` is, to a very good approximation, unchanged
    by isotopic substitution, since it depends on the electronic
    structure, not the nuclear mass -- the Born-Oppenheimer
    approximation). So

    .. math::

        \frac{B(\text{isotopologue }b)}{B(\text{isotopologue }a)}
        = \frac{\mu_a}{\mu_b}

    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 12.3). Every
    rotational transition wavenumber scales by this same ratio, since
    they are all proportional to `B`.

    Parameters
    ----------
    mass1a : float
        Mass of the lighter isotope of the substituted atom (e.g. H), in
        any consistent unit.
    mass2 : float
        Mass of the *other* (unsubstituted) atom, same units.
    mass1b : float
        Mass of the heavier isotope of the substituted atom (e.g. D),
        same units.

    Returns
    -------
    float
        :math:`B_b/B_a` (and :math:`\tilde\nu_b/\tilde\nu_a`), `< 1`
        since substituting a heavier isotope increases the reduced mass
        and moment of inertia, lowering `B`.

    Examples
    --------
    Substituting D for H in HCl roughly halves the rotational constant
    (deuterium is about twice as heavy as hydrogen, and Cl is much
    heavier than either, so the reduced mass roughly doubles):

    >>> import scipy.constants as sc
    >>> ratio = isotope_shift_ratio(mass1a=1.008 * sc.atomic_mass, mass2=34.97 * sc.atomic_mass, mass1b=2.014 * sc.atomic_mass)
    >>> bool(0.45 < ratio < 0.55)
    True

    Directly matches a from-scratch :class:`~chemistrykit.quantum.systems.rigid_rotor.RigidRotor`
    comparison:

    >>> r = 127.5e-12
    >>> rotor_hcl = RigidRotor.from_diatomic(1.008 * sc.atomic_mass, 34.97 * sc.atomic_mass, r)
    >>> rotor_dcl = RigidRotor.from_diatomic(2.014 * sc.atomic_mass, 34.97 * sc.atomic_mass, r)
    >>> ratio_direct = rotor_dcl.rotational_constant / rotor_hcl.rotational_constant
    >>> bool(abs(ratio_direct - ratio) < 1e-9)
    True
    """
    mu_a = mass1a * mass2 / (mass1a + mass2)
    mu_b = mass1b * mass2 / (mass1b + mass2)
    return float(mu_a / mu_b)
