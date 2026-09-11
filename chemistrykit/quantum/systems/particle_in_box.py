r"""The particle in a box (1D and 3D), and its free-electron-model application to conjugated-dye UV-Vis absorption.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 7.3 for the
particle-in-a-box energy levels and wavefunctions, and Ch. 7's "Impact on
nanoscience" (or Levine, *Quantum Chemistry*, 7th ed., Ch. 2.2c) for the
Kuhn free-electron model of a linear conjugated polyene/cyanine dye as a
1D box of the pi electrons (H. Kuhn, *J. Chem. Phys.* 17, 1198 (1949)).
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import ELECTRON_MASS, C, H
from chemistrykit.quantum.core.base_system import QuantumSystem

__all__ = ["ParticleInBox1D", "ParticleInBox3D", "conjugated_dye_absorption_wavelength"]


class ParticleInBox1D(QuantumSystem):
    r"""A particle of mass `m` confined to an infinite 1D box of length `L`, :math:`0\le x\le L`.

    .. math::

        E_n = \frac{n^2h^2}{8mL^2}, \qquad
        \psi_n(x) = \sqrt{\frac2L}\sin\!\left(\frac{n\pi x}{L}\right),
        \qquad n=1,2,3,\dots

    (Atkins & de Paula, *Physical Chemistry*, 11th ed., eq. 7.6a-b.)

    Parameters
    ----------
    length : float
        Box length `L`, in m.
    mass : float, default :data:`chemistrykit.constants.ELECTRON_MASS`
        Particle mass, in kg.

    Examples
    --------
    Doubling the box length cuts every energy level by a factor of 4
    (:math:`E\propto1/L^2`):

    >>> box1 = ParticleInBox1D(length=1.0e-9)
    >>> box2 = ParticleInBox1D(length=2.0e-9)
    >>> round(float(box1.energy(1) / box2.energy(1)), 6)
    4.0
    """

    def __init__(self, length: float, mass: float = ELECTRON_MASS):
        if length <= 0:
            raise ValueError("length must be positive")
        if mass <= 0:
            raise ValueError("mass must be positive")
        self.length = float(length)
        self.mass = float(mass)

    def energy(self, n):
        """Return :math:`E_n=n^2h^2/(8mL^2)`.

        Parameters
        ----------
        n : int or array-like of int
            Quantum number(s), each >= 1.

        Returns
        -------
        float or ndarray
            Energy, in J.

        Examples
        --------
        >>> box = ParticleInBox1D(length=1.0e-9)  # a 1 nm box
        >>> round(float(box.energy(1) / 1.602176634e-19), 4)  # in eV
        0.376
        """
        n = np.asarray(n, dtype=np.float64)
        if np.any(n < 1):
            raise ValueError("n must be >= 1")
        return n**2 * H**2 / (8.0 * self.mass * self.length**2)

    def wavefunction(self, x, n: int):
        r"""Return :math:`\psi_n(x)=\sqrt{2/L}\sin(n\pi x/L)`.

        Parameters
        ----------
        x : float or array-like of float
            Position(s), in m; should lie in ``[0, length]`` (zero
            outside, per the infinite-wall boundary condition, but this
            is not enforced -- callers restrict the domain themselves).
        n : int
            Quantum number, >= 1.

        Returns
        -------
        float or ndarray
            Amplitude, in m^-1/2.

        Examples
        --------
        The wavefunction vanishes at both walls (the boundary condition
        that quantizes `n` in the first place):

        >>> box = ParticleInBox1D(length=1.0e-9)
        >>> bool(abs(box.wavefunction(0.0, n=1)) < 1e-9)
        True
        >>> bool(abs(box.wavefunction(1.0e-9, n=1)) < 1e-9)
        True
        """
        if n < 1:
            raise ValueError("n must be >= 1")
        x = np.asarray(x, dtype=np.float64)
        return np.sqrt(2.0 / self.length) * np.sin(n * np.pi * x / self.length)

    def probability_density(self, x, n: int):
        r"""Return :math:`|\psi_n(x)|^2`.

        Parameters
        ----------
        x : float or array-like of float
        n : int

        Returns
        -------
        float or ndarray
            Probability density, in m^-1.
        """
        return self.wavefunction(x, n) ** 2


class ParticleInBox3D(QuantumSystem):
    r"""A particle confined to an infinite rectangular box of dimensions `(Lx, Ly, Lz)`.

    .. math::

        E_{n_x,n_y,n_z} = \frac{h^2}{8m}\left(\frac{n_x^2}{L_x^2}+\frac{n_y^2}{L_y^2}+\frac{n_z^2}{L_z^2}\right)

    (Atkins & de Paula, *Physical Chemistry*, 11th ed., eq. 7.10.) A cubic
    box (:math:`L_x=L_y=L_z`) shows genuine degeneracy -- e.g. the
    :math:`(2,1,1)`, :math:`(1,2,1)`, :math:`(1,1,2)` states are
    degenerate -- a textbook illustration of how symmetry produces
    degenerate levels (Atkins & de Paula, *Physical Chemistry*, 11th ed.,
    Ch. 7.3(c)).

    Parameters
    ----------
    Lx, Ly, Lz : float
        Box dimensions, in m.
    mass : float, default :data:`chemistrykit.constants.ELECTRON_MASS`
        Particle mass, in kg.

    Examples
    --------
    >>> box = ParticleInBox3D(Lx=1.0e-9, Ly=1.0e-9, Lz=1.0e-9)
    >>> e211 = box.energy(2, 1, 1)
    >>> e121 = box.energy(1, 2, 1)
    >>> e112 = box.energy(1, 1, 2)
    >>> round(e211, 30) == round(e121, 30) == round(e112, 30)
    True
    """

    def __init__(self, Lx: float, Ly: float, Lz: float, mass: float = ELECTRON_MASS):
        for L in (Lx, Ly, Lz):
            if L <= 0:
                raise ValueError("box dimensions must be positive")
        if mass <= 0:
            raise ValueError("mass must be positive")
        self.Lx, self.Ly, self.Lz = float(Lx), float(Ly), float(Lz)
        self.mass = float(mass)

    def energy(self, nx, ny, nz):
        """Return :math:`E_{n_x,n_y,n_z}`.

        Parameters
        ----------
        nx, ny, nz : int
            Quantum numbers, each >= 1.

        Returns
        -------
        float
            Energy, in J.
        """
        if nx < 1 or ny < 1 or nz < 1:
            raise ValueError("quantum numbers must be >= 1")
        return (H**2 / (8.0 * self.mass)) * (nx**2 / self.Lx**2 + ny**2 / self.Ly**2 + nz**2 / self.Lz**2)

    def degeneracy(self, n_max: int) -> dict:
        """Enumerate energies (rounded) and their degeneracies up to ``nx, ny, nz <= n_max``.

        Parameters
        ----------
        n_max : int
            Largest quantum number to include along each axis.

        Returns
        -------
        dict
            Maps a (relative-tolerance-rounded) energy value, in J, to the
            number of ``(nx, ny, nz)`` triples sharing it.

        Examples
        --------
        A cubic box has a 3-fold-degenerate first excited level:

        >>> box = ParticleInBox3D(Lx=1.0e-9, Ly=1.0e-9, Lz=1.0e-9)
        >>> degeneracies = box.degeneracy(n_max=2)
        >>> sorted(degeneracies.values())[-2]
        3
        """
        energies = {}
        for nx in range(1, n_max + 1):
            for ny in range(1, n_max + 1):
                for nz in range(1, n_max + 1):
                    e = round(self.energy(nx, ny, nz), 30)
                    key = next((k for k in energies if abs(k - e) < 1e-6 * abs(e) + 1e-30), e)
                    energies[key] = energies.get(key, 0) + 1
        return energies


def conjugated_dye_absorption_wavelength(box_length: float, n_pi_electrons: int, mass: float = ELECTRON_MASS) -> float:
    r"""Predict a linear conjugated dye's UV-Vis absorption wavelength via Kuhn's free-electron model.

    Models the delocalized pi electrons of a linear conjugated chain
    (e.g. a cyanine dye) as free particles in a 1D box spanning the
    conjugated system (H. Kuhn, *J. Chem. Phys.* 17, 1198 (1949); Levine,
    *Quantum Chemistry*, 7th ed., Ch. 2.2c). With `n_pi_electrons`
    electrons filling the box levels two at a time (Pauli exclusion), the
    HOMO is level :math:`n_{HOMO}=n_{\pi}/2` (only defined here for even
    `n_pi_electrons`, the closed-shell case) and the LUMO is
    :math:`n_{LUMO}=n_{HOMO}+1`. The HOMO->LUMO transition energy is then

    .. math::

        \Delta E = E_{n_{HOMO}+1}-E_{n_{HOMO}}
                 = \frac{h^2}{8mL^2}\left(2n_{HOMO}+1\right)
                 = \frac{h^2}{8mL^2}(n_\pi+1)

    and the predicted absorption wavelength follows from
    :math:`\Delta E=hc/\lambda`. This is a crude model (it ignores the
    actual sigma-bond framework's effect on the box walls and treats the
    pi electrons as fully independent/non-interacting), but it correctly
    predicts the qualitative trend that longer conjugated chains (more
    pi electrons, longer effective box) absorb at longer wavelength -- the
    basis of dye chemistry's color tuning by conjugation length.

    Parameters
    ----------
    box_length : float
        Effective 1D box length spanned by the conjugated pi system, in m
        (conventionally taken as the number of conjugated bonds times a
        typical C-C bond length, plus one bond length of "overhang" at
        each end).
    n_pi_electrons : int
        Number of pi electrons, must be a positive even integer
        (closed-shell filling).
    mass : float, default :data:`chemistrykit.constants.ELECTRON_MASS`

    Returns
    -------
    float
        Predicted absorption wavelength, in m.

    Raises
    ------
    ValueError
        If `n_pi_electrons` is not a positive even integer.

    Examples
    --------
    A cyanine-dye-like chain with 8 pi electrons over a ~1.2 nm box
    absorbs in the visible, as observed for real cyanine dyes:

    >>> wavelength = conjugated_dye_absorption_wavelength(box_length=1.2e-9, n_pi_electrons=8)
    >>> bool(3.8e-7 < wavelength < 8.0e-7)  # visible range
    True

    Adding more conjugation (more pi electrons over a longer box) red-shifts the absorption:

    >>> short = conjugated_dye_absorption_wavelength(box_length=0.8e-9, n_pi_electrons=6)
    >>> long = conjugated_dye_absorption_wavelength(box_length=1.4e-9, n_pi_electrons=10)
    >>> bool(long > short)
    True
    """
    if n_pi_electrons <= 0 or n_pi_electrons % 2 != 0:
        raise ValueError("n_pi_electrons must be a positive even integer (closed-shell filling)")
    delta_e = (H**2 / (8.0 * mass * box_length**2)) * (n_pi_electrons + 1)
    return H * C / delta_e
