r"""First-order NMR multiplet simulation: chemical shifts and J-coupling splitting patterns.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 13.7-13.8, or
Silverstein, Webster & Kiemle, *Spectrometric Identification of Organic
Compounds*, 8th ed., Ch. 4, throughout.

The first-order (weak-coupling) approximation applies when a coupling
constant `J` is much smaller than the chemical-shift difference between
the coupled nuclei (in Hz, i.e. after multiplying a ppm difference by the
spectrometer frequency); under that approximation, `n` magnetically
equivalent spin-1/2 neighbors split a resonance into `n+1` lines with
relative intensities given by the binomial coefficients (Pascal's
triangle), symmetric about the unperturbed chemical shift (N. F. Ramsey
& E. M. Purcell, *Phys. Rev.* 85, 143 (1952), for the coupling
mechanism; Pascal's triangle itself long predates NMR). Coupling to
several *inequivalent* sets of neighbors (different `J` values) splits
the line independently and multiplicatively -- each new splitting acts
on every line already present.
"""

from __future__ import annotations

import numpy as np
from scipy.special import comb

from chemistrykit.spectro.core.base_system import Spectrum

__all__ = ["multiplicity", "pascals_triangle_intensities", "first_order_multiplet", "multi_coupling_multiplet"]


def multiplicity(n_equivalent_neighbors: int) -> int:
    r"""The n+1 rule: `n` equivalent spin-1/2 neighbors split a resonance into `n+1` lines.

    Parameters
    ----------
    n_equivalent_neighbors : int
        Number of magnetically equivalent, spin-1/2, first-order-coupled
        neighboring nuclei.

    Returns
    -------
    int

    Examples
    --------
    Ethanol's -CH2- protons, coupled to the 3 equivalent -CH3- protons,
    form a quartet:

    >>> multiplicity(3)
    4
    """
    if n_equivalent_neighbors < 0:
        raise ValueError("n_equivalent_neighbors must be >= 0")
    return n_equivalent_neighbors + 1


def pascals_triangle_intensities(n_equivalent_neighbors: int) -> np.ndarray:
    r"""Relative line intensities for `n` equivalent spin-1/2 neighbors: the binomial coefficients :math:`\binom{n}{k}`.

    Parameters
    ----------
    n_equivalent_neighbors : int

    Returns
    -------
    ndarray, shape (n_equivalent_neighbors + 1,)

    Examples
    --------
    A quartet's classic 1:3:3:1 intensity pattern (e.g. ethanol's CH2,
    coupled to 3 equivalent CH3 protons):

    >>> pascals_triangle_intensities(3)
    array([1., 3., 3., 1.])
    """
    n = int(n_equivalent_neighbors)
    return np.array([comb(n, k) for k in range(n + 1)], dtype=np.float64)


def first_order_multiplet(chemical_shift_ppm: float, j_coupling_hz: float, n_neighbors: int, spectrometer_frequency_mhz: float) -> Spectrum:
    r"""Simulate a simple multiplet from coupling to one set of `n` equivalent neighbors.

    Line positions are evenly spaced by :math:`J/\nu_0` in ppm (`J` in Hz
    converted to ppm by dividing by the spectrometer frequency in MHz,
    since 1 ppm on the chemical-shift scale corresponds to
    :math:`\nu_0\,[\text{Hz}]/10^6` at operating frequency
    :math:`\nu_0`), symmetric about `chemical_shift_ppm`, with the
    Pascal's-triangle relative intensities from
    :func:`pascals_triangle_intensities`.

    Parameters
    ----------
    chemical_shift_ppm : float
        Unperturbed chemical shift, in ppm.
    j_coupling_hz : float
        Coupling constant `J`, in Hz.
    n_neighbors : int
        Number of equivalent coupled neighbors.
    spectrometer_frequency_mhz : float
        Spectrometer (proton, or whichever nucleus) operating frequency,
        in MHz -- needed to convert `J` from Hz to the ppm scale.

    Returns
    -------
    Spectrum
        `positions` in ppm, `intensities` the relative (unnormalized)
        binomial line intensities.

    Examples
    --------
    A triplet, symmetric about the chemical shift, spaced by `J`:

    >>> spectrum = first_order_multiplet(chemical_shift_ppm=1.2, j_coupling_hz=7.0, n_neighbors=2, spectrometer_frequency_mhz=400.0)
    >>> len(spectrum.positions)
    3
    >>> round(float(np.mean(spectrum.positions)), 6) == round(1.2, 6)
    True
    >>> spectrum.intensities.tolist()
    [1.0, 2.0, 1.0]
    """
    j_ppm = j_coupling_hz / spectrometer_frequency_mhz
    n = multiplicity(n_neighbors) - 1
    offsets = j_ppm * (np.arange(0, n + 1) - n / 2.0)
    positions = chemical_shift_ppm + offsets
    intensities = pascals_triangle_intensities(n)
    return Spectrum(positions=positions, intensities=intensities)


def multi_coupling_multiplet(chemical_shift_ppm: float, couplings, spectrometer_frequency_mhz: float) -> Spectrum:
    r"""Simulate a multiplet from independent first-order coupling to several *inequivalent* neighbor sets.

    Each coupling constant splits every line already present (the
    couplings act independently and multiplicatively), giving
    :math:`\prod_k(n_k+1)` lines in the general case (fewer if some
    coincide) -- e.g. a doublet of triplets (dt) from one neighbor with
    :math:`J_1` and two equivalent neighbors with :math:`J_2\ne J_1`.

    Parameters
    ----------
    chemical_shift_ppm : float
        Unperturbed chemical shift, in ppm.
    couplings : sequence of tuple(float, int)
        `(J_hz, n_equivalent_neighbors)` for each independent coupling
        partner set.
    spectrometer_frequency_mhz : float

    Returns
    -------
    Spectrum
        `positions` in ppm (sorted ascending), `intensities` the combined
        relative intensities (coincident positions within numerical
        tolerance are merged, summing their intensities).

    Examples
    --------
    A doublet of doublets (dd) with equal J's degenerates into a triplet
    with 1:2:1 intensities (two different coupling partners with the
    same `J` are indistinguishable from one set of two equivalent
    partners):

    >>> spectrum = multi_coupling_multiplet(chemical_shift_ppm=5.0, couplings=[(7.0, 1), (7.0, 1)], spectrometer_frequency_mhz=400.0)
    >>> len(spectrum.positions)
    3
    >>> spectrum.intensities.tolist()
    [1.0, 2.0, 1.0]

    A genuine doublet of triplets (:math:`J_1=12` Hz to one neighbor,
    :math:`J_2=5` Hz to two equivalent neighbors) has 6 distinct lines:

    >>> dt = multi_coupling_multiplet(chemical_shift_ppm=3.0, couplings=[(12.0, 1), (5.0, 2)], spectrometer_frequency_mhz=400.0)
    >>> len(dt.positions)
    6
    >>> round(float(np.sum(dt.intensities)), 6)
    8.0
    """
    positions = np.array([chemical_shift_ppm])
    intensities = np.array([1.0])
    for j_hz, n in couplings:
        j_ppm = j_hz / spectrometer_frequency_mhz
        offsets = j_ppm * (np.arange(0, n + 1) - n / 2.0)
        sub_intensities = pascals_triangle_intensities(n)
        new_positions = (positions[:, None] + offsets[None, :]).ravel()
        new_intensities = (intensities[:, None] * sub_intensities[None, :]).ravel()
        positions, intensities = new_positions, new_intensities

    order = np.argsort(positions)
    positions, intensities = positions[order], intensities[order]
    merged_positions: list = []
    merged_intensities: list = []
    tol = 1e-9
    for pos, inten in zip(positions, intensities, strict=True):
        if merged_positions and abs(pos - merged_positions[-1]) < tol:
            merged_intensities[-1] += inten
        else:
            merged_positions.append(pos)
            merged_intensities.append(inten)
    return Spectrum(positions=np.array(merged_positions), intensities=np.array(merged_intensities))
