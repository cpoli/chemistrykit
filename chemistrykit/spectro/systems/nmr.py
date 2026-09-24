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

__all__ = [
    "multiplicity",
    "pascals_triangle_intensities",
    "first_order_multiplet",
    "multi_coupling_multiplet",
    "larmor_frequency",
    "chemical_shift_ppm",
    "karplus_coupling",
    "free_induction_decay",
    "fid_to_spectrum",
]


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


def larmor_frequency(gyromagnetic_ratio: float, field_tesla):
    r"""Nuclear Larmor (resonance) frequency :math:`\nu_0=\gamma B_0/(2\pi)`, in Hz.

    The frequency at which Bloch's and Purcell's 1946 experiments found a
    sharp resonant response from nuclear spins in a static field
    :math:`B_0` (F. Bloch, *Phys. Rev.* 70, 460 (1946); Atkins & de
    Paula, *Physical Chemistry*, 11th ed., Ch. 12A).

    Parameters
    ----------
    gyromagnetic_ratio : float
        Nuclear gyromagnetic ratio :math:`\gamma`, in rad s^-1 T^-1
        (e.g. ``scipy.constants.physical_constants["proton gyromag. ratio"][0]``).
    field_tesla : float or array-like of float
        Static field :math:`B_0`, in tesla.

    Returns
    -------
    float or ndarray
        Resonance frequency, in Hz.

    Examples
    --------
    A "400 MHz" spectrometer's magnet is about 9.4 T:

    >>> import scipy.constants as sc
    >>> gamma_h = sc.physical_constants["proton gyromag. ratio"][0]
    >>> round(larmor_frequency(gamma_h, 9.3947) / 1e6, 1)
    400.0
    """
    result = gyromagnetic_ratio * np.asarray(field_tesla, dtype=np.float64) / (2.0 * np.pi)
    return float(result) if result.ndim == 0 else result


def chemical_shift_ppm(frequency_hz, reference_frequency_hz: float):
    r"""Chemical shift :math:`\delta=10^6(\nu-\nu_\text{ref})/\nu_\text{ref}`, in ppm.

    Dividing by the reference frequency removes the (field-proportional)
    Larmor frequency, so :math:`\delta` is the same on every
    spectrometer -- the reason NMR spectra are reported in ppm.

    Parameters
    ----------
    frequency_hz : float or array-like of float
        Resonance frequency of the nucleus of interest, in Hz.
    reference_frequency_hz : float
        Resonance frequency of the reference compound (e.g. TMS), in Hz.

    Returns
    -------
    float or ndarray

    Examples
    --------
    A proton resonating 2800 Hz above TMS on a 400 MHz instrument sits
    at 7.0 ppm:

    >>> round(chemical_shift_ppm(400.0e6 + 2800.0, 400.0e6), 6)
    7.0
    """
    if reference_frequency_hz <= 0:
        raise ValueError("reference_frequency_hz must be positive")
    result = 1e6 * (np.asarray(frequency_hz, dtype=np.float64) - reference_frequency_hz) / reference_frequency_hz
    return float(result) if result.ndim == 0 else result


def karplus_coupling(dihedral_degrees, coefficients=None):
    r"""Vicinal H-C-C-H coupling constant :math:`^3J` from the dihedral angle: the Karplus relation.

    With ``coefficients=None`` this is Karplus's original 1959
    valence-bond result (M. Karplus, *J. Chem. Phys.* 30, 11 (1959)),

    .. math::

        ^3J(\phi) = \begin{cases}
        8.5\cos^2\phi - 0.28, & 0^\circ\le\phi\le90^\circ \\
        9.5\cos^2\phi - 0.28, & 90^\circ<\phi\le180^\circ
        \end{cases}\ \text{Hz};

    otherwise the general three-term form
    :math:`^3J=A\cos^2\phi+B\cos\phi+C` (M. Karplus, *J. Am. Chem.
    Soc.* 85, 2870 (1963)) with ``coefficients=(A, B, C)`` in Hz.

    Parameters
    ----------
    dihedral_degrees : float or array-like of float
        H-C-C-H dihedral angle :math:`\phi`, in degrees.
    coefficients : tuple of float, optional
        ``(A, B, C)`` in Hz for the three-term form.

    Returns
    -------
    float or ndarray
        :math:`^3J`, in Hz.

    Examples
    --------
    Anti (180 degrees) protons couple strongly, gauche (60 degrees)
    weakly, and perpendicular (90 degrees) protons hardly at all:

    >>> round(karplus_coupling(180.0), 2), round(karplus_coupling(60.0), 3), round(karplus_coupling(90.0), 2)
    (9.22, 1.845, -0.28)
    """
    phi = np.radians(np.abs(np.asarray(dihedral_degrees, dtype=np.float64)) % 360.0)
    phi = np.where(phi > np.pi, 2.0 * np.pi - phi, phi)
    cos_phi = np.cos(phi)
    if coefficients is None:
        prefactor = np.where(phi <= np.pi / 2.0 + 1e-12, 8.5, 9.5)
        result = prefactor * cos_phi**2 - 0.28
    else:
        A, B, C = coefficients
        result = A * cos_phi**2 + B * cos_phi + C
    return float(result) if np.ndim(result) == 0 else result


def free_induction_decay(t, offsets_hz, amplitudes, t2: float) -> np.ndarray:
    r"""Complex free-induction decay (FID) after a single 90-degree pulse.

    Each resonance, at offset :math:`\Delta\nu_k` from the rotating-frame
    reference, contributes a decaying complex exponential,

    .. math::

        s(t)=\sum_k a_k\,e^{2\pi i\Delta\nu_k t}\,e^{-t/T_2},

    the time-domain signal whose Fourier transform is the spectrum
    (R. R. Ernst & W. A. Anderson, *Rev. Sci. Instrum.* 37, 93 (1966)).

    Parameters
    ----------
    t : array-like of float
        Sample times, in s.
    offsets_hz : array-like of float
        Resonance offsets, in Hz.
    amplitudes : array-like of float
        Relative amplitude of each resonance.
    t2 : float
        Transverse relaxation time :math:`T_2`, in s (shared by all lines).

    Returns
    -------
    ndarray of complex, shape matching `t`

    Examples
    --------
    At :math:`t=0` every component is in phase, so the signal is the sum
    of the amplitudes:

    >>> fid = free_induction_decay([0.0, 0.1], offsets_hz=[10.0, -25.0], amplitudes=[1.0, 2.0], t2=0.5)
    >>> complex(fid[0])
    (3+0j)
    """
    if t2 <= 0:
        raise ValueError("t2 must be positive")
    t = np.asarray(t, dtype=np.float64)
    offsets = np.asarray(offsets_hz, dtype=np.float64)
    amps = np.asarray(amplitudes, dtype=np.float64)
    if offsets.shape != amps.shape:
        raise ValueError("offsets_hz and amplitudes must have the same shape")
    phases = np.exp(2j * np.pi * offsets[None, :] * t[:, None])
    return (phases @ amps) * np.exp(-t / t2)


def fid_to_spectrum(fid, dwell_time: float) -> tuple:
    r"""Fourier-transform an FID into a frequency-domain (absorption-mode) spectrum.

    Uses the discrete Fourier transform with the first point halved (the
    standard correction that removes a constant baseline offset) and
    returns the real (absorption) part, frequency axis centered on zero.
    A line of transverse relaxation time :math:`T_2` becomes a Lorentzian
    of FWHM :math:`1/(\pi T_2)` Hz.

    Parameters
    ----------
    fid : array-like of complex
        Uniformly sampled FID, starting at :math:`t=0`.
    dwell_time : float
        Sampling interval, in s.

    Returns
    -------
    frequencies_hz : ndarray
        Ascending frequency axis, in Hz.
    spectrum : ndarray
        Real (absorption-mode) spectrum, in units of FID amplitude times seconds.

    Examples
    --------
    One resonance at +50 Hz Fourier-transforms into a single peak at +50 Hz:

    >>> t = np.arange(4096) * 1e-3
    >>> freqs, spec = fid_to_spectrum(free_induction_decay(t, [50.0], [1.0], t2=0.2), 1e-3)
    >>> round(float(freqs[np.argmax(spec)]), 1)
    50.0
    """
    data = np.array(fid, dtype=np.complex128)
    data[0] *= 0.5
    n = data.size
    spectrum = np.fft.fftshift(np.fft.fft(data)) * dwell_time
    frequencies = np.fft.fftshift(np.fft.fftfreq(n, d=dwell_time))
    return frequencies, spectrum.real
