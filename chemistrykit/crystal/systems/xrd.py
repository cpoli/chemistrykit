r"""Bragg's law, cubic d-spacings, and structure factors (with systematic absences) for powder XRD.

See West, *Solid State Chemistry and its Applications*, 2nd ed. (2014),
Ch. 5, or Ashcroft & Mermin, *Solid State Physics* (1976), Ch. 6, for
Bragg diffraction and the kinematic structure factor.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations

import numpy as np

__all__ = [
    "bragg_angle",
    "d_spacing_cubic",
    "structure_factor",
    "XRDPeak",
    "powder_xrd_peaks",
    "scherrer_crystallite_size",
]

#: dict: Fractional atomic positions (basis) of the lattice point per
#: conventional cubic unit cell, by Bravais lattice type.
_CUBIC_BASES: dict = {
    "SC": [(0.0, 0.0, 0.0)],
    "BCC": [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
    "FCC": [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)],
}


def bragg_angle(d: float, wavelength: float, order: int = 1):
    r"""Solve Bragg's law :math:`n\lambda=2d\sin\theta` for the diffraction angle :math:`\theta`.

    Parameters
    ----------
    d : float or array-like of float
        Interplanar (`d`-)spacing, in the same length unit as `wavelength`.
    wavelength : float
        X-ray wavelength.
    order : int, default 1
        Diffraction order `n`.

    Returns
    -------
    float or ndarray
        Bragg angle :math:`\theta`, in radians. ``nan`` where
        :math:`n\lambda/(2d) > 1` (no diffraction possible at that
        spacing/wavelength/order).

    Examples
    --------
    Cu-Kalpha (:math:`\lambda=154.18\,\text{pm}`) reflecting off a 3.0 A
    (300 pm) `d`-spacing:

    >>> import numpy as np
    >>> theta = bragg_angle(d=300.0, wavelength=154.18)
    >>> round(float(np.degrees(theta)), 3)
    14.89
    >>> bool(round(float(2.0 * 300.0 * np.sin(theta)), 4) == round(154.18, 4))
    True
    """
    d = np.asarray(d, dtype=np.float64)
    sin_theta = order * wavelength / (2.0 * d)
    sin_theta = np.where(np.abs(sin_theta) <= 1.0, sin_theta, np.nan)
    result = np.arcsin(sin_theta)
    return float(result) if result.ndim == 0 else result


def d_spacing_cubic(a: float, h: int, k: int, l: int):
    r"""Interplanar spacing for the :math:`(hkl)` plane of a cubic lattice, :math:`d=a/\sqrt{h^2+k^2+l^2}`.

    Parameters
    ----------
    a : float
        Cubic lattice constant.
    h, k, l : int
        Miller indices (not all zero).

    Returns
    -------
    float

    Examples
    --------
    >>> round(float(d_spacing_cubic(a=4.0, h=1, k=1, l=1)), 6)
    2.309401
    """
    if h == 0 and k == 0 and l == 0:
        raise ValueError("(h, k, l) = (0, 0, 0) is not a diffracting plane")
    return float(a / np.sqrt(h * h + k * k + l * l))


def structure_factor(hkl, basis, scattering_factors=None) -> complex:
    r"""The kinematic structure factor :math:`F_{hkl}=\sum_j f_j\exp[2\pi i(hx_j+ky_j+lz_j)]`.

    Summed over the atoms in one unit cell at fractional coordinates
    `basis`; :math:`|F_{hkl}|^2` is proportional to a reflection's
    diffracted intensity, and :math:`F_{hkl}=0` identically for certain
    :math:`(hkl)` combinations depending on the lattice centering --
    "systematic absences" (West, *Solid State Chemistry and its
    Applications*, 2nd ed., Ch. 5.4).

    Parameters
    ----------
    hkl : tuple of (int, int, int)
        Miller indices.
    basis : sequence of (float, float, float)
        Fractional atomic coordinates of every atom in the unit cell.
    scattering_factors : sequence of float, optional
        Per-atom scattering factor/amplitude `f_j`, same length as
        `basis`. Defaults to 1.0 for every atom (structure of the
        lattice alone, ignoring atomic form-factor variation).

    Returns
    -------
    complex

    Examples
    --------
    A single atom per cell (simple cubic) never has a systematically
    absent reflection -- :math:`F_{hkl}=f` for every `(hkl)`:

    >>> abs(structure_factor((1, 1, 0), [(0.0, 0.0, 0.0)]))
    1.0

    Body-centered cubic: :math:`F_{hkl}=0` whenever :math:`h+k+l` is odd:

    >>> bcc_basis = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]
    >>> round(abs(structure_factor((1, 0, 0), bcc_basis)), 9)
    0.0
    >>> round(abs(structure_factor((1, 1, 0), bcc_basis)), 9)
    2.0
    """
    h, k, l = hkl
    basis = np.asarray(basis, dtype=np.float64)
    f = np.ones(len(basis)) if scattering_factors is None else np.asarray(scattering_factors, dtype=np.float64)
    phase = 2.0 * np.pi * (h * basis[:, 0] + k * basis[:, 1] + l * basis[:, 2])
    return complex(np.sum(f * np.exp(1j * phase)))


@dataclass
class XRDPeak:
    """A single allowed powder-XRD reflection."""

    hkl: tuple
    """tuple of int: Miller indices :math:`(h,k,l)`."""

    d_spacing: float
    """float: Interplanar spacing, same length unit as the input lattice constant."""

    two_theta: float
    """float: Diffraction angle :math:`2\\theta`, in degrees."""

    relative_intensity: float
    """float: :math:`|F_{hkl}|^2`, not multiplicity- or Lorentz-polarization-corrected."""

    multiplicity: int = 1
    """int: Number of symmetry-equivalent :math:`(hkl)` planes (all sign flips and
    permutations) in the family this peak represents, e.g. 6 for {100}, 12
    for {110}, 8 for {111}; multiply by `relative_intensity` for the
    multiplicity-corrected intensity."""


def powder_xrd_peaks(lattice_type: str, a: float, wavelength: float, hkl_max: int = 3) -> list:
    r"""Enumerate the allowed powder-XRD peaks of a cubic Bravais lattice up to a Miller-index cutoff.

    Returns one peak per family of symmetry-equivalent planes
    :math:`\{hkl\}`, labelled by its canonical member
    :math:`h\ge k\ge l\ge 0` (so BCC's first line is a single (110)
    peak with ``multiplicity=12``, not separate (110)/(101)/(011) peaks).
    Each family's structure factor comes from :func:`structure_factor`
    (systematically-absent families, :math:`|F_{hkl}|^2<10^{-9}`, are
    dropped) and its Bragg angle from :func:`bragg_angle` (families with no
    real solution are dropped); peaks are sorted by :math:`2\theta`.
    Distinct families that happen to share :math:`h^2+k^2+l^2`, and so
    :math:`2\theta` (e.g. {300} and {221}), are kept as separate,
    coincident peaks.

    Parameters
    ----------
    lattice_type : {"SC", "BCC", "FCC"}
        Cubic Bravais lattice type (see :data:`_CUBIC_BASES`).
    a : float
        Cubic lattice constant, in the same length unit as `wavelength`.
    wavelength : float
        X-ray wavelength.
    hkl_max : int, default 3
        Largest Miller index to search. Only canonical
        :math:`h\ge k\ge l\ge 0` indices are enumerated, since
        :math:`d_{hkl}` and :math:`|F_{hkl}|` depend only on
        :math:`h^2+k^2+l^2` and the mixed-parity pattern, both invariant
        under sign flips and permutations of the indices for these
        cubic lattices.

    Returns
    -------
    list of XRDPeak

    Examples
    --------
    BCC iron's first powder line is (110), not (100) -- (100) is
    systematically absent since :math:`1+0+0=1` is odd:

    >>> peaks = powder_xrd_peaks("BCC", a=286.65, wavelength=154.18, hkl_max=2)
    >>> peaks[0].hkl, peaks[0].multiplicity
    ((1, 1, 0), 12)

    FCC's first line is (111) -- (100) and (110) are both absent (mixed
    parity):

    >>> peaks = powder_xrd_peaks("FCC", a=408.6, wavelength=154.18, hkl_max=2)
    >>> peaks[0].hkl
    (1, 1, 1)
    """
    basis = _CUBIC_BASES[lattice_type]
    peaks = []
    for h in range(0, hkl_max + 1):
        for k in range(0, h + 1):
            for l in range(0, k + 1):
                if h == 0:
                    continue
                F = structure_factor((h, k, l), basis)
                intensity = abs(F) ** 2
                if intensity < 1e-9:
                    continue
                d = d_spacing_cubic(a, h, k, l)
                theta = bragg_angle(d, wavelength)
                if np.isnan(theta):
                    continue
                peaks.append(
                    XRDPeak(
                        hkl=(h, k, l),
                        d_spacing=float(d),
                        two_theta=float(np.degrees(2.0 * theta)),
                        relative_intensity=float(intensity),
                        multiplicity=_cubic_multiplicity(h, k, l),
                    )
                )
    # Distinct families that tie exactly on 2*theta (e.g. {300}/{221}) are
    # ordered highest-index-first, purely so results are deterministic.
    peaks.sort(key=lambda p: (p.two_theta, tuple(-x for x in p.hkl)))
    return peaks


def _cubic_multiplicity(h: int, k: int, l: int) -> int:
    """Number of distinct (hkl) related to ``(h, k, l)`` by permutations and sign flips (point group m-3m)."""
    n_permutations = len(set(permutations((h, k, l))))
    n_sign_flips = 2 ** sum(1 for x in (h, k, l) if x != 0)
    return n_permutations * n_sign_flips


def scherrer_crystallite_size(fwhm_deg, two_theta_deg, wavelength: float, shape_factor: float = 0.9):
    r"""Mean crystallite size from powder-XRD line broadening, via the Scherrer equation.

    .. math::

        \tau = \frac{K\lambda}{\beta\cos\theta}

    with :math:`\beta` the peak's full width at half maximum in radians of
    :math:`2\theta` (instrumental broadening already subtracted),
    :math:`\theta` the Bragg angle, and `K` a dimensionless shape factor
    (about 0.9 for roughly spherical crystallites) (P. Scherrer, *Nachr.
    Ges. Wiss. Göttingen*, 1918, 98-100; A. L. Patterson, *Phys. Rev.* 56
    (1939), 978).

    Parameters
    ----------
    fwhm_deg : float or array-like of float
        Full width at half maximum of the peak, in degrees of :math:`2\theta`.
    two_theta_deg : float or array-like of float
        Peak position :math:`2\theta`, in degrees.
    wavelength : float
        X-ray wavelength; the result is in the same length unit.
    shape_factor : float, default 0.9
        Scherrer constant `K`.

    Returns
    -------
    float or ndarray
        Mean crystallite dimension perpendicular to the diffracting planes.

    Examples
    --------
    A Cu-Kalpha peak at :math:`2\theta=38.2°` that is 0.5° wide comes from
    crystallites about 17 nm across:

    >>> round(scherrer_crystallite_size(0.5, 38.2, wavelength=0.15418), 1)
    16.8

    Size is inversely proportional to width -- halve the width, double the size:

    >>> tau_1 = scherrer_crystallite_size(0.5, 38.2, 0.15418)
    >>> tau_2 = scherrer_crystallite_size(0.25, 38.2, 0.15418)
    >>> round(tau_2 / tau_1, 12)
    2.0
    """
    beta = np.radians(np.asarray(fwhm_deg, dtype=np.float64))
    theta = np.radians(np.asarray(two_theta_deg, dtype=np.float64)) / 2.0
    tau = shape_factor * wavelength / (beta * np.cos(theta))
    return float(tau) if tau.ndim == 0 else tau
