r"""Bragg's law, cubic d-spacings, and structure factors (with systematic absences) for powder XRD.

See West, *Solid State Chemistry and its Applications*, 2nd ed. (2014),
Ch. 5, or Ashcroft & Mermin, *Solid State Physics* (1976), Ch. 6, for
Bragg diffraction and the kinematic structure factor.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "bragg_angle",
    "d_spacing_cubic",
    "structure_factor",
    "XRDPeak",
    "powder_xrd_peaks",
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


def powder_xrd_peaks(lattice_type: str, a: float, wavelength: float, hkl_max: int = 3) -> list:
    r"""Enumerate the allowed powder-XRD peaks of a cubic Bravais lattice up to a Miller-index cutoff.

    Iterates every :math:`(h,k,l)` with :math:`0\le h,k,l\le` `hkl_max`
    (not all zero), computes each reflection's structure factor via
    :func:`structure_factor` (dropping systematically-absent
    reflections, :math:`|F_{hkl}|^2<10^{-9}`) and Bragg angle via
    :func:`bragg_angle` (dropping reflections with no real solution),
    and returns the surviving peaks sorted by :math:`2\theta`.

    Parameters
    ----------
    lattice_type : {"SC", "BCC", "FCC"}
        Cubic Bravais lattice type (see :data:`_CUBIC_BASES`).
    a : float
        Cubic lattice constant, in the same length unit as `wavelength`.
    wavelength : float
        X-ray wavelength.
    hkl_max : int, default 3
        Largest Miller index to search (non-negative indices only, since
        :math:`d_{hkl}` and :math:`|F_{hkl}|` depend only on
        :math:`h^2+k^2+l^2` and the mixed-parity pattern, both invariant
        under sign flips of any index).

    Returns
    -------
    list of XRDPeak

    Examples
    --------
    BCC iron's first powder line is (110), not (100) -- (100) is
    systematically absent since :math:`1+0+0=1` is odd:

    >>> peaks = powder_xrd_peaks("BCC", a=286.65, wavelength=154.18, hkl_max=2)
    >>> peaks[0].hkl
    (1, 1, 0)

    FCC's first line is (111) -- (100) and (110) are both absent (mixed
    parity):

    >>> peaks = powder_xrd_peaks("FCC", a=408.6, wavelength=154.18, hkl_max=2)
    >>> peaks[0].hkl
    (1, 1, 1)
    """
    basis = _CUBIC_BASES[lattice_type]
    peaks = []
    for h in range(0, hkl_max + 1):
        for k in range(0, hkl_max + 1):
            for l in range(0, hkl_max + 1):
                if h == 0 and k == 0 and l == 0:
                    continue
                F = structure_factor((h, k, l), basis)
                intensity = abs(F) ** 2
                if intensity < 1e-9:
                    continue
                d = d_spacing_cubic(a, h, k, l)
                theta = bragg_angle(d, wavelength)
                if np.isnan(theta):
                    continue
                peaks.append(XRDPeak(hkl=(h, k, l), d_spacing=float(d), two_theta=float(np.degrees(2.0 * theta)), relative_intensity=float(intensity)))
    # Peaks that tie exactly on 2*theta (e.g. BCC's (110)/(101)/(011), all
    # equivalent by cubic symmetry) are broken by a canonical
    # highest-index-first ordering, purely so results are deterministic.
    peaks.sort(key=lambda p: (p.two_theta, tuple(-x for x in p.hkl)))
    return peaks
