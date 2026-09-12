r"""Fluorescence/phosphorescence quantum yields, and the photochemical quantum yield via Beer-Lambert.

See Lakowicz, *Principles of Fluorescence Spectroscopy*, 3rd ed., Ch. 1,
for the radiative/total-decay-rate quantum yields, and Turro, Ramamurthy
& Scaiano, *Modern Molecular Photochemistry of Organic Molecules* (2010),
Ch. 7, for the photochemical (photon-in, product-out) quantum yield.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.spectro.systems.beer_lambert import transmittance

__all__ = [
    "fluorescence_quantum_yield",
    "intersystem_crossing_yield",
    "phosphorescence_quantum_yield",
    "photons_absorbed",
    "photochemical_quantum_yield",
]


def fluorescence_quantum_yield(kf: float, kic: float, kisc: float) -> float:
    r"""Fluorescence quantum yield :math:`\Phi_f = k_f/(k_f+k_{ic}+k_{isc})`.

    The fraction of excited :math:`S_1` molecules that relax by emitting
    a photon (fluorescence) rather than through either nonradiative decay
    channel -- the ratio of the radiative rate constant to the *total*
    :math:`S_1` decay rate (Lakowicz, *Principles of Fluorescence
    Spectroscopy*, 3rd ed., Ch. 1, eq. 1.1). See
    :mod:`chemistrykit.photochem.systems.jablonski` for the underlying
    kinetic model this ratio comes from.

    Parameters
    ----------
    kf : float
        Fluorescence (radiative) rate constant.
    kic : float
        Internal conversion (nonradiative, :math:`S_1\to S_0`) rate
        constant.
    kisc : float
        Intersystem crossing (:math:`S_1\to T_1`) rate constant.

    Returns
    -------
    float
        In :math:`[0, 1]`.

    Examples
    --------
    With no competing nonradiative pathway, every excited molecule
    fluoresces:

    >>> round(fluorescence_quantum_yield(kf=1.0, kic=0.0, kisc=0.0), 6)
    1.0

    Adding a nonradiative channel of the same magnitude as `kf` exactly
    halves the yield:

    >>> round(fluorescence_quantum_yield(kf=1.0, kic=1.0, kisc=0.0), 6)
    0.5
    """
    return kf / (kf + kic + kisc)


def intersystem_crossing_yield(kisc: float, kf: float, kic: float) -> float:
    r"""Intersystem-crossing yield :math:`\Phi_{isc} = k_{isc}/(k_f+k_{ic}+k_{isc})`.

    The fraction of excited :math:`S_1` population that crosses to the
    triplet manifold :math:`T_1` rather than returning directly to
    :math:`S_0` (Turro, Ramamurthy & Scaiano, *Modern Molecular
    Photochemistry of Organic Molecules*, Ch. 5).

    Parameters
    ----------
    kisc : float
        Intersystem crossing rate constant.
    kf : float
        Fluorescence rate constant.
    kic : float
        Internal conversion rate constant.

    Returns
    -------
    float
        In :math:`[0, 1]`.

    Examples
    --------
    Together with :func:`fluorescence_quantum_yield` and the (implicit)
    internal-conversion yield, the three :math:`S_1` decay-channel yields
    sum to exactly 1:

    >>> kf, kic, kisc = 2.0, 1.0, 0.5
    >>> phi_f = fluorescence_quantum_yield(kf, kic, kisc)
    >>> phi_isc = intersystem_crossing_yield(kisc, kf, kic)
    >>> phi_ic = kic / (kf + kic + kisc)
    >>> round(phi_f + phi_isc + phi_ic, 9)
    1.0
    """
    return kisc / (kf + kic + kisc)


def phosphorescence_quantum_yield(kisc: float, kf: float, kic: float, kp: float, kic_T: float) -> float:
    r"""Phosphorescence quantum yield :math:`\Phi_p = \Phi_{isc}\cdot\frac{k_p}{k_p+k_{ic,T}}`.

    Phosphorescence requires *two* successive branching events to
    succeed: the molecule must first cross to the triplet manifold
    (probability :math:`\Phi_{isc}`), and the resulting triplet must then
    decay radiatively rather than nonradiatively (probability
    :math:`k_p/(k_p+k_{ic,T})`) -- the product of these two independent
    branching ratios (Lakowicz, *Principles of Fluorescence
    Spectroscopy*, 3rd ed., Ch. 1; Turro, Ramamurthy & Scaiano, *Modern
    Molecular Photochemistry of Organic Molecules*, Ch. 5).

    Parameters
    ----------
    kisc : float
        Intersystem crossing rate constant.
    kf : float
        Fluorescence rate constant.
    kic : float
        :math:`S_1` internal conversion rate constant.
    kp : float
        Phosphorescence (radiative :math:`T_1\to S_0`) rate constant.
    kic_T : float
        Triplet nonradiative decay rate constant.

    Returns
    -------
    float
        In :math:`[0, 1]`.

    Examples
    --------
    With no intersystem crossing at all, no phosphorescence is possible:

    >>> round(phosphorescence_quantum_yield(kisc=0.0, kf=1.0, kic=0.0, kp=1.0, kic_T=0.0), 6)
    0.0

    With every :math:`S_1` crossing to the triplet, and every triplet
    decaying radiatively, the phosphorescence yield is exactly 1:

    >>> round(phosphorescence_quantum_yield(kisc=1.0, kf=0.0, kic=0.0, kp=1.0, kic_T=0.0), 6)
    1.0
    """
    phi_isc = intersystem_crossing_yield(kisc, kf, kic)
    return phi_isc * (kp / (kp + kic_T))


def photons_absorbed(photon_flux_incident, absorbance):
    r"""Photon flux (or count) actually absorbed by a sample, from the Beer-Lambert absorbance.

    :math:`I_{abs} = I_0(1-10^{-A}) = I_0(1-T)`, using
    :func:`chemistrykit.spectro.systems.beer_lambert.transmittance` for
    :math:`T=10^{-A}` -- the fraction of incident light *not*
    transmitted is, by energy conservation (neglecting reflection and
    scattering losses), the fraction absorbed (Turro, Ramamurthy &
    Scaiano, *Modern Molecular Photochemistry of Organic Molecules*, Ch.
    7).

    Parameters
    ----------
    photon_flux_incident : float or array-like of float
        Incident photon flux (or photon count over an exposure time), in
        any consistent unit (e.g. mol photons / s, or einstein/s).
    absorbance : float or array-like of float
        Beer-Lambert absorbance of the sample at the excitation
        wavelength (see :mod:`chemistrykit.spectro.systems.beer_lambert`).

    Returns
    -------
    float or ndarray
        Absorbed photon flux, same unit as `photon_flux_incident`.

    Examples
    --------
    At zero absorbance, nothing is absorbed:

    >>> round(float(photons_absorbed(1.0, absorbance=0.0)), 6)
    0.0

    At very high absorbance, essentially all incident light is absorbed:

    >>> round(float(photons_absorbed(1.0, absorbance=10.0)), 6)
    1.0
    """
    photon_flux_incident = np.asarray(photon_flux_incident, dtype=np.float64)
    result = photon_flux_incident * (1.0 - transmittance(absorbance))
    return float(result) if result.ndim == 0 else result


def photochemical_quantum_yield(moles_product_formed, moles_photons_absorbed):
    r"""Photochemical quantum yield :math:`\Phi = \frac{\text{moles of product formed}}{\text{moles of photons absorbed}}`.

    The photochemistry analogue of the excited-state quantum yields
    above: the efficiency of converting absorbed photons into chemical
    product (Turro, Ramamurthy & Scaiano, *Modern Molecular
    Photochemistry of Organic Molecules*, Ch. 7). **Note on range**: for
    a simple, non-chain photoreaction :math:`0\le\Phi\le1` (each absorbed
    photon produces at most one product molecule), but a radical-chain
    photoreaction can propagate after the initiating absorption event and
    give :math:`\Phi\gg1` -- this function does not itself enforce an
    upper bound, since chain reactions are a real and common exception,
    not implemented here.

    Parameters
    ----------
    moles_product_formed : float or array-like of float
    moles_photons_absorbed : float or array-like of float

    Returns
    -------
    float or ndarray

    Examples
    --------
    One photon absorbed producing exactly one product molecule (a
    typical simple photoisomerization) gives :math:`\Phi=1`:

    >>> round(float(photochemical_quantum_yield(1.0, 1.0)), 6)
    1.0
    """
    result = np.asarray(moles_product_formed, dtype=np.float64) / np.asarray(moles_photons_absorbed, dtype=np.float64)
    return float(result) if result.ndim == 0 else result
