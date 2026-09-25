"""Tests for chemistrykit.crystal.systems.xrd: Bragg's law, d-spacings, and combinatorial systematic absences."""

import itertools

import numpy as np
import pytest

from chemistrykit.crystal.systems.xrd import bragg_angle, d_spacing_cubic, powder_xrd_peaks, structure_factor

_SC_BASIS = [(0.0, 0.0, 0.0)]
_BCC_BASIS = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]
_FCC_BASIS = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]


def test_bragg_law_round_trip():
    d = 250.0
    wavelength = 154.18
    theta = bragg_angle(d, wavelength)
    assert 2.0 * d * np.sin(theta) == pytest.approx(wavelength)


def test_bragg_angle_returns_nan_when_no_solution():
    assert np.isnan(bragg_angle(d=10.0, wavelength=154.18))


def test_d_spacing_cubic_matches_pythagorean_formula():
    a = 4.0
    assert d_spacing_cubic(a, 1, 1, 1) == pytest.approx(a / np.sqrt(3.0))
    assert d_spacing_cubic(a, 2, 0, 0) == pytest.approx(a / 2.0)


def test_d_spacing_cubic_rejects_000():
    with pytest.raises(ValueError):
        d_spacing_cubic(4.0, 0, 0, 0)


def test_sc_structure_factor_never_systematically_absent():
    for h, k, l in itertools.product(range(-3, 4), repeat=3):
        if (h, k, l) == (0, 0, 0):
            continue
        assert abs(structure_factor((h, k, l), _SC_BASIS)) == pytest.approx(1.0)


def test_bcc_systematic_absence_combinatorial():
    """F_hkl = 0 iff h+k+l is odd, verified over a broad combinatorial sweep."""
    for h, k, l in itertools.product(range(-4, 5), repeat=3):
        if (h, k, l) == (0, 0, 0):
            continue
        F = structure_factor((h, k, l), _BCC_BASIS)
        odd_sum = (h + k + l) % 2 != 0
        if odd_sum:
            assert abs(F) < 1e-9
        else:
            assert abs(F) == pytest.approx(2.0)


def test_fcc_systematic_absence_combinatorial():
    """F_hkl = 0 iff h, k, l have mixed parity, verified over a broad combinatorial sweep."""
    for h, k, l in itertools.product(range(-4, 5), repeat=3):
        if (h, k, l) == (0, 0, 0):
            continue
        F = structure_factor((h, k, l), _FCC_BASIS)
        mixed_parity = len({h % 2, k % 2, l % 2}) > 1
        if mixed_parity:
            assert abs(F) < 1e-9
        else:
            assert abs(F) == pytest.approx(4.0)


def test_powder_xrd_peaks_bcc_first_line_is_110():
    peaks = powder_xrd_peaks("BCC", a=286.65, wavelength=154.18, hkl_max=2)
    assert peaks[0].hkl == (1, 1, 0)


def test_powder_xrd_peaks_fcc_first_line_is_111():
    peaks = powder_xrd_peaks("FCC", a=408.6, wavelength=154.18, hkl_max=2)
    assert peaks[0].hkl == (1, 1, 1)


@pytest.mark.parametrize("lattice_type", ["SC", "BCC", "FCC"])
def test_powder_xrd_peaks_one_peak_per_hkl_family(lattice_type):
    # Previously (110), (101), and (011) etc. came back as separate, coincident peaks.
    peaks = powder_xrd_peaks(lattice_type, a=408.6, wavelength=154.18, hkl_max=3)
    families = [tuple(sorted(p.hkl, reverse=True)) for p in peaks]
    assert len(families) == len(set(families))
    assert all(p.hkl == family for p, family in zip(peaks, families, strict=True))


def test_powder_xrd_peaks_cubic_multiplicities():
    peaks = {p.hkl: p.multiplicity for p in powder_xrd_peaks("SC", a=408.6, wavelength=154.18, hkl_max=3)}
    assert peaks[(1, 0, 0)] == 6
    assert peaks[(1, 1, 0)] == 12
    assert peaks[(1, 1, 1)] == 8
    assert peaks[(2, 1, 0)] == 24
    assert peaks[(3, 2, 1)] == 48


def test_powder_xrd_peaks_sorted_by_two_theta():
    peaks = powder_xrd_peaks("FCC", a=408.6, wavelength=154.18, hkl_max=3)
    two_thetas = [p.two_theta for p in peaks]
    assert two_thetas == sorted(two_thetas)


def test_powder_xrd_peaks_no_100_or_110_for_fcc():
    peaks = powder_xrd_peaks("FCC", a=408.6, wavelength=154.18, hkl_max=3)
    hkls = {p.hkl for p in peaks}
    assert (1, 0, 0) not in hkls
    assert (1, 1, 0) not in hkls


def test_scherrer_closed_form_and_inverse_width_scaling():
    import numpy as np
    import pytest

    from chemistrykit.crystal.systems.xrd import scherrer_crystallite_size

    fwhm, two_theta, lam = 0.4, 44.0, 0.15418
    expected = 0.9 * lam / (np.radians(fwhm) * np.cos(np.radians(two_theta / 2.0)))
    assert scherrer_crystallite_size(fwhm, two_theta, lam) == pytest.approx(expected)
    sizes = scherrer_crystallite_size(np.array([0.2, 0.4, 0.8]), two_theta, lam)
    np.testing.assert_allclose(sizes[0] / sizes, [1.0, 2.0, 4.0])
