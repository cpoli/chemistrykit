"""Tests for chemistrykit.spectro.systems.nmr against Pascal's-triangle/known multiplet patterns."""

import numpy as np
import pytest

from chemistrykit.spectro.systems.nmr import first_order_multiplet, multi_coupling_multiplet, multiplicity, pascals_triangle_intensities


@pytest.mark.parametrize(("n", "expected"), [(0, 1), (1, 2), (2, 3), (3, 4), (6, 7)])
def test_multiplicity_n_plus_one_rule(n, expected):
    assert multiplicity(n) == expected


def test_multiplicity_rejects_negative():
    with pytest.raises(ValueError):
        multiplicity(-1)


@pytest.mark.parametrize(
    ("n", "expected"),
    [
        (0, [1.0]),
        (1, [1.0, 1.0]),
        (2, [1.0, 2.0, 1.0]),
        (3, [1.0, 3.0, 3.0, 1.0]),
        (4, [1.0, 4.0, 6.0, 4.0, 1.0]),
    ],
)
def test_pascals_triangle_rows(n, expected):
    assert pascals_triangle_intensities(n).tolist() == pytest.approx(expected)


def test_first_order_multiplet_line_count():
    spectrum = first_order_multiplet(chemical_shift_ppm=2.5, j_coupling_hz=7.5, n_neighbors=3, spectrometer_frequency_mhz=300.0)
    assert len(spectrum.positions) == 4


def test_first_order_multiplet_symmetric_about_shift():
    spectrum = first_order_multiplet(chemical_shift_ppm=2.5, j_coupling_hz=7.5, n_neighbors=4, spectrometer_frequency_mhz=300.0)
    assert np.mean(spectrum.positions) == pytest.approx(2.5)


def test_first_order_multiplet_spacing_matches_j_in_ppm():
    freq_mhz = 400.0
    j_hz = 7.0
    spectrum = first_order_multiplet(chemical_shift_ppm=1.0, j_coupling_hz=j_hz, n_neighbors=2, spectrometer_frequency_mhz=freq_mhz)
    spacing_ppm = np.diff(spectrum.positions)
    assert spacing_ppm == pytest.approx((j_hz / freq_mhz) * np.ones_like(spacing_ppm))


def test_first_order_multiplet_intensities_are_binomial():
    spectrum = first_order_multiplet(chemical_shift_ppm=1.0, j_coupling_hz=7.0, n_neighbors=4, spectrometer_frequency_mhz=400.0)
    assert spectrum.intensities.tolist() == pytest.approx([1.0, 4.0, 6.0, 4.0, 1.0])


def test_multi_coupling_equal_j_degenerates_to_larger_equivalent_set():
    combined = multi_coupling_multiplet(chemical_shift_ppm=5.0, couplings=[(7.0, 1), (7.0, 1)], spectrometer_frequency_mhz=400.0)
    equivalent = first_order_multiplet(chemical_shift_ppm=5.0, j_coupling_hz=7.0, n_neighbors=2, spectrometer_frequency_mhz=400.0)
    assert combined.positions == pytest.approx(equivalent.positions)
    assert combined.intensities == pytest.approx(equivalent.intensities)


def test_multi_coupling_doublet_of_triplets_has_six_lines_and_conserves_total_intensity():
    dt = multi_coupling_multiplet(chemical_shift_ppm=3.0, couplings=[(12.0, 1), (5.0, 2)], spectrometer_frequency_mhz=400.0)
    assert len(dt.positions) == 6
    assert np.sum(dt.intensities) == pytest.approx(2.0 * 4.0)


def test_multi_coupling_no_couplings_gives_single_line():
    singlet = multi_coupling_multiplet(chemical_shift_ppm=7.26, couplings=[], spectrometer_frequency_mhz=400.0)
    assert len(singlet.positions) == 1
    assert singlet.positions[0] == pytest.approx(7.26)
    assert singlet.intensities[0] == pytest.approx(1.0)
