"""Tests for chemistrykit.spectro.systems.nmr against Pascal's-triangle/known multiplet patterns."""

import numpy as np
import pytest
import scipy.constants as sc

from chemistrykit.spectro.systems.nmr import (
    chemical_shift_ppm,
    fid_to_spectrum,
    first_order_multiplet,
    free_induction_decay,
    karplus_coupling,
    larmor_frequency,
    multi_coupling_multiplet,
    multiplicity,
    pascals_triangle_intensities,
)


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


def test_larmor_frequency_closed_form():
    gamma_h = sc.physical_constants["proton gyromag. ratio"][0]
    assert larmor_frequency(gamma_h, 11.7) == pytest.approx(gamma_h * 11.7 / (2.0 * np.pi))
    assert larmor_frequency(gamma_h, 2.0 * 11.7) == pytest.approx(2.0 * larmor_frequency(gamma_h, 11.7))


def test_chemical_shift_is_field_independent():
    gamma_h = sc.physical_constants["proton gyromag. ratio"][0]
    shifts = []
    for field in (7.05, 9.4, 14.1):
        ref = larmor_frequency(gamma_h, field)
        shifts.append(chemical_shift_ppm(ref * (1.0 + 3.7e-6), ref))
    assert shifts == pytest.approx([3.7, 3.7, 3.7])


def test_karplus_1959_values_and_general_form():
    assert karplus_coupling(0.0) == pytest.approx(8.5 - 0.28)
    assert karplus_coupling(180.0) == pytest.approx(9.5 - 0.28)
    assert karplus_coupling(90.0) == pytest.approx(-0.28)
    assert karplus_coupling(-60.0) == pytest.approx(karplus_coupling(60.0))
    phi = np.array([0.0, 45.0, 120.0])
    expected = 7.0 * np.cos(np.radians(phi)) ** 2 - 1.0 * np.cos(np.radians(phi)) + 1.5
    assert karplus_coupling(phi, coefficients=(7.0, -1.0, 1.5)) == pytest.approx(expected)


def test_fid_transforms_to_lorentzians_at_offsets():
    dt, n, t2 = 1e-3, 16384, 0.1
    t = np.arange(n) * dt
    fid = free_induction_decay(t, offsets_hz=[-120.0, 80.0], amplitudes=[1.0, 2.0], t2=t2)
    freqs, spec = fid_to_spectrum(fid, dt)
    for offset in (-120.0, 80.0):
        window = np.abs(freqs - offset) < 20.0
        assert freqs[window][np.argmax(spec[window])] == pytest.approx(offset, abs=freqs[1] - freqs[0])
    # Lorentzian absorption peak height a*T2, FWHM 1/(pi*T2).
    i80 = np.argmin(np.abs(freqs - 80.0))
    assert spec[i80] == pytest.approx(2.0 * t2, rel=0.02)
    half = spec[i80] / 2.0
    region = (freqs > 60.0) & (freqs < 100.0) & (spec > half)
    width = freqs[region][-1] - freqs[region][0]
    assert width == pytest.approx(1.0 / (np.pi * t2), abs=2 * (freqs[1] - freqs[0]))


def test_fid_with_seeded_noise_still_recovers_peak():
    rng = np.random.default_rng(0)
    dt, n = 1e-3, 4096
    t = np.arange(n) * dt
    fid = free_induction_decay(t, [50.0], [1.0], t2=0.2) + 0.05 * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    freqs, spec = fid_to_spectrum(fid, dt)
    assert freqs[np.argmax(spec)] == pytest.approx(50.0, abs=0.5)
