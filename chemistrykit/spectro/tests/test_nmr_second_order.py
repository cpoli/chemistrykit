"""Tests for exact second-order NMR spectra (AB, ABX, general spin-1/2 systems)."""

import numpy as np
import pytest

from chemistrykit.spectro.systems.nmr import ab_quartet, abx_spectrum, multi_coupling_multiplet, second_order_spectrum


@pytest.mark.parametrize("delta_ppm, J", [(0.5, 7.0), (0.05, 10.0), (0.01, 15.0)])
def test_exact_two_spin_matches_ab_formula(delta_ppm, J):
    exact = second_order_spectrum([1.0, 1.0 + delta_ppm], [[0, J], [J, 0]], 400.0)
    analytic = ab_quartet(1.0, 1.0 + delta_ppm, J, 400.0)
    np.testing.assert_allclose(exact.positions, analytic.positions, atol=1e-12)
    np.testing.assert_allclose(exact.intensities, analytic.intensities, atol=1e-10)


def test_equivalent_spins_give_a_single_line():
    s = second_order_spectrum([2.0, 2.0], [[0, 7.0], [7.0, 0]], 400.0)
    assert s.positions == pytest.approx([2.0])
    assert s.intensities == pytest.approx([4.0])  # N * 2**(N-1), same total as two uncoupled spins


def test_weak_coupling_limit_reproduces_first_order_multiplet():
    s = abx_spectrum([1.0, 3.0, 5.0], j_ab_hz=4.0, j_ax_hz=8.0, j_bx_hz=0.0, spectrometer_frequency_mhz=60000.0)
    first_order = multi_coupling_multiplet(1.0, [(4.0, 1), (8.0, 1)], 60000.0)
    a_lines = s.positions < 2.0
    np.testing.assert_allclose(np.sort(s.positions[a_lines]), np.sort(first_order.positions), atol=1e-9)
    np.testing.assert_allclose(s.intensities[a_lines], 1.0, atol=1e-3)


def test_abx_total_intensity_and_symmetry():
    s = abx_spectrum([2.50, 2.55, 4.80], 16.0, 5.0, 8.0, 400.0)
    assert s.intensities.sum() == pytest.approx(12.0)
    assert np.all(np.diff(s.positions) > 0)
    x_part = s.positions > 4.0
    assert s.intensities[x_part].sum() == pytest.approx(4.0, abs=0.1)


def test_bad_coupling_shape_raises():
    with pytest.raises(ValueError):
        second_order_spectrum([1.0, 2.0], [[0.0]], 400.0)
