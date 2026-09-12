"""Tests for chemistrykit.analytical.systems.calibration: LOD/LOQ per IUPAC convention."""

import numpy as np
import pytest

from chemistrykit.analytical.systems.calibration import fit_calibration


def test_calibration_recovers_known_slope_and_intercept_from_exact_data():
    conc = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    signal = 2.5 * conc + 0.3
    cal = fit_calibration(conc, signal)
    assert cal.slope == pytest.approx(2.5, rel=1e-9)
    assert cal.intercept == pytest.approx(0.3, rel=1e-9)
    assert cal.r_squared == pytest.approx(1.0, abs=1e-9)


def test_calibration_exact_data_has_zero_lod_and_loq():
    conc = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    signal = 2.5 * conc + 0.3
    cal = fit_calibration(conc, signal)
    assert cal.lod() == pytest.approx(0.0, abs=1e-9)
    assert cal.loq() == pytest.approx(0.0, abs=1e-9)


def test_loq_to_lod_ratio_is_exactly_10_over_3_3():
    rng = np.random.default_rng(1)
    conc = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    signal = 2.0 * conc + 0.1 + rng.normal(scale=0.05, size=conc.shape)
    cal = fit_calibration(conc, signal)
    assert cal.loq() / cal.lod() == pytest.approx(10.0 / 3.3, rel=1e-9)


def test_predict_signal_and_predict_concentration_are_inverses():
    conc = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    signal = 2.5 * conc + 0.3
    cal = fit_calibration(conc, signal)
    test_conc = 2.7
    predicted_signal = cal.predict_signal(test_conc)
    recovered_conc = cal.predict_concentration(predicted_signal)
    assert recovered_conc == pytest.approx(test_conc, rel=1e-9)


def test_noisier_data_gives_larger_lod():
    rng = np.random.default_rng(2)
    conc = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    low_noise = fit_calibration(conc, 2.0 * conc + 0.1 + rng.normal(scale=0.01, size=conc.shape))
    high_noise = fit_calibration(conc, 2.0 * conc + 0.1 + rng.normal(scale=0.5, size=conc.shape))
    assert high_noise.lod() > low_noise.lod()
