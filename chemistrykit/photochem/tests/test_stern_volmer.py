"""Tests for chemistrykit.photochem.systems.stern_volmer against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.photochem.systems.stern_volmer import (
    classify_quenching_mechanism,
    dynamic_quenching_constant,
    fit_stern_volmer,
    stern_volmer_ratio,
)


def test_stern_volmer_ratio_is_one_at_zero_quencher():
    assert stern_volmer_ratio(Ksv=25.0, Q=0.0) == pytest.approx(1.0)


def test_stern_volmer_ratio_is_linear_in_quencher_concentration():
    Ksv = 30.0
    Q = np.linspace(0.0, 0.05, 20)
    ratio = stern_volmer_ratio(Ksv, Q)
    # Linear: constant first difference / constant slope.
    diffs = np.diff(ratio)
    np.testing.assert_allclose(diffs, diffs[0])
    slope = np.polyfit(Q, ratio, 1)[0]
    assert slope == pytest.approx(Ksv)


def test_dynamic_quenching_constant_is_product_of_rate_and_lifetime():
    assert dynamic_quenching_constant(kq=2.0e9, tau0=1.0e-8) == pytest.approx(20.0)


def test_fit_stern_volmer_recovers_known_ksv():
    Ksv_true = 42.0
    Q = np.array([0.0, 0.002, 0.004, 0.008, 0.016])
    ratio = stern_volmer_ratio(Ksv_true, Q)
    fit = fit_stern_volmer(Q, ratio)
    assert fit.Ksv == pytest.approx(Ksv_true, rel=1e-6)
    assert fit.r_squared == pytest.approx(1.0, abs=1e-9)


def test_classify_quenching_mechanism_dynamic():
    assert classify_quenching_mechanism(intensity_ratio_slope=40.0, lifetime_ratio_slope=40.0) == "dynamic"


def test_classify_quenching_mechanism_static():
    assert classify_quenching_mechanism(intensity_ratio_slope=40.0, lifetime_ratio_slope=0.0) == "static"


def test_classify_quenching_mechanism_mixed():
    assert classify_quenching_mechanism(intensity_ratio_slope=40.0, lifetime_ratio_slope=20.0) == "mixed"
