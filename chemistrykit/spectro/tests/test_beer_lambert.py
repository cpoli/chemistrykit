"""Tests for chemistrykit.spectro.systems.beer_lambert against closed-form results."""

import numpy as np
import pytest

from chemistrykit.spectro.systems.beer_lambert import (
    absorbance,
    apparent_absorbance_with_stray_light,
    concentration_from_absorbance,
    transmittance,
)


def test_absorbance_linear_in_concentration():
    a1 = absorbance(5000.0, 1.0e-5, 1.0)
    a2 = absorbance(5000.0, 2.0e-5, 1.0)
    assert a2 == pytest.approx(2.0 * a1)


def test_absorbance_linear_in_path_length():
    a1 = absorbance(5000.0, 1.0e-5, 1.0)
    a2 = absorbance(5000.0, 1.0e-5, 2.0)
    assert a2 == pytest.approx(2.0 * a1)


def test_transmittance_absorbance_are_consistent():
    a = 0.301
    t = transmittance(a)
    assert t == pytest.approx(10.0**-0.301)


def test_transmittance_zero_absorbance_is_full_transmission():
    assert transmittance(0.0) == pytest.approx(1.0)


def test_concentration_from_absorbance_is_inverse():
    c = 3.7e-5
    a = absorbance(5000.0, c, 1.0)
    assert concentration_from_absorbance(a, 5000.0, 1.0) == pytest.approx(c)


def test_apparent_absorbance_equals_true_absorbance_with_no_stray_light():
    a_true = absorbance(5000.0, 5.0e-5, 1.0)
    a_apparent = apparent_absorbance_with_stray_light(5000.0, 5.0e-5, 1.0, stray_light_fraction=0.0)
    assert a_apparent == pytest.approx(a_true)


def test_apparent_absorbance_deviates_more_at_high_concentration():
    kwargs = dict(molar_absorptivity=5000.0, path_length=1.0, stray_light_fraction=0.001)
    low_true = absorbance(5000.0, 1.0e-6, 1.0)
    low_apparent = apparent_absorbance_with_stray_light(concentration=1.0e-6, **kwargs)
    high_true = absorbance(5000.0, 1.0e-3, 1.0)
    high_apparent = apparent_absorbance_with_stray_light(concentration=1.0e-3, **kwargs)
    assert abs(high_apparent - high_true) > abs(low_apparent - low_true)


def test_apparent_absorbance_never_exceeds_true_absorbance():
    a_true = absorbance(5000.0, 2.0e-4, 1.0)
    a_apparent = apparent_absorbance_with_stray_light(5000.0, 2.0e-4, 1.0, stray_light_fraction=0.002)
    assert a_apparent <= a_true


def test_apparent_absorbance_of_blank_is_zero_with_stray_light():
    assert apparent_absorbance_with_stray_light(5000.0, 0.0, 1.0, stray_light_fraction=0.01) == pytest.approx(0.0, abs=1e-15)


def test_apparent_absorbance_saturates_at_stray_light_limit():
    # As the sample becomes opaque, only stray light reaches the detector: A -> log10((1 + s)/s).
    s = 0.01
    a_apparent = apparent_absorbance_with_stray_light(5000.0, 1.0, 1.0, stray_light_fraction=s)
    assert a_apparent == pytest.approx(np.log10((1.0 + s) / s), rel=1e-12)


def test_apparent_absorbance_rejects_invalid_stray_light_fraction():
    with pytest.raises(ValueError):
        apparent_absorbance_with_stray_light(5000.0, 1.0e-5, 1.0, stray_light_fraction=1.0)
    with pytest.raises(ValueError):
        apparent_absorbance_with_stray_light(5000.0, 1.0e-5, 1.0, stray_light_fraction=-0.1)
