"""Tests for chemistrykit.spectro.core.base_system.Spectrum."""

import numpy as np
import pytest

from chemistrykit.spectro.core.base_system import Spectrum


def test_rejects_mismatched_shapes():
    with pytest.raises(ValueError):
        Spectrum(positions=np.array([1.0, 2.0]), intensities=np.array([1.0]))


def test_broaden_integral_equals_sum_of_intensities():
    spectrum = Spectrum(positions=np.array([10.0, 20.0, 30.0]), intensities=np.array([1.0, 2.0, 0.5]))
    x = np.linspace(0.0, 40.0, 40001)
    broadened = spectrum.broaden(x, shape="gaussian", fwhm=1.0)
    assert np.trapezoid(broadened, x) == pytest.approx(3.5, abs=1e-2)


def test_normalized_scales_to_unit_max():
    spectrum = Spectrum(positions=np.array([1.0, 2.0, 3.0]), intensities=np.array([2.0, 8.0, 4.0]))
    normalized = spectrum.normalized()
    assert np.max(normalized.intensities) == pytest.approx(1.0)
    assert normalized.intensities.tolist() == pytest.approx([0.25, 1.0, 0.5])


def test_normalized_does_not_mutate_original():
    spectrum = Spectrum(positions=np.array([1.0]), intensities=np.array([5.0]))
    spectrum.normalized()
    assert spectrum.intensities[0] == pytest.approx(5.0)
