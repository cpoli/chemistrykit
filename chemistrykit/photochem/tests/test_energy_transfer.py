"""Tests for chemistrykit.photochem.systems.energy_transfer."""

import numpy as np
import pytest

from chemistrykit.photochem.systems.energy_transfer import dexter_rate, forster_efficiency, forster_radius, forster_rate


def test_forster_efficiency_half_at_R0_and_rate_consistency():
    R0, tau_D = 5.0, 4.0
    r = np.linspace(2.0, 10.0, 9)
    assert forster_efficiency(R0, R0) == pytest.approx(0.5)
    kT = forster_rate(r, R0, tau_D)
    assert forster_efficiency(r, R0) == pytest.approx(kT / (kT + 1.0 / tau_D))


def test_forster_radius_sixth_root_scaling():
    base = forster_radius(2 / 3, 1.4, 0.5, 1e15)
    assert forster_radius(2 / 3, 1.4, 0.5, 64e15) == pytest.approx(2.0 * base)
    assert base == pytest.approx(0.211 * ((2 / 3) * 1.4**-4 * 0.5 * 1e15) ** (1 / 6))


def test_dexter_rate_exponential():
    r = np.linspace(0.0, 5.0, 6)
    k = dexter_rate(r, K=3.0, J=2.0, L=1.0)
    assert k[0] == pytest.approx(6.0)
    assert np.log(k) == pytest.approx(np.log(6.0) - 2.0 * r)
