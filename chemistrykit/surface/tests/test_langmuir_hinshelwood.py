"""Tests for chemistrykit.surface.systems.langmuir_hinshelwood against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.surface.systems.langmuir_hinshelwood import lh_rate_dual_site, lh_rate_single_site


def test_lh_single_site_half_at_half_saturation_pressure():
    k, K_A = 4.0, 2.0
    assert lh_rate_single_site(k, K_A, P_A=1.0 / K_A) == pytest.approx(k / 2.0)


def test_lh_single_site_saturates_to_k_at_high_pressure():
    assert lh_rate_single_site(k=5.0, K_A=1.0, P_A=1.0e8) == pytest.approx(5.0, rel=1e-6)


def test_lh_single_site_low_pressure_is_first_order():
    k, K_A = 4.0, 2.0
    P = np.array([1e-4, 2e-4, 4e-4])
    rate = lh_rate_single_site(k, K_A, P)
    np.testing.assert_allclose(rate, k * K_A * P, rtol=1e-3)


def test_lh_dual_site_vanishes_without_either_reactant():
    assert lh_rate_dual_site(k=4.0, K_A=2.0, P_A=0.0, K_B=5.0, P_B=1.0) == pytest.approx(0.0)
    assert lh_rate_dual_site(k=4.0, K_A=2.0, P_A=1.0, K_B=5.0, P_B=0.0) == pytest.approx(0.0)


def test_lh_dual_site_low_coverage_matches_mass_action_product():
    k, K_A, K_B = 4.0, 2.0, 5.0
    P_A = np.array([1e-4, 2e-4])
    P_B = np.array([1e-4, 3e-4])
    rate = lh_rate_dual_site(k, K_A, P_A, K_B, P_B)
    np.testing.assert_allclose(rate, k * K_A * P_A * K_B * P_B, rtol=1e-2)


def test_lh_dual_site_rate_is_symmetric_in_A_and_B():
    k, K_A, K_B, P_A, P_B = 3.0, 2.0, 5.0, 0.3, 0.7
    rate_AB = lh_rate_dual_site(k, K_A, P_A, K_B, P_B)
    rate_BA = lh_rate_dual_site(k, K_B, P_B, K_A, P_A)
    assert rate_AB == pytest.approx(rate_BA)


def test_lh_dual_site_rate_is_nonmonotonic_in_P_A():
    P_A = np.linspace(0.01, 50.0, 200)
    rate = lh_rate_dual_site(k=4.0, K_A=2.0, P_A=P_A, K_B=5.0, P_B=0.2)
    peak = int(np.argmax(rate))
    assert 0 < peak < len(P_A) - 1
