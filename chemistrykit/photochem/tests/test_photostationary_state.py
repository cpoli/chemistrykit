"""Tests for chemistrykit.photochem.systems.photostationary_state against the exact algebraic solution."""

import pytest

from chemistrykit.photochem.systems.photostationary_state import (
    photostationary_ratio,
    photostationary_state,
    photoswitch_network,
    photoswitch_rate_constants,
)


def test_photostationary_ratio_matches_direct_algebra():
    k_AB, k_BA = 1.5, 0.4
    assert photostationary_ratio(k_AB, k_BA) == pytest.approx(k_AB / k_BA)


def test_photostationary_state_matches_long_time_numerical_integration():
    """The exact algebraic photostationary-state solution should agree
    with direct numerical integration of the two-state ODE system to
    long time -- the required cross-check per this package's convention
    of not trusting closed-form derivations without independent
    numerical verification."""
    k_AB, k_BA, total = 1.5, 0.4, 1.0
    net = photoswitch_network(k_AB, k_BA, A0=total)
    result = net.integrate((0.0, 200.0), dt=1e-2, method="dopri5", rtol=1e-10, atol=1e-12)
    A_numeric, B_numeric = result.concentration("A")[-1], result.concentration("B")[-1]

    pss = photostationary_state(k_AB, k_BA, total_concentration=total)
    assert pss.A_pss == pytest.approx(A_numeric, abs=1e-6)
    assert pss.B_pss == pytest.approx(B_numeric, abs=1e-6)
    assert pss.ratio_B_over_A == pytest.approx(B_numeric / A_numeric, rel=1e-6)


def test_photostationary_state_conserves_total_population():
    pss = photostationary_state(k_AB=2.0, k_BA=0.5, total_concentration=3.0)
    assert pss.A_pss + pss.B_pss == pytest.approx(3.0)


def test_photostationary_ratio_equal_rate_constants_gives_one_to_one():
    pss = photostationary_state(k_AB=1.0, k_BA=1.0)
    assert pss.ratio_B_over_A == pytest.approx(1.0)
    assert pss.A_pss == pytest.approx(pss.B_pss)


def test_photoswitch_rate_constants_ratio_matches_quantum_yield_absorptivity_ratio():
    phi_AB, eps_A, phi_BA, eps_B = 0.6, 2000.0, 0.3, 500.0
    k_AB, k_BA = photoswitch_rate_constants(phi_AB, eps_A, phi_BA, eps_B)
    expected_ratio = (phi_AB * eps_A) / (phi_BA * eps_B)
    assert (k_AB / k_BA) == pytest.approx(expected_ratio)


def test_photoswitch_rate_constants_scale_with_intensity_without_changing_ratio():
    phi_AB, eps_A, phi_BA, eps_B = 0.6, 2000.0, 0.3, 500.0
    k_AB_1, k_BA_1 = photoswitch_rate_constants(phi_AB, eps_A, phi_BA, eps_B, I0=1.0)
    k_AB_2, k_BA_2 = photoswitch_rate_constants(phi_AB, eps_A, phi_BA, eps_B, I0=5.0)
    assert k_AB_2 == pytest.approx(5.0 * k_AB_1)
    assert k_BA_2 == pytest.approx(5.0 * k_BA_1)
    assert (k_AB_1 / k_BA_1) == pytest.approx(k_AB_2 / k_BA_2)


def test_photoswitch_network_approaches_photostationary_ratio():
    k_AB, k_BA = 1.5, 0.4
    net = photoswitch_network(k_AB, k_BA, A0=1.0)
    result = net.integrate((0.0, 100.0), dt=1e-2, method="rk4")
    final_ratio = result.concentration("B")[-1] / result.concentration("A")[-1]
    assert final_ratio == pytest.approx(photostationary_ratio(k_AB, k_BA), rel=1e-4)


def test_photostationary_state_scales_with_total_concentration():
    pss_1 = photostationary_state(k_AB=1.0, k_BA=2.0, total_concentration=1.0)
    pss_2 = photostationary_state(k_AB=1.0, k_BA=2.0, total_concentration=5.0)
    assert pss_2.A_pss == pytest.approx(5.0 * pss_1.A_pss)
    assert pss_2.B_pss == pytest.approx(5.0 * pss_1.B_pss)
    assert pss_2.ratio_B_over_A == pytest.approx(pss_1.ratio_B_over_A)
