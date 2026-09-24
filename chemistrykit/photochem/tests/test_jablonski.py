"""Tests for chemistrykit.photochem.systems.jablonski against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.photochem.systems.jablonski import jablonski_network, jablonski_populations_analytic


def test_jablonski_network_conserves_total_population():
    net = jablonski_network(kf=2.0, kic=1.0, kisc=0.5, kp=0.3, kic_T=0.2, S1_0=1.0)
    result = net.integrate((0.0, 30.0), dt=1e-3, method="rk4")
    total = result.concentration("S1") + result.concentration("T1") + result.concentration("S0")
    np.testing.assert_allclose(total, 1.0, atol=1e-6)


def test_jablonski_populations_analytic_matches_numerical_integration():
    kf, kic, kisc, kp, kic_T, S1_0 = 2.0, 1.0, 0.5, 0.3, 0.2, 1.0
    net = jablonski_network(kf, kic, kisc, kp, kic_T, S1_0)
    result = net.integrate((0.0, 15.0), dt=1e-4, method="rk4")
    S1, T1, S0 = jablonski_populations_analytic(kf, kic, kisc, kp, kic_T, S1_0, result.t)
    np.testing.assert_allclose(result.concentration("S1"), S1, atol=1e-6)
    np.testing.assert_allclose(result.concentration("T1"), T1, atol=1e-6)
    np.testing.assert_allclose(result.concentration("S0"), S0, atol=1e-6)


def test_jablonski_analytic_matches_at_degenerate_rates():
    """When k_S1 == k_T1, consecutive_analytic's degenerate branch is
    exercised; check the network integration still agrees."""
    kf, kic, kisc = 1.0, 0.5, 0.5  # k_S1 = 2.0
    kp, kic_T = 1.2, 0.8  # k_T1 = 2.0, equal to k_S1
    S1_0 = 1.0
    net = jablonski_network(kf, kic, kisc, kp, kic_T, S1_0)
    result = net.integrate((0.0, 10.0), dt=1e-4, method="rk4")
    S1, T1, S0 = jablonski_populations_analytic(kf, kic, kisc, kp, kic_T, S1_0, result.t)
    np.testing.assert_allclose(result.concentration("T1"), T1, atol=1e-6)


def test_s1_decays_as_single_exponential_with_total_rate():
    kf, kic, kisc, kp, kic_T, S1_0 = 3.0, 1.0, 1.0, 0.5, 0.5, 2.0
    t = np.linspace(0.0, 5.0, 50)
    S1, _, _ = jablonski_populations_analytic(kf, kic, kisc, kp, kic_T, S1_0, t)
    k_S1 = kf + kic + kisc
    np.testing.assert_allclose(S1, S1_0 * np.exp(-k_S1 * t), rtol=1e-10)


def test_triplet_population_vanishes_with_no_intersystem_crossing():
    kf, kic, kisc, kp, kic_T, S1_0 = 1.0, 1.0, 0.0, 0.5, 0.5, 1.0
    t = np.linspace(0.0, 10.0, 20)
    _, T1, _ = jablonski_populations_analytic(kf, kic, kisc, kp, kic_T, S1_0, t)
    np.testing.assert_allclose(T1, 0.0, atol=1e-12)


def test_all_population_reaches_ground_state_at_long_time():
    net = jablonski_network(kf=2.0, kic=1.0, kisc=0.5, kp=0.3, kic_T=0.2, S1_0=1.0)
    result = net.integrate((0.0, 100.0), dt=1e-2, method="rk4")
    assert result.concentration("S0")[-1] == pytest.approx(1.0, abs=1e-6)


def test_kasha_emission_yields_branching():
    from chemistrykit.photochem.systems.jablonski import kasha_emission_yields

    phi2, phi1 = kasha_emission_yields(kf2=1.0, k_ic21=3.0, kf1=1.0, knr1=1.0)
    assert phi2 == pytest.approx(0.25)
    assert phi1 == pytest.approx(0.75 * 0.5)
    phi2, phi1 = kasha_emission_yields(kf2=1e8, k_ic21=1e13, kf1=1e8, knr1=1e8)
    assert phi2 == pytest.approx(1e-5, rel=1e-4)
    assert phi1 == pytest.approx(0.5, rel=1e-4)
    assert kasha_emission_yields(1.0, 1.0, 1.0, 3.0, excite="S1") == (0.0, 0.25)
    with pytest.raises(ValueError):
        kasha_emission_yields(1.0, 1.0, 1.0, 1.0, excite="S3")
