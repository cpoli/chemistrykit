"""Tests for chemistrykit.kinetics.systems.networks against closed-form results."""

import numpy as np
import pytest

from chemistrykit.kinetics.systems.networks import (
    StoichiometricNetwork,
    consecutive_analytic,
    reversible_analytic,
    ssa_intermediate_concentration,
)


def test_parallel_network_conserves_mass():
    net = StoichiometricNetwork.parallel(k1=2.0, k2=0.5, A0=1.0)
    result = net.integrate((0.0, 10.0), dt=1e-3, method="rk4")
    total = result.concentration("A") + result.concentration("B") + result.concentration("C")
    np.testing.assert_allclose(total, 1.0, atol=1e-6)


def test_parallel_network_product_ratio_equals_rate_constant_ratio():
    """For parallel first-order channels A->B (k1), A->C (k2), [B]/[C] = k1/k2 at all t > 0."""
    k1, k2 = 3.0, 0.7
    net = StoichiometricNetwork.parallel(k1=k1, k2=k2, A0=1.0)
    result = net.integrate((0.0, 8.0), dt=1e-3, method="rk4")
    ratio = result.concentration("B")[1:] / result.concentration("C")[1:]
    np.testing.assert_allclose(ratio, k1 / k2, rtol=1e-4)


def test_consecutive_network_matches_bateman_solution():
    k1, k2, A0 = 1.2, 0.4, 1.0
    net = StoichiometricNetwork.consecutive(k1=k1, k2=k2, A0=A0)
    result = net.integrate((0.0, 15.0), dt=1e-3, method="rk4")
    A_exact, B_exact, C_exact = consecutive_analytic(A0, k1, k2, result.t)
    np.testing.assert_allclose(result.concentration("A"), A_exact, atol=1e-6)
    np.testing.assert_allclose(result.concentration("B"), B_exact, atol=1e-6)
    np.testing.assert_allclose(result.concentration("C"), C_exact, atol=1e-6)


def test_consecutive_network_conserves_mass():
    net = StoichiometricNetwork.consecutive(k1=0.9, k2=2.3, A0=1.0)
    result = net.integrate((0.0, 20.0), dt=1e-3, method="rk4")
    total = result.concentration("A") + result.concentration("B") + result.concentration("C")
    np.testing.assert_allclose(total, 1.0, atol=1e-6)


def test_consecutive_analytic_handles_degenerate_equal_rate_constants():
    """The k1 != k2 Bateman formula has a removable singularity at k1 ==
    k2; check the degenerate closed form (A0*k1*t*exp(-k1*t)) matches the
    limit approached from k2 slightly away from k1, and that the network
    integration agrees with it directly."""
    A0, k = 1.0, 1.0
    t = np.linspace(0.1, 10.0, 50)
    _, B_degenerate, _ = consecutive_analytic(A0, k, k, t)
    _, B_near, _ = consecutive_analytic(A0, k, k * 1.0001, t)
    np.testing.assert_allclose(B_degenerate, B_near, rtol=2e-3)

    net = StoichiometricNetwork.consecutive(k1=k, k2=k, A0=A0)
    result = net.integrate((0.0, 10.0), dt=1e-3, method="rk4")
    _, B_expected, _ = consecutive_analytic(A0, k, k, result.t)
    np.testing.assert_allclose(result.concentration("B"), B_expected, atol=1e-6)


def test_steady_state_approximation_improves_as_k2_grows():
    """The SSA estimate of [B] should get closer to the exact Bateman
    solution as k2/k1 grows, per the physical justification in
    ssa_intermediate_concentration's docstring."""
    A0, k1 = 1.0, 1.0
    t = np.linspace(0.5, 10.0, 100)
    errors = []
    for k2 in (2.0, 20.0, 200.0):
        _, B_exact, _ = consecutive_analytic(A0, k1, k2, t)
        B_ssa = ssa_intermediate_concentration(A0, k1, k2, t)
        errors.append(np.max(np.abs(B_exact - B_ssa)))
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 1e-2


def test_reversible_network_matches_analytic_relaxation():
    kf, kr, A0 = 1.5, 0.5, 1.0
    net = StoichiometricNetwork.reversible(kf=kf, kr=kr, A0=A0)
    result = net.integrate((0.0, 10.0), dt=1e-3, method="rk4")
    A_exact, B_exact = reversible_analytic(A0, kf, kr, result.t)
    np.testing.assert_allclose(result.concentration("A"), A_exact, atol=1e-6)
    np.testing.assert_allclose(result.concentration("B"), B_exact, atol=1e-6)


def test_reversible_network_reaches_expected_equilibrium():
    kf, kr, A0 = 3.0, 1.0, 2.0
    net = StoichiometricNetwork.reversible(kf=kf, kr=kr, A0=A0)
    result = net.integrate((0.0, 50.0), dt=1e-3, method="rk4")
    A_eq_expected = kr * A0 / (kf + kr)
    assert result.concentration("A")[-1] == pytest.approx(A_eq_expected, abs=1e-4)


def test_general_constructor_matches_dopri5_and_rk4():
    """Sanity check that the adaptive integrator agrees with fixed-step RK4."""
    net_rk4 = StoichiometricNetwork.consecutive(k1=1.0, k2=5.0, A0=1.0)
    net_adaptive = StoichiometricNetwork.consecutive(k1=1.0, k2=5.0, A0=1.0)
    result_rk4 = net_rk4.integrate((0.0, 5.0), dt=1e-4, method="rk4")
    result_adaptive = net_adaptive.integrate((0.0, 5.0), dt=1e-3, method="dopri5", rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(result_rk4.y[-1], result_adaptive.y[-1], atol=1e-5)
