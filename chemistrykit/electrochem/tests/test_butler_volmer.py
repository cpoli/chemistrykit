"""Tests for chemistrykit.electrochem.systems.butler_volmer against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.electrochem.systems.butler_volmer import (
    butler_volmer_current_density,
    exchange_current_density,
    fit_tafel_plot,
    tafel_overpotential,
    tafel_slope,
)


def test_zero_overpotential_gives_zero_net_current():
    assert butler_volmer_current_density(i0=1e-6, eta=0.0) == pytest.approx(0.0, abs=1e-15)


def test_current_direction_matches_overpotential_sign():
    assert butler_volmer_current_density(1e-6, eta=0.1) > 0.0
    assert butler_volmer_current_density(1e-6, eta=-0.1) < 0.0


def test_tafel_linearization_converges_to_full_equation_at_high_overpotential():
    """The Tafel approximation drops the Butler-Volmer back-reaction term,
    whose magnitude falls off exponentially as eta grows; the error
    should be large close to equilibrium and become small (<1%) by
    ~150-200 mV. (The error crosses zero somewhere in between as the
    approximation switches from under- to over-estimating, so it is not
    strictly monotonic pointwise -- only its overall envelope shrinks.)"""
    i0, alpha, n = 1e-6, 0.5, 1
    etas = [0.05, 0.10, 0.15, 0.20, 0.30]
    rel_errors = []
    for eta in etas:
        i_full = butler_volmer_current_density(i0, eta, alpha=alpha, n=n)
        eta_tafel = tafel_overpotential(i_full, i0, alpha=alpha, n=n, branch="anodic")
        rel_errors.append(abs(eta_tafel - eta) / eta)
    # The error at 50 mV (close to equilibrium) is much larger than at
    # 150-300 mV (solidly in the Tafel regime), and every high-eta point
    # is small.
    assert rel_errors[0] > 10.0 * rel_errors[2]
    assert all(err < 0.005 for err in rel_errors[2:])


def test_tafel_cathodic_branch_matches_full_equation_at_high_overpotential():
    i0, alpha, n = 2e-7, 0.4, 2
    eta = -0.35
    i_full = butler_volmer_current_density(i0, eta, alpha=alpha, n=n)
    eta_tafel = tafel_overpotential(i_full, i0, alpha=alpha, n=n, branch="cathodic")
    assert eta_tafel == pytest.approx(eta, rel=1e-3)


def test_tafel_slope_signs_and_symmetry():
    b_a = tafel_slope(alpha=0.5, n=1, branch="anodic")
    b_c = tafel_slope(alpha=0.5, n=1, branch="cathodic")
    assert b_a > 0.0
    assert b_c < 0.0
    assert b_a == pytest.approx(-b_c)  # symmetric charge-transfer coefficient


def test_tafel_slope_numeric_derivative_matches_formula():
    """b = d(eta)/d(log10 i) in the Tafel (high-overpotential) regime,
    verified by numerical differentiation of the full nonlinear BV curve."""
    i0, alpha, n = 1e-6, 0.35, 1
    eta = np.linspace(0.25, 0.35, 5)
    i = butler_volmer_current_density(i0, eta, alpha=alpha, n=n)
    numeric_slope = np.polyfit(np.log10(i), eta, 1)[0]
    assert numeric_slope == pytest.approx(tafel_slope(alpha=alpha, n=n, branch="anodic"), rel=1e-3)


def test_exchange_current_density_scales_linearly_with_rate_constant():
    i0_1 = exchange_current_density(k0=1e-5, C_ox=1.0, C_red=1.0)
    i0_2 = exchange_current_density(k0=2e-5, C_ox=1.0, C_red=1.0)
    assert i0_2 == pytest.approx(2.0 * i0_1)


def test_fit_tafel_plot_recovers_known_parameters():
    i0_true, alpha, n = 3e-6, 0.5, 1
    eta = np.linspace(0.2, 0.4, 15)
    i = butler_volmer_current_density(i0_true, eta, alpha=alpha, n=n)
    fit = fit_tafel_plot(eta, i)
    assert fit.exchange_current_density == pytest.approx(i0_true, rel=1e-3)
    assert fit.tafel_slope == pytest.approx(tafel_slope(alpha=alpha, n=n), rel=1e-3)
    assert fit.r_squared == pytest.approx(1.0, abs=1e-6)
