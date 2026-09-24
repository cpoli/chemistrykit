"""Tests for the polymer models added for the history page, against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.polymer import (
    IdealChain,
    azeotropic_feed_composition,
    bernoullian_triad_fractions,
    branching_number_average_DP,
    branching_weight_average_DP,
    carothers_gel_point,
    debye_Kc_over_R,
    fit_mark_houwink,
    flory_huggins_critical_point,
    flory_huggins_free_energy,
    flory_huggins_spinodal_chi,
    flory_schulz_weight_average_DP,
    flory_stockmayer_gel_point,
    freely_jointed_chain,
    mark_houwink_exponent_from_flory,
    mark_houwink_intrinsic_viscosity,
    mayo_lewis_copolymer_composition,
    mean_isotactic_run_length,
    osmotic_pressure_dilute_mixture,
    poisson_number_average_DP,
    poisson_number_fraction,
    poisson_pdi,
    poisson_weight_average_DP,
    rayleigh_ratio_dilute_mixture,
    sample_dyad_sequence,
    simulate_living_polymerization,
    staudinger_specific_viscosity,
    worm_like_chain_mean_square_end_to_end,
)
from chemistrykit.polymer.utils.moments import number_average, weight_average


def test_freely_jointed_chain_mean_square_end_to_end_matches_ideal_chain():
    n, b = 200, 0.7
    X = freely_jointed_chain(n, b, n_chains=20000, rng=42)
    R2 = np.sum(X[:, -1] ** 2, axis=-1).mean()
    assert R2 == pytest.approx(IdealChain().mean_square_end_to_end(n, b), rel=0.03)
    assert np.allclose(np.linalg.norm(np.diff(X, axis=1), axis=-1), b)


def test_worm_like_chain_limits_and_monotonicity():
    P = 50.0
    assert worm_like_chain_mean_square_end_to_end(1e-4, P) == pytest.approx(1e-8, rel=1e-5)
    assert worm_like_chain_mean_square_end_to_end(1e8, P) == pytest.approx(2 * P * 1e8, rel=1e-5)
    L = np.logspace(-1, 4, 50)
    assert np.all(worm_like_chain_mean_square_end_to_end(L, P) <= L**2)


def test_flory_huggins_spinodal_is_zero_curvature():
    N, h = 100.0, 1e-4
    phi = np.linspace(0.02, 0.9, 20)
    chi_s = flory_huggins_spinodal_chi(phi, N)
    curv = (flory_huggins_free_energy(phi + h, N, chi_s) - 2 * flory_huggins_free_energy(phi, N, chi_s) + flory_huggins_free_energy(phi - h, N, chi_s)) / h**2
    assert np.allclose(curv, 0.0, atol=1e-4)
    phi_c, chi_c = flory_huggins_critical_point(N)
    grid = np.linspace(0.01, 0.99, 20001)
    assert grid[np.argmin(flory_huggins_spinodal_chi(grid, N))] == pytest.approx(phi_c, abs=1e-3)
    assert flory_huggins_spinodal_chi(phi_c, N) == pytest.approx(chi_c)


def test_mark_houwink_fit_and_staudinger_limit():
    M = np.logspace(4, 6, 5)
    K, a = fit_mark_houwink(M, mark_houwink_intrinsic_viscosity(M, 0.03, 0.72))
    assert (K, a) == pytest.approx((0.03, 0.72))
    assert mark_houwink_intrinsic_viscosity(M, 2.0, 1.0) == pytest.approx(staudinger_specific_viscosity(1.0, M, 2.0))
    assert mark_houwink_exponent_from_flory(0.5) == pytest.approx(0.5)
    assert mark_houwink_exponent_from_flory(0.6) == pytest.approx(0.8)


def test_gelation_limits():
    assert flory_stockmayer_gel_point(3) == pytest.approx(0.5)
    assert carothers_gel_point(3) == pytest.approx(2 / 3)
    assert flory_stockmayer_gel_point(4) < carothers_gel_point(4)
    assert branching_weight_average_DP(0.8, 2) == pytest.approx(flory_schulz_weight_average_DP(0.8))
    assert branching_number_average_DP(0.5, 3) == pytest.approx(4.0)
    assert branching_weight_average_DP(0.5 - 1e-9, 3) > 1e8


def test_mayo_lewis_limits_and_azeotrope():
    f1 = np.linspace(0.01, 0.99, 50)
    assert mayo_lewis_copolymer_composition(f1, 1.0, 1.0) == pytest.approx(f1)
    assert mayo_lewis_copolymer_composition(f1, 0.0, 0.0) == pytest.approx(0.5)
    fstar = azeotropic_feed_composition(0.4, 0.2)
    assert mayo_lewis_copolymer_composition(fstar, 0.4, 0.2) == pytest.approx(fstar)
    # r1*r2 = 1 (ideal): F1 = r1 f1 / (r1 f1 + f2)
    assert mayo_lewis_copolymer_composition(0.3, 2.0, 0.5) == pytest.approx(2 * 0.3 / (2 * 0.3 + 0.7))


def test_poisson_distribution_moments_and_simulation():
    nu = 40.0
    x = np.arange(1, 400)
    N = poisson_number_fraction(x, nu)
    assert N.sum() == pytest.approx(1.0)
    assert number_average(x, N) == pytest.approx(poisson_number_average_DP(nu))
    assert weight_average(x, N) == pytest.approx(poisson_weight_average_DP(nu))
    assert poisson_weight_average_DP(nu) / poisson_number_average_DP(nu) == pytest.approx(poisson_pdi(nu))
    sim = simulate_living_polymerization(5000, 200000, rng=3)
    assert sim.mean() == pytest.approx(1 + nu)
    assert np.mean(sim**2) / sim.mean() ** 2 == pytest.approx(poisson_pdi(nu), rel=2e-3)


def test_tacticity_statistics():
    mm, mr, rr = bernoullian_triad_fractions(np.array([0.2, 0.5, 0.95]))
    assert mm + mr + rr == pytest.approx(1.0)
    assert mean_isotactic_run_length(0.9) == pytest.approx(10.0)
    s = sample_dyad_sequence(200000, 0.8, rng=5)
    triad_mm = np.mean(s[:-1] & s[1:])
    assert triad_mm == pytest.approx(0.64, abs=0.005)


def test_light_scattering_gives_Mw_and_osmometry_gives_Mn():
    c = np.array([2.0, 1.0, 1.0])
    M = np.array([10.0, 50.0, 200.0])
    moles = c / M
    Mw_true = weight_average(M, moles)
    Mn_true = number_average(M, moles)
    assert rayleigh_ratio_dilute_mixture(c, M, K=3.0) / (3.0 * c.sum()) == pytest.approx(Mw_true)
    T = 298.15
    assert 8.314462618 * T * c.sum() / osmotic_pressure_dilute_mixture(c, M, T) == pytest.approx(Mn_true)
    assert debye_Kc_over_R(0.0, Mw_true, A2=1.0) == pytest.approx(1 / Mw_true)
