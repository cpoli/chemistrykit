"""Tests for chemistrykit.thermo.systems.equilibrium against closed-form results."""

import numpy as np
import pytest

from chemistrykit.constants import R
from chemistrykit.thermo.systems.equilibrium import (
    fit_van_t_hoff,
    kc_from_kp,
    kp_from_kc,
    reaction_quotient,
    solve_equilibrium_composition,
    van_t_hoff_equilibrium_constant,
)


def test_reaction_quotient_of_pure_products_and_reactants():
    # 2 A -> B, all activities = 1: Q = 1^1 / 1^2 = 1
    assert reaction_quotient([1.0, 1.0], [-2.0, 1.0]) == pytest.approx(1.0)


def test_reaction_quotient_matches_target_k_at_equilibrium_mole_fractions():
    """N2O4 <-> 2 NO2 at the mole fractions for Q = K = 4 (xi = sqrt(K/(4+K)))."""
    xi = np.sqrt(4.0 / 8.0)
    x_n2o4, x_no2 = (1.0 - xi) / (1.0 + xi), 2.0 * xi / (1.0 + xi)
    assert reaction_quotient([x_n2o4, x_no2], [-1.0, 2.0]) == pytest.approx(4.0)


def test_kp_kc_round_trip():
    Kc, delta_n, T = 3.7, 2.0, 400.0
    Kp = kp_from_kc(Kc, delta_n, T)
    assert kc_from_kp(Kp, delta_n, T) == pytest.approx(Kc)


def test_kp_equals_kc_when_delta_n_zero():
    assert kp_from_kc(Kc=5.0, delta_n=0.0, T=298.15) == pytest.approx(5.0)


def test_van_t_hoff_reduces_to_reference_at_reference_temperature():
    K = van_t_hoff_equilibrium_constant(T=310.0, T_ref=310.0, K_ref=2.5, delta_h=50e3)
    assert K == pytest.approx(2.5)


def test_van_t_hoff_endothermic_reaction_favored_by_heat():
    T = np.array([280.0, 300.0, 320.0])
    K = van_t_hoff_equilibrium_constant(T, T_ref=298.15, K_ref=1.0, delta_h=60e3)
    assert np.all(np.diff(K) > 0)


def test_van_t_hoff_exothermic_reaction_disfavored_by_heat():
    T = np.array([280.0, 300.0, 320.0])
    K = van_t_hoff_equilibrium_constant(T, T_ref=298.15, K_ref=1.0, delta_h=-60e3)
    assert np.all(np.diff(K) < 0)


def test_fit_van_t_hoff_recovers_exact_enthalpy():
    T = np.linspace(280.0, 360.0, 8)
    delta_h_true = 55_000.0
    K = van_t_hoff_equilibrium_constant(T, T_ref=300.0, K_ref=1.0, delta_h=delta_h_true)
    fit = fit_van_t_hoff(T, K)
    assert fit.delta_h == pytest.approx(delta_h_true, rel=1e-6)
    assert fit.r_squared == pytest.approx(1.0, abs=1e-9)


def test_solve_equilibrium_composition_matches_algebraic_extent():
    r"""N2O4 <-> 2 NO2: Q(xi) = 4*xi^2/(1-xi^2) at P = P deg, so Q=K gives xi = sqrt(K/(4+K))."""
    T = 298.15
    K_target = 4.0
    delta_g_rxn = -R * T * np.log(K_target)
    gibbs_formation = [0.0, delta_g_rxn / 2.0]
    result = solve_equilibrium_composition(
        species=("N2O4", "NO2"),
        stoich_matrix=[[-1.0], [2.0]],
        n0=[1.0, 0.0],
        gibbs_formation=gibbs_formation,
        T=T,
    )
    xi_expected = np.sqrt(K_target / (4.0 + K_target))
    assert result.converged
    assert result.extents[0] == pytest.approx(xi_expected, rel=1e-4)
    assert result.moles("N2O4") == pytest.approx(1.0 - xi_expected, rel=1e-4)
    assert result.moles("NO2") == pytest.approx(2.0 * xi_expected, rel=1e-4)


def test_solve_equilibrium_composition_conserves_atoms():
    """Total N atoms (1 per N2O4 counts as 2 N, 1 per NO2 counts as 1 N) is conserved."""
    T = 298.15
    delta_g_rxn = -R * T * np.log(2.5)
    gibbs_formation = [0.0, delta_g_rxn / 2.0]
    result = solve_equilibrium_composition(
        species=("N2O4", "NO2"),
        stoich_matrix=[[-1.0], [2.0]],
        n0=[2.0, 0.0],
        gibbs_formation=gibbs_formation,
        T=T,
    )
    total_N = 2.0 * result.moles("N2O4") + result.moles("NO2")
    assert total_N == pytest.approx(4.0, rel=1e-4)


def test_solve_equilibrium_composition_no_reaction_when_products_disfavored():
    """A strongly positive Gibbs energy of reaction should leave the system almost entirely reactant."""
    T = 298.15
    gibbs_formation = [0.0, 200_000.0]  # NO2 formation strongly disfavored
    result = solve_equilibrium_composition(
        species=("N2O4", "NO2"),
        stoich_matrix=[[-1.0], [2.0]],
        n0=[1.0, 0.0],
        gibbs_formation=gibbs_formation,
        T=T,
    )
    assert result.extents[0] == pytest.approx(0.0, abs=1e-3)
