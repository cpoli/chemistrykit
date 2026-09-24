"""Tests for chemistrykit.quantum.systems.helium (effective-nuclear-charge variational helium)."""

import numpy as np
import pytest
import scipy.constants as sc
from scipy.optimize import minimize_scalar

from chemistrykit.quantum.systems.helium import (
    HARTREE_ENERGY,
    helium_like_variational_energy,
    optimize_helium_like_effective_charge,
)


def test_hartree_matches_codata():
    assert HARTREE_ENERGY == pytest.approx(sc.physical_constants["Hartree energy"][0], rel=1e-8)


def test_helium_closed_form_optimum():
    result = optimize_helium_like_effective_charge(2)
    assert result.effective_charge == pytest.approx(27.0 / 16.0)
    assert result.energy / HARTREE_ENERGY == pytest.approx(-((27.0 / 16.0) ** 2))
    assert result.screening_constant == pytest.approx(5.0 / 16.0)
    assert result.first_order_energy / HARTREE_ENERGY == pytest.approx(-2.75)
    assert result.independent_electron_energy / HARTREE_ENERGY == pytest.approx(-4.0)


@pytest.mark.parametrize("Z", [1, 2, 3, 4, 5])
def test_numerical_minimum_matches_closed_form(Z):
    res = minimize_scalar(lambda z: helium_like_variational_energy(z, Z=Z) / HARTREE_ENERGY, bounds=(0.1, 10.0), method="bounded", options={"xatol": 1e-10})
    assert res.x == pytest.approx(Z - 5.0 / 16.0, abs=1e-6)


def test_variational_bound_above_exact_helium_energy():
    exact = -2.903724  # nonrelativistic helium ground state, hartree
    result = optimize_helium_like_effective_charge(2)
    energy = result.energy / HARTREE_ENERGY
    assert exact < energy < result.first_order_energy / HARTREE_ENERGY
    assert abs(energy - exact) / abs(exact) < 0.02


def test_hydride_predicted_unbound():
    assert optimize_helium_like_effective_charge(1).ionization_energy < 0
    assert optimize_helium_like_effective_charge(2).ionization_energy > 0


def test_array_input_and_validation():
    energies = helium_like_variational_energy(np.array([1.5, 1.6875, 2.0]))
    assert energies.shape == (3,)
    assert np.argmin(energies) == 1
    with pytest.raises(ValueError):
        helium_like_variational_energy(-1.0)
    with pytest.raises(ValueError):
        optimize_helium_like_effective_charge(0.2)
