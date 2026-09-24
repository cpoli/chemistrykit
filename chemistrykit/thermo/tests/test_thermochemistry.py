"""Tests for chemistrykit.thermo.systems.thermochemistry (Hess's law)."""

import numpy as np
import pytest

from chemistrykit.thermo.systems.thermochemistry import hess_law_enthalpy, reaction_enthalpy_from_formation


def test_two_step_carbon_combustion_equals_direct_route():
    assert hess_law_enthalpy([1, 1], [-110.5, -283.0]) == pytest.approx(-393.5)


def test_reversing_a_step_flips_its_sign():
    assert hess_law_enthalpy([-1], [57.2]) == pytest.approx(-57.2)


def test_path_independence_for_random_cycles():
    """Any combination of formation reactions reproduces the formation-enthalpy route."""
    rng = np.random.default_rng(0)
    dHf = rng.normal(size=5)
    nu = rng.integers(-3, 4, size=5).astype(float)
    assert hess_law_enthalpy(nu, dHf) == pytest.approx(reaction_enthalpy_from_formation(nu, dHf))


def test_methane_combustion_from_formation_enthalpies():
    assert reaction_enthalpy_from_formation([-1, -2, 1, 2], [-74.8, 0.0, -393.5, -285.8]) == pytest.approx(-890.3)


def test_shape_mismatch_raises():
    with pytest.raises(ValueError):
        hess_law_enthalpy([1, 2], [1.0])
    with pytest.raises(ValueError):
        reaction_enthalpy_from_formation([1, 2], [1.0])
