"""Tests for chemistrykit.surface.systems.gibbs_adsorption against closed-form results."""

import numpy as np
import pytest

from chemistrykit.constants import R
from chemistrykit.surface.systems.gibbs_adsorption import gibbs_surface_excess, szyszkowski_surface_tension
from chemistrykit.surface.systems.langmuir import langmuir_coverage


def test_szyszkowski_pure_solvent_is_gamma0():
    assert szyszkowski_surface_tension(0.0, gamma0=0.072, Gamma_max=5e-6, K=1.0, T=298.15) == pytest.approx(0.072)


def test_gibbs_of_szyszkowski_is_langmuir():
    c = np.logspace(-3, 2, 4001)
    gamma = szyszkowski_surface_tension(c, gamma0=0.072, Gamma_max=4e-6, K=3.0, T=293.0)
    Gamma = gibbs_surface_excess(c, gamma, T=293.0)
    np.testing.assert_allclose(Gamma[1:-1], 4e-6 * langmuir_coverage(3.0, c[1:-1]), rtol=1e-5)


def test_gibbs_linear_in_ln_c_gives_constant_excess():
    c = np.logspace(-2, 1, 50)
    T, Gamma0 = 300.0, 2e-6
    gamma = 0.07 - R * T * Gamma0 * np.log(c)
    np.testing.assert_allclose(gibbs_surface_excess(c, gamma, T), Gamma0, rtol=1e-10)


def test_gibbs_rejects_nonpositive_concentration():
    with pytest.raises(ValueError):
        gibbs_surface_excess([0.0, 1.0], [0.07, 0.06], T=300.0)
