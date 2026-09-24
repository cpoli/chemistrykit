"""Tests for chemistrykit.statmech.systems.virial against closed-form second virial coefficients."""

import numpy as np
import pytest
from scipy.optimize import brentq

from chemistrykit.constants import K_B, NA
from chemistrykit.md.systems.lj_fluid import LennardJones
from chemistrykit.statmech.systems.virial import second_virial_coefficient

SIGMA = 3.4e-10
B0 = 2.0 * np.pi / 3.0 * NA * SIGMA**3


def test_hard_spheres():
    hs = lambda r: np.inf if r < SIGMA else 0.0
    assert second_virial_coefficient(hs, 100.0, SIGMA, [SIGMA]) == pytest.approx(B0, rel=1e-10)


@pytest.mark.parametrize("T", [80.0, 200.0, 1000.0])
def test_square_well_closed_form(T):
    eps, lam = 100.0 * K_B, 1.5
    sw = lambda r: np.inf if r < SIGMA else (-eps if r < lam * SIGMA else 0.0)
    expected = B0 * (1.0 - (lam**3 - 1.0) * np.expm1(eps / (K_B * T)))
    assert second_virial_coefficient(sw, T, SIGMA, [SIGMA, lam * SIGMA]) == pytest.approx(expected, rel=1e-9)


def test_lennard_jones_boyle_temperature():
    eps = 120.0 * K_B
    lj = LennardJones(eps, SIGMA)
    T_boyle = brentq(lambda T: second_virial_coefficient(lj.energy, T, SIGMA), 100.0, 1000.0)
    assert T_boyle / 120.0 == pytest.approx(3.4179, abs=1e-3)
