"""Tests for chemistrykit.statmech.systems.lattice_gas against closed-form results."""

import numpy as np
import pytest

from chemistrykit.constants import K_B
from chemistrykit.statmech.systems.lattice_gas import LatticeGasAdsorption

NITROGEN_MASS = 28.0 * 1.66053906660e-27


def test_coverage_is_half_at_p_half_by_construction():
    model = LatticeGasAdsorption(adsorption_energy=3.0e-20, mass=NITROGEN_MASS, T=300.0)
    assert float(model.coverage(model.p_half())) == pytest.approx(0.5)


def test_coverage_limits_at_vacuum_and_high_pressure():
    model = LatticeGasAdsorption(adsorption_energy=3.0e-20, mass=NITROGEN_MASS, T=300.0)
    assert float(model.coverage(0.0)) == pytest.approx(0.0)
    assert float(model.coverage(1.0e15)) > 0.9999


def test_coverage_monotonically_increasing_with_pressure():
    model = LatticeGasAdsorption(adsorption_energy=3.0e-20, mass=NITROGEN_MASS, T=300.0)
    P = np.logspace(0, 10, 50)
    theta = model.coverage(P)
    assert np.all(np.diff(theta) > 0)


def test_stronger_binding_lowers_p_half():
    weak = LatticeGasAdsorption(adsorption_energy=1.0e-20, mass=NITROGEN_MASS, T=300.0)
    strong = LatticeGasAdsorption(adsorption_energy=5.0e-20, mass=NITROGEN_MASS, T=300.0)
    assert strong.p_half() < weak.p_half()


def test_canonical_entropy_vanishes_at_the_extremes():
    assert LatticeGasAdsorption.canonical_entropy(0, 100) == pytest.approx(0.0)
    assert LatticeGasAdsorption.canonical_entropy(100, 100) == pytest.approx(0.0)


def test_canonical_entropy_maximal_at_half_filling():
    M = 200
    entropies = [LatticeGasAdsorption.canonical_entropy(N, M) for N in range(0, M + 1, 10)]
    assert np.argmax(entropies) == np.argmin(np.abs(np.arange(0, M + 1, 10) - M / 2))


def test_canonical_entropy_approaches_stirling_large_m_limit():
    M = 100_000
    S_half = LatticeGasAdsorption.canonical_entropy(M // 2, M)
    S_stirling = M * K_B * np.log(2.0)
    assert S_half == pytest.approx(S_stirling, rel=1e-4)


def test_rejects_nonpositive_mass_or_temperature():
    with pytest.raises(ValueError):
        LatticeGasAdsorption(adsorption_energy=1.0, mass=0.0, T=300.0)
    with pytest.raises(ValueError):
        LatticeGasAdsorption(adsorption_energy=1.0, mass=1.0, T=-1.0)
