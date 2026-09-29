"""Tests for the Semenov chain-branching explosion model."""

import numpy as np
import pytest

from chemistrykit.constants import R
from chemistrykit.kinetics.systems.explosions import ChainBranchingExplosion


def test_limits_are_roots_of_branching_factor():
    model = ChainBranchingExplosion(k_branch=1.0e3, k_wall=50.0, k_termination=1.0e3)
    P1, P2 = model.explosion_limits(800.0)
    assert P1 < P2
    assert model.branching_factor(P1, 800.0) == pytest.approx(0.0, abs=1e-9)
    assert model.branching_factor(P2, 800.0) == pytest.approx(0.0, abs=1e-9)
    assert not model.explodes(0.5 * P1, 800.0)
    assert model.explodes(np.sqrt(P1 * P2), 800.0)
    assert not model.explodes(2 * P2, 800.0)


def test_second_limit_tends_to_2kb_over_kt_without_wall_loss():
    model = ChainBranchingExplosion(k_branch=1.0e3, k_wall=1e-12, k_termination=4.0e3)
    _, P2 = model.explosion_limits(700.0)
    assert P2 / (R * 700.0) == pytest.approx(2 * 1.0e3 / 4.0e3, rel=1e-9)


def test_arrhenius_branching_gives_a_peninsula():
    model = ChainBranchingExplosion(k_branch=lambda T: 1e8 * np.exp(-70000.0 / (R * T)), k_wall=30.0, k_termination=2.0e2)
    assert model.explosion_limits(500.0) is None
    P1_lo, P2_lo = model.explosion_limits(750.0)
    P1_hi, P2_hi = model.explosion_limits(850.0)
    assert P1_hi < P1_lo and P2_hi > P2_lo


def test_carrier_growth_matches_branching_factor():
    model = ChainBranchingExplosion(k_branch=1.0e3, k_wall=50.0, k_termination=1.0e3)
    P1, P2 = model.explosion_limits(800.0)
    t = np.array([0.0, 1e-3])
    grow = model.carrier_concentration(t, np.sqrt(P1 * P2), 800.0, w0=0.0, n0=1.0)
    assert grow[1] == pytest.approx(np.exp(model.branching_factor(np.sqrt(P1 * P2), 800.0) * 1e-3))
    steady = model.carrier_concentration(1e3, 0.5 * P1, 800.0, w0=2.0)
    assert steady == pytest.approx(-2.0 / model.branching_factor(0.5 * P1, 800.0))
