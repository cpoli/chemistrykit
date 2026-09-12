"""Tests for chemistrykit.polymer.systems.step_growth against closed-form/known results."""

import pytest

from chemistrykit.polymer.systems.step_growth import (
    degree_of_polymerization,
    degree_of_polymerization_stoichiometric_imbalance,
    extent_of_reaction_for_DP,
)


def test_degree_of_polymerization_exact_carothers_relation():
    assert degree_of_polymerization(0.9) == pytest.approx(10.0)
    assert degree_of_polymerization(0.99) == pytest.approx(100.0)
    assert degree_of_polymerization(0.0) == pytest.approx(1.0)


def test_extent_of_reaction_for_DP_inverts_degree_of_polymerization():
    for p in (0.1, 0.5, 0.9, 0.99):
        Xn = degree_of_polymerization(p)
        assert extent_of_reaction_for_DP(Xn) == pytest.approx(p, rel=1e-9)


def test_degree_of_polymerization_stoichiometric_imbalance_reduces_to_carothers_at_r_equals_one():
    for p in (0.1, 0.5, 0.9, 0.99):
        assert degree_of_polymerization_stoichiometric_imbalance(p, r=1.0) == pytest.approx(degree_of_polymerization(p))


def test_degree_of_polymerization_stoichiometric_imbalance_caps_DP_at_full_conversion():
    Xn = degree_of_polymerization_stoichiometric_imbalance(p=1.0, r=0.98)
    assert Xn == pytest.approx((1.0 + 0.98) / (1.0 - 0.98))
    assert Xn < degree_of_polymerization(0.999)  # far below the unbounded Carothers value near p=1


def test_degree_of_polymerization_diverges_as_p_approaches_one():
    assert degree_of_polymerization(0.999999) > degree_of_polymerization(0.99)
