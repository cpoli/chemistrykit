"""Tests for chemistrykit.polymer.systems.chain_statistics against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.polymer.systems.chain_statistics import (
    FLORY_EXPONENT_GOOD_SOLVENT_RENORMALIZATION_GROUP,
    FLORY_EXPONENTS,
    IdealChain,
    RealChain,
    flory_exponent,
)


def test_ideal_chain_mean_square_end_to_end_is_exactly_n_times_b_squared():
    chain = IdealChain()
    n, b = 500, 0.7
    assert chain.mean_square_end_to_end(n, b) == pytest.approx(n * b**2)


def test_ideal_chain_scales_linearly_with_n():
    chain = IdealChain()
    b = 0.5
    assert chain.mean_square_end_to_end(400, b) / chain.mean_square_end_to_end(200, b) == pytest.approx(2.0)


def test_ideal_chain_radius_of_gyration_ratio_is_exactly_six():
    chain = IdealChain()
    n, b = 300, 0.4
    ratio = chain.mean_square_end_to_end(n, b) / chain.mean_square_radius_of_gyration(n, b)
    assert ratio == pytest.approx(6.0)


def test_flory_exponent_values():
    assert flory_exponent("theta") == pytest.approx(0.5)
    assert flory_exponent("good") == pytest.approx(3.0 / 5.0)
    assert flory_exponent("poor") == pytest.approx(1.0 / 3.0)
    assert set(FLORY_EXPONENTS) == {"theta", "good", "poor"}


def test_flory_exponent_rejects_unknown_solvent():
    with pytest.raises(ValueError):
        flory_exponent("unknown")


def test_real_chain_theta_solvent_matches_ideal_chain_end_to_end():
    ideal = IdealChain()
    theta = RealChain.theta_solvent()
    n, b = 200, 0.6
    assert theta.end_to_end_distance(n, b) == pytest.approx(ideal.end_to_end_distance(n, b))


def test_real_chain_good_solvent_is_more_swollen_than_ideal():
    ideal = IdealChain()
    good = RealChain.good_solvent()
    n, b = 1000, 0.5
    assert good.end_to_end_distance(n, b) > ideal.end_to_end_distance(n, b)


def test_real_chain_poor_solvent_is_more_collapsed_than_ideal():
    ideal = IdealChain()
    poor = RealChain.poor_solvent()
    n, b = 1000, 0.5
    assert poor.end_to_end_distance(n, b) < ideal.end_to_end_distance(n, b)


def test_real_chain_scaling_exponent_matches_flory_exponent():
    good = RealChain.good_solvent()
    n = np.array([100.0, 10000.0])
    b = 0.5
    R = good.end_to_end_distance(n, b)
    # R ~ n^nu => log(R2/R1)/log(n2/n1) == nu
    measured_nu = np.log(R[1] / R[0]) / np.log(n[1] / n[0])
    assert measured_nu == pytest.approx(good.nu, rel=1e-8)


def test_real_chain_rejects_out_of_range_nu():
    with pytest.raises(ValueError):
        RealChain(nu=0.0)
    with pytest.raises(ValueError):
        RealChain(nu=1.0)


def test_good_solvent_renormalization_group_matches_documented_constant():
    de_gennes = RealChain.good_solvent_renormalization_group()
    assert de_gennes.nu == pytest.approx(FLORY_EXPONENT_GOOD_SOLVENT_RENORMALIZATION_GROUP)


def test_good_solvent_renormalization_group_differs_from_flory_by_less_than_two_percent():
    flory = RealChain.good_solvent()
    de_gennes = RealChain.good_solvent_renormalization_group()
    relative_difference = abs(flory.nu - de_gennes.nu) / flory.nu
    assert relative_difference == pytest.approx(0.02, abs=1e-9)
    # both are still "good solvent" -- more swollen than the ideal chain
    ideal = IdealChain()
    n, b = 1000, 0.5
    assert de_gennes.end_to_end_distance(n, b) > ideal.end_to_end_distance(n, b)
