"""Tests for chemistrykit.md.utils.pbc and neighbor_list against closed-form/brute-force results."""

import numpy as np

from chemistrykit.md.utils.neighbor_list import VerletNeighborList, build_neighbor_list
from chemistrykit.md.utils.pbc import minimum_image_displacement, wrap_positions


def test_minimum_image_displacement_wraps_across_boundary():
    r_i = np.array([0.5, 0.0, 0.0])
    r_j = np.array([9.5, 0.0, 0.0])
    diff = minimum_image_displacement(r_i, r_j, box_length=10.0)
    np.testing.assert_allclose(diff, [1.0, 0.0, 0.0])


def test_minimum_image_displacement_no_box_is_plain_difference():
    r_i = np.array([0.5, 0.0, 0.0])
    r_j = np.array([9.5, 0.0, 0.0])
    diff = minimum_image_displacement(r_i, r_j, box_length=None)
    np.testing.assert_allclose(diff, r_i - r_j)


def test_wrap_positions_brings_negative_and_overflowing_coordinates_into_cell():
    wrapped = wrap_positions(np.array([-0.5, 10.2, 5.0]), box_length=10.0)
    np.testing.assert_allclose(wrapped, [9.5, 0.2, 5.0])
    assert np.all(wrapped >= 0.0) and np.all(wrapped < 10.0)


def test_wrap_positions_no_box_returns_unchanged():
    positions = np.array([-1.0, 20.0])
    np.testing.assert_allclose(wrap_positions(positions, None), positions)


def test_build_neighbor_list_matches_brute_force():
    rng = np.random.default_rng(0)
    positions = rng.uniform(0.0, 10.0, size=(30, 3))
    box_length = 10.0
    cutoff = 3.0
    pairs_i, pairs_j = build_neighbor_list(positions, box_length, cutoff)
    n = positions.shape[0]
    brute = set()
    for i in range(n):
        for j in range(i + 1, n):
            diff = positions[i] - positions[j]
            diff -= box_length * np.round(diff / box_length)
            if np.linalg.norm(diff) < cutoff:
                brute.add((i, j))
    found = set(zip(pairs_i.tolist(), pairs_j.tolist(), strict=True))
    assert found == brute


def test_build_neighbor_list_all_pairs_i_less_than_j():
    rng = np.random.default_rng(1)
    positions = rng.uniform(0.0, 5.0, size=(10, 3))
    pairs_i, pairs_j = build_neighbor_list(positions, box_length=5.0, cutoff=2.0)
    assert np.all(pairs_i < pairs_j)


def test_verlet_neighbor_list_rebuild_cadence():
    nlist = VerletNeighborList(cutoff=1.0, skin=0.3, rebuild_every=5)
    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    first = nlist.pairs(positions, box_length=None)
    for _ in range(4):
        again = nlist.pairs(positions, box_length=None)
        assert again[0] is first[0]  # cached, same array objects
    nlist.pairs(positions, box_length=None)  # 6th call triggers a rebuild
    rebuilt = nlist.pairs(positions, box_length=None)
    np.testing.assert_array_equal(rebuilt[0], first[0])


def test_verlet_neighbor_list_reset_forces_rebuild():
    nlist = VerletNeighborList(cutoff=1.0, rebuild_every=1000)
    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    first = nlist.pairs(positions, box_length=None)
    nlist.reset()
    second = nlist.pairs(positions, box_length=None)
    assert second[0] is not first[0]
