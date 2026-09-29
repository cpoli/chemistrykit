"""Tests for chemistrykit.md.utils.xyz_io."""

import numpy as np
import pytest

from chemistrykit.md.systems.lj_fluid import LJFluid
from chemistrykit.md.utils.xyz_io import read_xyz, write_xyz


def test_round_trip_with_scale_and_symbols(tmp_path):
    frames = np.random.default_rng(1).random((2, 3, 3))
    path = tmp_path / "t.xyz"
    write_xyz(path, frames, symbols=["O", "H", "H"], scale=3.4, comment="water")
    traj = read_xyz(path, scale=3.4)
    assert traj.symbols == ["O", "H", "H"]
    np.testing.assert_allclose(traj.positions, frames, atol=1e-9)
    assert all("water" in c for c in traj.comments)


def test_2d_positions_are_padded_with_zero_z(tmp_path):
    path = tmp_path / "t.xyz"
    write_xyz(path, np.ones((4, 2)))
    traj = read_xyz(path)
    assert traj.positions.shape == (1, 4, 3)
    np.testing.assert_allclose(traj.positions[0, :, 2], 0.0)


def test_mdresult_writes_lattice_and_time(tmp_path):
    fluid = LJFluid.from_lattice(n_per_side=6, density=0.6, temperature=1.0)
    result = fluid.run(dt=0.005, n_steps=4, sample_every=2)
    path = tmp_path / "lj.xyz"
    write_xyz(path, result, scale=3.405)
    traj = read_xyz(path, scale=3.405)
    assert traj.positions.shape == result.positions.shape
    np.testing.assert_allclose(traj.positions, result.positions, atol=1e-9)
    assert traj.comments[0].startswith('Lattice="')
    assert "Time=" in traj.comments[-1]


def test_comment_is_a_key_value_pair_in_extended_xyz(tmp_path):
    path = tmp_path / "c.xyz"
    write_xyz(path, np.zeros((1, 3)), box_length=2.0, comment='say "hi"')
    assert "Comment=\"say 'hi'\"" in read_xyz(path).comments[0]


def test_errors(tmp_path):
    with pytest.raises(ValueError):
        write_xyz(tmp_path / "a.xyz", np.zeros((3, 3)), symbols=["H", "H"])
    bad = tmp_path / "bad.xyz"
    bad.write_text("3\ncomment\nH 0 0 0\n")
    with pytest.raises(ValueError):
        read_xyz(bad)
