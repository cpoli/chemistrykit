"""Tests for chemistrykit.structure.systems.vsepr against exact, closed-form VSEPR angles."""

import numpy as np
import pytest

from chemistrykit.structure.core.base_system import angle_between
from chemistrykit.structure.systems.vsepr import VSEPRGeometry, build_vsepr_molecule, domain_positions

TETRAHEDRAL_ANGLE = 109.4712206


@pytest.mark.parametrize("steric_number", [2, 3, 4, 5, 6])
def test_domain_positions_are_unit_vectors(steric_number):
    positions = domain_positions(steric_number)
    assert positions.shape == (steric_number, 3)
    norms = np.linalg.norm(positions, axis=1)
    assert norms == pytest.approx(np.ones(steric_number))


def test_domain_positions_rejects_invalid_steric_number():
    with pytest.raises(ValueError):
        domain_positions(7)


def test_linear_angle_is_180():
    verts = domain_positions(2)
    assert angle_between(verts[0], verts[1]) == pytest.approx(180.0)


def test_trigonal_planar_angles_are_120():
    verts = domain_positions(3)
    angles = [angle_between(verts[i], verts[j]) for i in range(3) for j in range(i + 1, 3)]
    assert angles == pytest.approx([120.0] * 3)


def test_tetrahedral_angles_are_exactly_109_47():
    verts = domain_positions(4)
    angles = [angle_between(verts[i], verts[j]) for i in range(4) for j in range(i + 1, 4)]
    assert angles == pytest.approx([TETRAHEDRAL_ANGLE] * 6)


def test_octahedral_has_90_and_180_degree_angles():
    verts = domain_positions(6)
    angles = sorted({round(angle_between(verts[i], verts[j]), 3) for i in range(6) for j in range(i + 1, 6)})
    assert angles == pytest.approx([90.0, 180.0])


def test_methane_shape_and_angle():
    geometry = VSEPRGeometry(steric_number=4, lone_pairs=0)
    assert geometry.shape_name == "tetrahedral"
    methane = build_vsepr_molecule(4, 0, bond_length=1.09, central_symbol="C", ligand_symbol="H")
    for i in range(1, 4):
        for j in range(i + 1, 5):
            assert methane.bond_angle(i, 0, j) == pytest.approx(TETRAHEDRAL_ANGLE)


def test_ammonia_shape_is_trigonal_pyramidal():
    geometry = VSEPRGeometry(steric_number=4, lone_pairs=1)
    assert geometry.shape_name == "trigonal pyramidal"
    assert geometry.n_bonding_domains == 3


def test_water_shape_is_bent():
    geometry = VSEPRGeometry(steric_number=4, lone_pairs=2)
    assert geometry.shape_name == "bent"
    water = build_vsepr_molecule(4, 2, bond_length=0.96, central_symbol="O", ligand_symbol="H")
    assert water.bond_angle(1, 0, 2) == pytest.approx(TETRAHEDRAL_ANGLE)


def test_sf4_seesaw_from_trigonal_bipyramidal():
    geometry = VSEPRGeometry(steric_number=5, lone_pairs=1)
    assert geometry.shape_name == "seesaw"
    assert geometry.n_bonding_domains == 4


def test_clf3_t_shaped():
    geometry = VSEPRGeometry(steric_number=5, lone_pairs=2)
    assert geometry.shape_name == "T-shaped"


def test_xef2_linear_from_five_domains():
    geometry = VSEPRGeometry(steric_number=5, lone_pairs=3)
    assert geometry.shape_name == "linear"
    molecule = build_vsepr_molecule(5, 3, bond_length=1.98, central_symbol="Xe", ligand_symbol="F")
    assert molecule.bond_angle(1, 0, 2) == pytest.approx(180.0)


def test_xef4_square_planar():
    geometry = VSEPRGeometry(steric_number=6, lone_pairs=2)
    assert geometry.shape_name == "square planar"
    molecule = build_vsepr_molecule(6, 2, bond_length=1.95, central_symbol="Xe", ligand_symbol="F")
    angles = sorted(round(molecule.bond_angle(i, 0, j), 3) for i in range(1, 5) for j in range(i + 1, 5))
    assert angles == pytest.approx([90.0, 90.0, 90.0, 90.0, 180.0, 180.0])


def test_geometry_rejects_invalid_lone_pair_count():
    with pytest.raises(ValueError):
        VSEPRGeometry(steric_number=4, lone_pairs=5)


def test_geometry_rejects_invalid_steric_number():
    with pytest.raises(ValueError):
        VSEPRGeometry(steric_number=1, lone_pairs=0)
