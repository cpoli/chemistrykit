"""Tests for chemistrykit.structure.systems.point_group against known point-group assignments."""

import numpy as np
import pytest

from chemistrykit.structure.core.base_system import Molecule
from chemistrykit.structure.systems.point_group import determine_point_group, get_character_table


def _water():
    angle = np.radians(104.5)
    r = 0.958
    return Molecule(
        symbols=["O", "H", "H"],
        coordinates=[
            [0.0, 0.0, 0.0],
            [r * np.sin(angle / 2), 0.0, r * np.cos(angle / 2)],
            [-r * np.sin(angle / 2), 0.0, r * np.cos(angle / 2)],
        ],
    )


def _ammonia():
    r = 1.012
    # polar angle from the C3 axis chosen so the H-N-H angle is 106.7 degrees
    target = np.cos(np.radians(106.7))

    def f(phi):
        return np.cos(phi) ** 2 + np.sin(phi) ** 2 * np.cos(np.radians(120.0)) - target

    lo, hi = 0.01, np.pi / 2
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    phi = 0.5 * (lo + hi)
    hydrogens = [r * np.array([np.sin(phi) * np.cos(a), np.sin(phi) * np.sin(a), np.cos(phi)]) for a in (0.0, 2 * np.pi / 3, 4 * np.pi / 3)]
    return Molecule(symbols=["N", "H", "H", "H"], coordinates=np.vstack([[0.0, 0.0, 0.0]] + hydrogens))


def _methane():
    verts = np.array([[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]])
    verts = verts / np.linalg.norm(verts[0]) * 1.09
    return Molecule(symbols=["C", "H", "H", "H", "H"], coordinates=np.vstack([[0.0, 0.0, 0.0], verts]))


def _co2():
    return Molecule(symbols=["O", "C", "O"], coordinates=[[0.0, 0.0, -1.163], [0.0, 0.0, 0.0], [0.0, 0.0, 1.163]])


def _ocs():
    """Carbonyl sulfide: linear, but not centrosymmetric (asymmetric ends)."""
    return Molecule(symbols=["O", "C", "S"], coordinates=[[0.0, 0.0, -1.16], [0.0, 0.0, 0.0], [0.0, 0.0, 1.56]])


def _bf3():
    r = 1.30
    positions = [[r * np.cos(2 * np.pi * k / 3), r * np.sin(2 * np.pi * k / 3), 0.0] for k in range(3)]
    return Molecule(symbols=["B", "F", "F", "F"], coordinates=[[0.0, 0.0, 0.0]] + positions)


def test_water_is_c2v():
    result = determine_point_group(_water())
    assert result.group_name == "C2v"
    assert result.principal_axis_order == 2
    assert result.has_sigma_v is True
    assert result.has_sigma_h is False
    assert result.n_perpendicular_c2 == 0


def test_ammonia_is_c3v():
    result = determine_point_group(_ammonia())
    assert result.group_name == "C3v"
    assert result.principal_axis_order == 3


def test_methane_is_td():
    result = determine_point_group(_methane())
    assert result.group_name == "Td"
    assert result.n_c3_axes >= 4


def test_carbon_dioxide_is_d_inf_h():
    result = determine_point_group(_co2())
    assert result.group_name == "D_inf_h"
    assert result.is_linear is True
    assert result.has_inversion_center is True


def test_ocs_is_c_inf_v():
    result = determine_point_group(_ocs())
    assert result.group_name == "C_inf_v"
    assert result.is_linear is True
    assert result.has_inversion_center is False


def test_boron_trifluoride_is_d3h():
    result = determine_point_group(_bf3())
    assert result.group_name == "D3h"
    assert result.has_sigma_h is True


def test_isolated_atom_is_asymmetric_top_c1_or_higher_not_crashing():
    # A single (trivially "symmetric") atom should not crash the detector,
    # even though point-group theory doesn't meaningfully apply to it.
    atom = Molecule(symbols=["Ar"], coordinates=[[0.0, 0.0, 0.0]])
    result = determine_point_group(atom)
    assert isinstance(result.group_name, str)


def test_character_table_c2v_identity_characters_equal_dimension():
    table = get_character_table("C2v")
    for irrep in table.irreps:
        assert table.character(irrep, "E") == pytest.approx(1.0)


def test_character_table_td_e_and_t_dimensions():
    table = get_character_table("Td")
    assert table.character("A1", "E") == pytest.approx(1.0)
    assert table.character("E", "E") == pytest.approx(2.0)
    assert table.character("T1", "E") == pytest.approx(3.0)
    assert table.character("T2", "E") == pytest.approx(3.0)


def test_character_table_orthogonality_c2v():
    """The great orthogonality theorem: distinct irreps' characters are orthogonal (weighted by class size)."""
    table = get_character_table("C2v")  # all classes have size 1 in C2v
    a1 = table.characters[table.irreps.index("A1")]
    b1 = table.characters[table.irreps.index("B1")]
    dot = sum(x * y for x, y in zip(a1, b1, strict=True))
    assert dot == pytest.approx(0.0)


def test_unknown_character_table_raises():
    with pytest.raises(KeyError):
        get_character_table("not_a_real_point_group")
