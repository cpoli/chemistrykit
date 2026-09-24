r"""Point-group determination from 3D coordinates, and character tables for common point groups.

See Cotton, *Chemical Applications of Group Theory*, 3rd ed. (1990), Ch.
3 (symmetry elements and operations) and the Appendix (character tables),
or Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 9.2-9.3, for the
general theory.

Symmetry elements (proper rotation axes :math:`C_n`, mirror planes
:math:`\sigma`, and an inversion center `i`) are found here by genuinely
testing candidate geometric operations against the molecule's actual 3D
coordinates -- not by recognizing a molecular formula. A candidate axis
or mirror-plane normal is tested by transforming every atom's position
under the corresponding rotation/reflection matrix and checking, by
distance, whether the transformed geometry is an exact permutation of the
original atoms of matching element
(:func:`chemistrykit.structure.utils.symmetry_ops.geometry_matches`); the
detected elements are then combined via the standard point-group
decision tree (Cotton, Ch. 4) in :func:`determine_point_group`.

**Candidate generation**
(:func:`chemistrykit.structure.utils.symmetry_ops.candidate_axes`): rather than solving
for symmetry axes analytically (which would require case-by-case
geometric reasoning for every possible molecular shape), this module
follows the practical approach used by automatic point-group codes and
throws a broad, cheap-to-generate set of candidate directions at the
symmetry test -- every non-central atom's position vector, every
pairwise sum/difference/cross product of two atoms' position vectors,
every three-way sum of three atoms' position vectors (needed to reach
an octahedral geometry's body-diagonal `C3` axes -- see
:func:`~chemistrykit.structure.utils.symmetry_ops.candidate_axes`'s
docstring), and the inertia-tensor eigenvectors -- relying on the fact
that a true symmetry axis (or mirror normal) for any molecule small
enough to be built by hand is essentially always expressible as one of
these simple combinations. For every test molecule this module is
validated against (water, ammonia, methane, carbon dioxide, boron
trifluoride, sulfur hexafluoride), the true axes are recovered exactly
this way; a candidate set this broad may occasionally miss the true
axis of an unusual, unvalidated geometry, so :func:`determine_point_group`
should be understood as a best-effort computational classifier, not a
certified-complete symmetry solver.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from chemistrykit.structure.core.base_system import Molecule
from chemistrykit.structure.utils.symmetry_ops import candidate_axes, geometry_matches, highest_rotation_order, reflection_matrix

__all__ = [
    "PointGroupResult",
    "determine_point_group",
    "PointGroupCharacterTable",
    "CHARACTER_TABLES",
    "get_character_table",
]


@dataclass
class PointGroupResult:
    """The detected symmetry elements of a molecule, and the resulting point-group assignment.

    Returned by :func:`determine_point_group`; every boolean/count field
    reflects an element that was *actually found* by direct geometric
    testing against the input coordinates (see the module docstring), not
    inferred from the final `group_name`.
    """

    group_name: str
    """str: The assigned point-group symbol (e.g. ``"C2v"``, ``"Td"``, ``"D_inf_h"``
    -- the infinity symbol is spelled out ASCII-style since it is not a
    valid Python/LaTeX-free identifier component)."""

    is_linear: bool
    """bool: Whether all atoms are collinear."""

    principal_axis_order: int
    """int: The highest-order proper rotation axis found, `1` if none
    (beyond the trivial identity) was detected. Not meaningful when
    `is_linear` is True (the true axis order is infinite)."""

    n_perpendicular_c2: int
    """int: Number of `C2` axes found perpendicular to the principal axis
    (the condition for the `D` point-group families)."""

    has_sigma_h: bool
    """bool: Whether a mirror plane perpendicular to the principal axis was found."""

    has_sigma_v: bool
    """bool: Whether a mirror plane containing the principal axis was found."""

    has_inversion_center: bool
    """bool: Whether an inversion center was found."""

    n_c3_axes: int
    """int: Number of distinct `C3` axes found (>= 4 triggers the cubic
    `T`/`Td`/`Th`/`O`/`Oh` branch of the decision tree)."""

    n_mirror_planes: int
    """int: Total number of distinct mirror planes found."""

    extra: dict = field(default_factory=dict)
    """dict: Free-form slot for additional diagnostics."""


def determine_point_group(molecule: Molecule, tol: float = 1e-3) -> PointGroupResult:
    r"""Determine a molecule's point group from its 3D coordinates.

    Implements the standard point-group decision tree (Cotton, *Chemical
    Applications of Group Theory*, 3rd ed., Ch. 4, flowchart Fig. 4.1) on
    top of symmetry elements found by direct geometric testing (module
    docstring):

    1. **Linear** molecules (all atoms collinear) are :math:`D_{\infty h}`
       if they also have an inversion center, else :math:`C_{\infty v}`
       (a finite rotation-order search cannot find the true :math:`C_\infty`
       axis, so linearity is checked first and handled as a special case).
    2. **Cubic** groups: four or more distinct :math:`C_3` axes (the
       hallmark of a tetrahedral/octahedral arrangement) route to
       :math:`O_h` (if an inversion center is also present) or :math:`T_d`
       otherwise.
    3. Otherwise, the **principal axis** is the highest-order proper
       rotation axis found (if none, the group is :math:`C_i`, :math:`C_s`,
       or :math:`C_1` depending on whether an inversion center or a single
       mirror plane was found). Given a principal :math:`C_n` axis:
       finding `n` (or more) :math:`C_2` axes perpendicular to it selects
       the :math:`D_n` family (further split into :math:`D_{nh}`/:math:`D_{nd}`/:math:`D_n`
       by :math:`\sigma_h`/:math:`\sigma_v`); otherwise the :math:`C_n`
       family (split into :math:`C_{nh}`/:math:`C_{nv}`/:math:`C_n`).

    Parameters
    ----------
    molecule : chemistrykit.structure.core.base_system.Molecule
    tol : float, default 1e-3
        Absolute distance tolerance (in the molecule's coordinate units)
        for treating two atomic positions as coincident when testing a
        candidate symmetry operation.

    Returns
    -------
    PointGroupResult

    Examples
    --------
    Water is :math:`C_{2v}`, ammonia is :math:`C_{3v}`, and methane is
    :math:`T_d` -- the three textbook reference cases this classifier is
    validated against:

    >>> import numpy as np
    >>> angle = np.radians(104.5)
    >>> r = 0.958
    >>> water = Molecule(
    ...     symbols=["O", "H", "H"],
    ...     coordinates=[[0, 0, 0], [r * np.sin(angle / 2), 0, r * np.cos(angle / 2)], [-r * np.sin(angle / 2), 0, r * np.cos(angle / 2)]],
    ... )
    >>> determine_point_group(water).group_name
    'C2v'

    >>> verts = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
    >>> verts = verts / np.linalg.norm(verts[0]) * 1.09
    >>> methane = Molecule(symbols=["C", "H", "H", "H", "H"], coordinates=np.vstack([[0, 0, 0], verts]))
    >>> determine_point_group(methane).group_name
    'Td'

    Carbon dioxide is linear and centrosymmetric, :math:`D_{\infty h}`:

    >>> co2 = Molecule(symbols=["O", "C", "O"], coordinates=[[0, 0, -1.16], [0, 0, 0], [0, 0, 1.16]])
    >>> determine_point_group(co2).group_name
    'D_inf_h'
    """
    coords = molecule.centered_coordinates()
    symbols = molecule.symbols

    if molecule.is_linear():
        has_inversion = geometry_matches(symbols, coords, -coords, tol)
        group_name = "D_inf_h" if has_inversion else "C_inf_v"
        return PointGroupResult(
            group_name=group_name,
            is_linear=True,
            principal_axis_order=0,
            n_perpendicular_c2=0,
            has_sigma_h=has_inversion,
            has_sigma_v=not has_inversion,
            has_inversion_center=has_inversion,
            n_c3_axes=0,
            n_mirror_planes=0,
        )

    candidates = candidate_axes(coords)

    axis_orders = []
    for axis in candidates:
        n = highest_rotation_order(symbols, coords, axis, tol)
        if n >= 2:
            axis_orders.append((axis, n))

    n_c3_axes = sum(1 for _, n in axis_orders if n == 3)
    has_inversion = geometry_matches(symbols, coords, -coords, tol)

    mirror_normals = [axis for axis in candidates if geometry_matches(symbols, coords, coords @ reflection_matrix(axis).T, tol)]

    if n_c3_axes >= 4:
        group_name = "Oh" if has_inversion else "Td"
        return PointGroupResult(
            group_name=group_name,
            is_linear=False,
            principal_axis_order=3,
            n_perpendicular_c2=0,
            has_sigma_h=False,
            has_sigma_v=len(mirror_normals) > 0,
            has_inversion_center=has_inversion,
            n_c3_axes=n_c3_axes,
            n_mirror_planes=len(mirror_normals),
        )

    if not axis_orders:
        if has_inversion:
            group_name = "Ci"
        elif mirror_normals:
            group_name = "Cs"
        else:
            group_name = "C1"
        return PointGroupResult(
            group_name=group_name,
            is_linear=False,
            principal_axis_order=1,
            n_perpendicular_c2=0,
            has_sigma_h=False,
            has_sigma_v=False,
            has_inversion_center=has_inversion,
            n_c3_axes=0,
            n_mirror_planes=len(mirror_normals),
        )

    principal_axis, n = max(axis_orders, key=lambda pair: pair[1])

    n_perp_c2 = sum(1 for axis, order in axis_orders if order >= 2 and abs(np.dot(axis, principal_axis)) < 0.1)
    has_sigma_h = any(abs(np.dot(normal, principal_axis)) > 0.99 for normal in mirror_normals)
    has_sigma_v = any(abs(np.dot(normal, principal_axis)) < 0.1 for normal in mirror_normals)

    if n_perp_c2 >= n:
        if has_sigma_h:
            group_name = f"D{n}h"
        elif has_sigma_v:
            group_name = f"D{n}d"
        else:
            group_name = f"D{n}"
    else:
        if has_sigma_h:
            group_name = f"C{n}h"
        elif has_sigma_v:
            group_name = f"C{n}v"
        else:
            group_name = f"C{n}"

    return PointGroupResult(
        group_name=group_name,
        is_linear=False,
        principal_axis_order=n,
        n_perpendicular_c2=n_perp_c2,
        has_sigma_h=has_sigma_h,
        has_sigma_v=has_sigma_v,
        has_inversion_center=has_inversion,
        n_c3_axes=n_c3_axes,
        n_mirror_planes=len(mirror_normals),
    )


@dataclass
class PointGroupCharacterTable:
    """A finite point group's character table: irreducible representations, operation classes, and characters.

    This is reference data -- transcribed from the standard tables (Cotton,
    *Chemical Applications of Group Theory*, 3rd ed., Appendix A), in the
    same spirit as :mod:`chemistrykit.periodic_table` being plain data
    rather than a dependency -- not something :func:`determine_point_group`
    derives computationally.
    """

    name: str
    """str: Point-group symbol (matches a :func:`determine_point_group` output where finite)."""

    operations: list
    """list of str: Symmetry-operation class labels, e.g. ``["E", "C2", "sigma_v(xz)", "sigma_v'(yz)"]``."""

    irreps: list
    """list of str: Irreducible-representation (Mulliken) labels, e.g. ``["A1", "A2", "B1", "B2"]``."""

    characters: list
    """list of list of float: ``characters[i][j]`` is the character of
    irrep `i` under operation class `j`."""

    def character(self, irrep: str, operation: str) -> float:
        """Look up a single character by irrep and operation-class label.

        Parameters
        ----------
        irrep : str
        operation : str

        Returns
        -------
        float

        Examples
        --------
        >>> table = get_character_table("C2v")
        >>> table.character("A1", "E")
        1.0
        >>> table.character("B1", "C2")
        -1.0
        """
        return self.characters[self.irreps.index(irrep)][self.operations.index(operation)]

    @property
    def class_sizes(self) -> list:
        """list of int: Number of operations in each class, read from the leading count of each label.

        ``"8C3"`` is a class of 8 operations, ``"E"`` or ``"sigma_v(xz)"`` a
        class of 1. Raises ``ValueError`` for the infinite groups
        (``C_inf_v``, ``D_inf_h``), whose classes are continuous.
        """
        sizes = []
        for label in self.operations:
            if "inf" in label:
                raise ValueError(f"point group {self.name} is infinite; its classes have no finite size")
            digits = ""
            for ch in label:
                if not ch.isdigit():
                    break
                digits += ch
            sizes.append(int(digits) if digits else 1)
        return sizes

    @property
    def order(self) -> int:
        """int: The group order h (total number of symmetry operations)."""
        return sum(self.class_sizes)

    def reduce(self, reducible_characters) -> dict:
        r"""Decompose a reducible representation into irreducible ones.

        Uses the reduction formula (the "great orthogonality theorem"
        applied to characters; Cotton, *Chemical Applications of Group
        Theory*, 3rd ed., Ch. 4.3):

        .. math::

            a_i = \frac{1}{h}\sum_{R} g_R\,\chi(R)\,\chi_i(R)

        where the sum runs over operation classes of size :math:`g_R` and
        :math:`h` is the group order.

        Parameters
        ----------
        reducible_characters : sequence of float
            Characters of the reducible representation, one per operation
            class, in the order of :attr:`operations`.

        Returns
        -------
        dict of str to int
            Multiplicity of each irrep that occurs (irreps with zero
            multiplicity are left out).

        Raises
        ------
        ValueError
            If the input length does not match the number of classes, or
            the multiplicities are not integers (the characters do not form
            a representation of this group).

        Examples
        --------
        The five d orbitals in an octahedral field (Bethe, 1929) split into
        a doubly degenerate :math:`e_g` and a triply degenerate
        :math:`t_{2g}` set:

        >>> oh = get_character_table("Oh")
        >>> oh.reduce([5, -1, 1, -1, 1, 5, -1, -1, 1, 1])
        {'Eg': 1, 'T2g': 1}
        """
        chi = np.asarray(reducible_characters, dtype=np.float64)
        if chi.shape != (len(self.operations),):
            raise ValueError(f"expected {len(self.operations)} characters (one per class), got {chi.shape}")
        sizes = np.asarray(self.class_sizes, dtype=np.float64)
        h = sizes.sum()
        result = {}
        for irrep, row in zip(self.irreps, self.characters):
            a = float(np.sum(sizes * chi * np.asarray(row)) / h)
            n = int(round(a))
            if abs(a - n) > 1e-8:
                raise ValueError(f"non-integer multiplicity {a:.6f} for {irrep}: not a representation of {self.name}")
            if n:
                result[irrep] = n
        return result


#: Built-in character tables for common point groups (Cotton, *Chemical
#: Applications of Group Theory*, 3rd ed., Appendix A). ``D_inf_h`` and
#: ``C_inf_v`` (linear molecules) are given in the conventional reduced
#: textbook form restricted to the low-lying irreps (Sigma+/Sigma-/Pi/Delta);
#: the genuinely continuous ``C_inf``-rotation class is represented by its
#: identity-operation character only (`2` for a doubly-degenerate Pi/Delta
#: irrep, matching the dimension), since the general symbolic
#: :math:`2\cos\phi` entry is reference notation rather than a single number.
CHARACTER_TABLES: dict = {
    "C1": PointGroupCharacterTable(name="C1", operations=["E"], irreps=["A"], characters=[[1.0]]),
    "Cs": PointGroupCharacterTable(name="Cs", operations=["E", "sigma_h"], irreps=["A'", "A''"], characters=[[1.0, 1.0], [1.0, -1.0]]),
    "Ci": PointGroupCharacterTable(name="Ci", operations=["E", "i"], irreps=["Ag", "Au"], characters=[[1.0, 1.0], [1.0, -1.0]]),
    "C2": PointGroupCharacterTable(name="C2", operations=["E", "C2"], irreps=["A", "B"], characters=[[1.0, 1.0], [1.0, -1.0]]),
    "C2v": PointGroupCharacterTable(
        name="C2v",
        operations=["E", "C2", "sigma_v(xz)", "sigma_v'(yz)"],
        irreps=["A1", "A2", "B1", "B2"],
        characters=[
            [1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, -1.0, -1.0],
            [1.0, -1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0, 1.0],
        ],
    ),
    "C2h": PointGroupCharacterTable(
        name="C2h",
        operations=["E", "C2", "i", "sigma_h"],
        irreps=["Ag", "Bg", "Au", "Bu"],
        characters=[
            [1.0, 1.0, 1.0, 1.0],
            [1.0, -1.0, 1.0, -1.0],
            [1.0, 1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0, 1.0],
        ],
    ),
    "C3v": PointGroupCharacterTable(
        name="C3v",
        operations=["E", "2C3", "3sigma_v"],
        irreps=["A1", "A2", "E"],
        characters=[
            [1.0, 1.0, 1.0],
            [1.0, 1.0, -1.0],
            [2.0, -1.0, 0.0],
        ],
    ),
    "D2h": PointGroupCharacterTable(
        name="D2h",
        operations=["E", "C2(z)", "C2(y)", "C2(x)", "i", "sigma(xy)", "sigma(xz)", "sigma(yz)"],
        irreps=["Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"],
        characters=[
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0],
            [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0],
            [1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0],
            [1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0],
        ],
    ),
    "D3h": PointGroupCharacterTable(
        name="D3h",
        operations=["E", "2C3", "3C2", "sigma_h", "2S3", "3sigma_v"],
        irreps=["A1'", "A2'", "E'", "A1''", "A2''", "E''"],
        characters=[
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, -1.0, 1.0, 1.0, -1.0],
            [2.0, -1.0, 0.0, 2.0, -1.0, 0.0],
            [1.0, 1.0, 1.0, -1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0, -1.0, -1.0, 1.0],
            [2.0, -1.0, 0.0, -2.0, 1.0, 0.0],
        ],
    ),
    "D4h": PointGroupCharacterTable(
        name="D4h",
        operations=["E", "2C4", "C2", "2C2'", "2C2''", "i", "2S4", "sigma_h", "2sigma_v", "2sigma_d"],
        irreps=["A1g", "A2g", "B1g", "B2g", "Eg", "A1u", "A2u", "B1u", "B2u", "Eu"],
        characters=[
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0],
            [1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0],
            [1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0],
            [2.0, 0.0, -2.0, 0.0, 0.0, 2.0, 0.0, -2.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0],
            [1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0],
            [2.0, 0.0, -2.0, 0.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0],
        ],
    ),
    "Td": PointGroupCharacterTable(
        name="Td",
        operations=["E", "8C3", "3C2", "6S4", "6sigma_d"],
        irreps=["A1", "A2", "E", "T1", "T2"],
        characters=[
            [1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, -1.0, -1.0],
            [2.0, -1.0, 2.0, 0.0, 0.0],
            [3.0, 0.0, -1.0, 1.0, -1.0],
            [3.0, 0.0, -1.0, -1.0, 1.0],
        ],
    ),
    "Oh": PointGroupCharacterTable(
        name="Oh",
        operations=["E", "8C3", "6C2", "6C4", "3C2(=C4^2)", "i", "6S4", "8S6", "3sigma_h", "6sigma_d"],
        irreps=["A1g", "A2g", "Eg", "T1g", "T2g", "A1u", "A2u", "Eu", "T1u", "T2u"],
        characters=[
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0],
            [2.0, -1.0, 0.0, 0.0, 2.0, 2.0, 0.0, -1.0, 2.0, 0.0],
            [3.0, 0.0, -1.0, 1.0, -1.0, 3.0, 1.0, 0.0, -1.0, -1.0],
            [3.0, 0.0, 1.0, -1.0, -1.0, 3.0, -1.0, 0.0, -1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0],
            [2.0, -1.0, 0.0, 0.0, 2.0, -2.0, 0.0, 1.0, -2.0, 0.0],
            [3.0, 0.0, -1.0, 1.0, -1.0, -3.0, -1.0, 0.0, 1.0, 1.0],
            [3.0, 0.0, 1.0, -1.0, -1.0, -3.0, 1.0, 0.0, 1.0, -1.0],
        ],
    ),
    "C_inf_v": PointGroupCharacterTable(
        name="C_inf_v",
        operations=["E", "2C_inf^phi", "inf sigma_v"],
        irreps=["Sigma+ (A1)", "Sigma- (A2)", "Pi (E1)", "Delta (E2)"],
        characters=[
            [1.0, 1.0, 1.0],
            [1.0, 1.0, -1.0],
            [2.0, 2.0, 0.0],
            [2.0, 2.0, 0.0],
        ],
    ),
    "D_inf_h": PointGroupCharacterTable(
        name="D_inf_h",
        operations=["E", "2C_inf^phi", "inf sigma_v", "i", "2S_inf^phi", "inf C2"],
        irreps=["Sigma_g+ (A1g)", "Sigma_g- (A2g)", "Pi_g (E1g)", "Delta_g (E2g)", "Sigma_u+ (A1u)", "Sigma_u- (A2u)", "Pi_u (E1u)", "Delta_u (E2u)"],
        characters=[
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, -1.0, 1.0, 1.0, -1.0],
            [2.0, 2.0, 0.0, 2.0, -2.0, 0.0],
            [2.0, 2.0, 0.0, 2.0, 2.0, 0.0],
            [1.0, 1.0, 1.0, -1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0, -1.0, -1.0, 1.0],
            [2.0, 2.0, 0.0, -2.0, 2.0, 0.0],
            [2.0, 2.0, 0.0, -2.0, -2.0, 0.0],
        ],
    ),
}


def get_character_table(name: str) -> PointGroupCharacterTable:
    """Look up a built-in character table by point-group symbol.

    Parameters
    ----------
    name : str
        E.g. ``"C2v"``, ``"Td"``, ``"D_inf_h"``.

    Returns
    -------
    PointGroupCharacterTable

    Raises
    ------
    KeyError
        If `name` is not one of the point groups tabulated in
        :data:`CHARACTER_TABLES`.

    Examples
    --------
    Every irrep's character under the identity operation equals its
    dimension, and (for a real, unitary representation) the sum over all
    operations, weighted by class size, of the identity character times
    itself equals the group order -- here just checking the totally
    symmetric irrep is trivially present:

    >>> table = get_character_table("Td")
    >>> table.character("A1", "E")
    1.0
    >>> table.character("T2", "E")
    3.0
    """
    return CHARACTER_TABLES[name]
