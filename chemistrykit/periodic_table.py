"""A minimal, self-contained periodic-table data table.

This is plain data (atomic number, symbol, name, standard atomic weight),
not a cheminformatics dependency -- the same spirit as
:mod:`chemistrykit.constants` holding plain numeric values rather than
depending on a units library. Standard atomic weights are the 2021 IUPAC
Commission on Isotopic Abundances and Atomic Weights (CIAAW) conventional
values, rounded to 4 significant figures (sufficient for the stoichiometry
and molar-mass calculations this package needs; an application requiring
CIAAW's full uncertainty intervals should consult the primary table
directly). Radioactive elements with no stable isotope (Tc, Pm, and
everything past Bi) are given the mass number of their longest-lived known
isotope in brackets, following the same IUPAC convention, as a plain float.

Covers Z=1 (hydrogen) through Z=54 (xenon) -- the first six periods'
main-group and first-row transition elements -- which includes every
element referenced elsewhere in chemistrykit's worked examples (e.g. Na,
Cl, Ag, I).
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "Element",
    "PERIODIC_TABLE",
    "get_element",
    "molar_mass",
    "PAULING_ELECTRONEGATIVITY",
    "MAIN_GROUP_VALENCE_ELECTRONS",
    "electronegativity",
    "valence_electrons",
]


@dataclass(frozen=True)
class Element:
    """A single periodic-table entry.

    Parameters
    ----------
    atomic_number : int
        Number of protons, Z.
    symbol : str
        One- or two-letter element symbol.
    name : str
        Full element name.
    atomic_mass : float
        Standard atomic weight, in unified atomic mass units (u) / g
        mol⁻¹ (numerically identical in either unit).
    """

    atomic_number: int
    symbol: str
    name: str
    atomic_mass: float


_ELEMENTS = [
    (1, "H", "Hydrogen", 1.008),
    (2, "He", "Helium", 4.0026),
    (3, "Li", "Lithium", 6.94),
    (4, "Be", "Beryllium", 9.0122),
    (5, "B", "Boron", 10.81),
    (6, "C", "Carbon", 12.011),
    (7, "N", "Nitrogen", 14.007),
    (8, "O", "Oxygen", 15.999),
    (9, "F", "Fluorine", 18.998),
    (10, "Ne", "Neon", 20.180),
    (11, "Na", "Sodium", 22.990),
    (12, "Mg", "Magnesium", 24.305),
    (13, "Al", "Aluminium", 26.982),
    (14, "Si", "Silicon", 28.085),
    (15, "P", "Phosphorus", 30.974),
    (16, "S", "Sulfur", 32.06),
    (17, "Cl", "Chlorine", 35.45),
    (18, "Ar", "Argon", 39.948),
    (19, "K", "Potassium", 39.098),
    (20, "Ca", "Calcium", 40.078),
    (21, "Sc", "Scandium", 44.956),
    (22, "Ti", "Titanium", 47.867),
    (23, "V", "Vanadium", 50.942),
    (24, "Cr", "Chromium", 51.996),
    (25, "Mn", "Manganese", 54.938),
    (26, "Fe", "Iron", 55.845),
    (27, "Co", "Cobalt", 58.933),
    (28, "Ni", "Nickel", 58.693),
    (29, "Cu", "Copper", 63.546),
    (30, "Zn", "Zinc", 65.38),
    (31, "Ga", "Gallium", 69.723),
    (32, "Ge", "Germanium", 72.630),
    (33, "As", "Arsenic", 74.922),
    (34, "Se", "Selenium", 78.971),
    (35, "Br", "Bromine", 79.904),
    (36, "Kr", "Krypton", 83.798),
    (37, "Rb", "Rubidium", 85.468),
    (38, "Sr", "Strontium", 87.62),
    (39, "Y", "Yttrium", 88.906),
    (40, "Zr", "Zirconium", 91.224),
    (41, "Nb", "Niobium", 92.906),
    (42, "Mo", "Molybdenum", 95.95),
    (43, "Tc", "Technetium", 98.0),  # no stable isotope; mass number of Tc-98
    (44, "Ru", "Ruthenium", 101.07),
    (45, "Rh", "Rhodium", 102.91),
    (46, "Pd", "Palladium", 106.42),
    (47, "Ag", "Silver", 107.87),
    (48, "Cd", "Cadmium", 112.41),
    (49, "In", "Indium", 114.82),
    (50, "Sn", "Tin", 118.71),
    (51, "Sb", "Antimony", 121.76),
    (52, "Te", "Tellurium", 127.60),
    (53, "I", "Iodine", 126.90),
    (54, "Xe", "Xenon", 131.29),
]

#: Mapping of element symbol -> :class:`Element`, Z=1 (H) through Z=54 (Xe).
PERIODIC_TABLE: dict[str, Element] = {symbol: Element(z, symbol, name, mass) for z, symbol, name, mass in _ELEMENTS}

#: The same entries, keyed by atomic number instead of symbol.
PERIODIC_TABLE_BY_NUMBER: dict[int, Element] = {e.atomic_number: e for e in PERIODIC_TABLE.values()}


def get_element(symbol_or_number: str | int) -> Element:
    """Look up an :class:`Element` by symbol or atomic number.

    Parameters
    ----------
    symbol_or_number : str or int
        Element symbol (e.g. ``"Na"``) or atomic number (e.g. ``11``).

    Returns
    -------
    Element

    Examples
    --------
    >>> get_element("Na").atomic_mass
    22.99
    >>> get_element(17).symbol
    'Cl'
    """
    if isinstance(symbol_or_number, int):
        return PERIODIC_TABLE_BY_NUMBER[symbol_or_number]
    return PERIODIC_TABLE[symbol_or_number]


def molar_mass(formula_counts: dict) -> float:
    """Sum atomic masses for a simple ``{symbol: count}`` formula dict.

    Only handles a flat formula (no nested groups/parentheses) -- this is
    a small stoichiometry helper, not a chemical-formula parser.

    Parameters
    ----------
    formula_counts : dict of str to float
        E.g. ``{"Na": 1, "Cl": 1}`` for NaCl, or ``{"H": 2, "O": 1}`` for
        water.

    Returns
    -------
    float
        Molar mass in g/mol.

    Examples
    --------
    >>> round(molar_mass({"Na": 1, "Cl": 1}), 2)
    58.44
    >>> round(molar_mass({"H": 2, "O": 1}), 3)
    18.015
    """
    return sum(get_element(symbol).atomic_mass * count for symbol, count in formula_counts.items())


#: Pauling-scale electronegativities, Z=1 (H) through Z=54 (Xe), as tabulated
#: by A. L. Allred, *J. Inorg. Nucl. Chem.* 17, 215 (1961) -- the standard
#: modern compilation of L. Pauling's original scale (*The Nature of the
#: Chemical Bond*, 3rd ed., 1960, Ch. 3) reproduced in most general-chemistry
#: textbooks and the CRC Handbook of Chemistry and Physics. The noble gases
#: He, Ne, and Ar have no conventionally tabulated Pauling value (they form
#: essentially no compounds on which to calibrate one) and are omitted;
#: :func:`electronegativity` raises `KeyError` for them, same as for any
#: element past Xe.
PAULING_ELECTRONEGATIVITY: dict[str, float] = {
    "H": 2.20,
    "Li": 0.98,
    "Be": 1.57,
    "B": 2.04,
    "C": 2.55,
    "N": 3.04,
    "O": 3.44,
    "F": 3.98,
    "Na": 0.93,
    "Mg": 1.31,
    "Al": 1.61,
    "Si": 1.90,
    "P": 2.19,
    "S": 2.58,
    "Cl": 3.16,
    "K": 0.82,
    "Ca": 1.00,
    "Sc": 1.36,
    "Ti": 1.54,
    "V": 1.63,
    "Cr": 1.66,
    "Mn": 1.55,
    "Fe": 1.83,
    "Co": 1.88,
    "Ni": 1.91,
    "Cu": 1.90,
    "Zn": 1.65,
    "Ga": 1.81,
    "Ge": 2.01,
    "As": 2.18,
    "Se": 2.55,
    "Br": 2.96,
    "Kr": 3.00,
    "Rb": 0.82,
    "Sr": 0.95,
    "Y": 1.22,
    "Zr": 1.33,
    "Nb": 1.6,
    "Mo": 2.16,
    "Tc": 1.9,
    "Ru": 2.2,
    "Rh": 2.28,
    "Pd": 2.20,
    "Ag": 1.93,
    "Cd": 1.69,
    "In": 1.78,
    "Sn": 1.96,
    "Sb": 2.05,
    "Te": 2.1,
    "I": 2.66,
    "Xe": 2.60,
}

#: Number of valence (outermost-shell) electrons in the *neutral, free* atom,
#: following the main-group (representative-element) octet-rule counting
#: convention of G. N. Lewis (*J. Am. Chem. Soc.* 38, 762 (1916)) and I.
#: Langmuir (*J. Am. Chem. Soc.* 41, 868 (1919)): equal to the old-style
#: "A group" number (1, 2, then 3-8 for groups 13-18). Restricted to
#: main-group elements Z=1-54 (groups 1, 2, and 13-18) -- d-block elements
#: do not have a similarly unambiguous single "valence electron count" for
#: Lewis-structure formal-charge/oxidation-state bookkeeping (which
#: (n-1)d electrons count as "valence" is itself a modeling choice) and are
#: not included; :func:`valence_electrons` raises `KeyError` for them.
#: Used by :mod:`chemistrykit.structure.systems.lewis` for formal-charge and
#: oxidation-state assignment.
MAIN_GROUP_VALENCE_ELECTRONS: dict[str, int] = {
    "H": 1,
    "He": 2,
    "Li": 1,
    "Be": 2,
    "B": 3,
    "C": 4,
    "N": 5,
    "O": 6,
    "F": 7,
    "Ne": 8,
    "Na": 1,
    "Mg": 2,
    "Al": 3,
    "Si": 4,
    "P": 5,
    "S": 6,
    "Cl": 7,
    "Ar": 8,
    "K": 1,
    "Ca": 2,
    "Ga": 3,
    "Ge": 4,
    "As": 5,
    "Se": 6,
    "Br": 7,
    "Kr": 8,
    "Rb": 1,
    "Sr": 2,
    "In": 3,
    "Sn": 4,
    "Sb": 5,
    "Te": 6,
    "I": 7,
    "Xe": 8,
}


def electronegativity(symbol: str) -> float:
    """Look up an element's Pauling-scale electronegativity.

    Parameters
    ----------
    symbol : str
        Element symbol (e.g. ``"O"``).

    Returns
    -------
    float

    Raises
    ------
    KeyError
        If `symbol` has no conventionally tabulated Pauling value (e.g. a
        noble gas, or an element past Xe).

    Examples
    --------
    >>> electronegativity("F")
    3.98
    >>> electronegativity("O") > electronegativity("H")
    True
    """
    return PAULING_ELECTRONEGATIVITY[symbol]


def valence_electrons(symbol: str) -> int:
    """Look up a main-group element's neutral-atom valence-electron count.

    Parameters
    ----------
    symbol : str
        Element symbol (e.g. ``"O"``).

    Returns
    -------
    int

    Raises
    ------
    KeyError
        If `symbol` is not a main-group element covered by
        :data:`MAIN_GROUP_VALENCE_ELECTRONS` (e.g. a d-block metal).

    Examples
    --------
    >>> valence_electrons("O")
    6
    >>> valence_electrons("Na")
    1
    """
    return MAIN_GROUP_VALENCE_ELECTRONS[symbol]
