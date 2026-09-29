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
isotope in brackets, following the same IUPAC convention, as a plain float
(Th, Pa and U, which have a characteristic terrestrial isotopic
composition, keep their CIAAW standard weights).

Covers every named element, Z=1 (hydrogen) through Z=118 (oganesson).
:func:`parse_formula` turns a condensed formula string such as
``"Ca(OH)2"`` or ``"CuSO4*5H2O"`` into the ``{symbol: count}`` dict the
rest of this module works with, and :func:`molar_mass` accepts either.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

__all__ = [
    "Element",
    "PERIODIC_TABLE",
    "get_element",
    "molar_mass",
    "parse_formula",
    "formula_charge",
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
    (55, "Cs", "Caesium", 132.91),
    (56, "Ba", "Barium", 137.33),
    (57, "La", "Lanthanum", 138.91),
    (58, "Ce", "Cerium", 140.12),
    (59, "Pr", "Praseodymium", 140.91),
    (60, "Nd", "Neodymium", 144.24),
    (61, "Pm", "Promethium", 145.0),  # no stable isotope; mass number of Pm-145
    (62, "Sm", "Samarium", 150.36),
    (63, "Eu", "Europium", 151.96),
    (64, "Gd", "Gadolinium", 157.25),
    (65, "Tb", "Terbium", 158.93),
    (66, "Dy", "Dysprosium", 162.5),
    (67, "Ho", "Holmium", 164.93),
    (68, "Er", "Erbium", 167.26),
    (69, "Tm", "Thulium", 168.93),
    (70, "Yb", "Ytterbium", 173.05),
    (71, "Lu", "Lutetium", 174.97),
    (72, "Hf", "Hafnium", 178.49),
    (73, "Ta", "Tantalum", 180.95),
    (74, "W", "Tungsten", 183.84),
    (75, "Re", "Rhenium", 186.21),
    (76, "Os", "Osmium", 190.23),
    (77, "Ir", "Iridium", 192.22),
    (78, "Pt", "Platinum", 195.08),
    (79, "Au", "Gold", 196.97),
    (80, "Hg", "Mercury", 200.59),
    (81, "Tl", "Thallium", 204.38),
    (82, "Pb", "Lead", 207.2),
    (83, "Bi", "Bismuth", 208.98),
    (84, "Po", "Polonium", 209.0),  # no stable isotope from here on except Th, Pa, U (CIAAW standard weights)
    (85, "At", "Astatine", 210.0),
    (86, "Rn", "Radon", 222.0),
    (87, "Fr", "Francium", 223.0),
    (88, "Ra", "Radium", 226.0),
    (89, "Ac", "Actinium", 227.0),
    (90, "Th", "Thorium", 232.04),
    (91, "Pa", "Protactinium", 231.04),
    (92, "U", "Uranium", 238.03),
    (93, "Np", "Neptunium", 237.0),
    (94, "Pu", "Plutonium", 244.0),
    (95, "Am", "Americium", 243.0),
    (96, "Cm", "Curium", 247.0),
    (97, "Bk", "Berkelium", 247.0),
    (98, "Cf", "Californium", 251.0),
    (99, "Es", "Einsteinium", 252.0),
    (100, "Fm", "Fermium", 257.0),
    (101, "Md", "Mendelevium", 258.0),
    (102, "No", "Nobelium", 259.0),
    (103, "Lr", "Lawrencium", 266.0),
    (104, "Rf", "Rutherfordium", 267.0),
    (105, "Db", "Dubnium", 268.0),
    (106, "Sg", "Seaborgium", 269.0),
    (107, "Bh", "Bohrium", 270.0),
    (108, "Hs", "Hassium", 269.0),
    (109, "Mt", "Meitnerium", 278.0),
    (110, "Ds", "Darmstadtium", 281.0),
    (111, "Rg", "Roentgenium", 282.0),
    (112, "Cn", "Copernicium", 285.0),
    (113, "Nh", "Nihonium", 286.0),
    (114, "Fl", "Flerovium", 289.0),
    (115, "Mc", "Moscovium", 290.0),
    (116, "Lv", "Livermorium", 293.0),
    (117, "Ts", "Tennessine", 294.0),
    (118, "Og", "Oganesson", 294.0),
]

#: Mapping of element symbol -> :class:`Element`, Z=1 (H) through Z=118 (Og).
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


_TOKEN = re.compile(r"([A-Z][a-z]?)|([(\[{])|([)\]}])|(\d+(?:\.\d+)?)|(\S)")
_CLOSING = {"(": ")", "[": "]", "{": "}"}
_CHARGE = re.compile(r"(?:\^(\d*)([+-])|(\d?)([+-]))$")


def _split_charge(formula: str) -> tuple[str, int]:
    formula = formula.strip()
    match = _CHARGE.search(formula)
    if match is None:
        return formula, 0
    if match.group(2):
        digits, sign, body = match.group(1), match.group(2), formula[: match.start()]
    elif match.group(3) and re.fullmatch(r"[A-Z][a-z]?", formula[: match.start()].strip()):
        digits, sign, body = match.group(3), match.group(4), formula[: match.start()]
    else:
        digits, sign, body = "", match.group(4), formula[: match.end() - 1]
    magnitude = int(digits) if digits else 1
    return body.strip(), magnitude if sign == "+" else -magnitude


def _parse_group(tokens: list, pos: int, closing: str | None) -> tuple[dict, int]:
    counts: dict = {}
    while pos < len(tokens):
        symbol, opening, close, number, junk = tokens[pos]
        if junk:
            raise ValueError(f"unexpected character {junk!r} in formula")
        if close:
            if close != closing:
                raise ValueError(f"unbalanced {close!r} in formula")
            return counts, pos + 1
        if number:
            raise ValueError(f"count {number} does not follow an element or group")
        if symbol:
            if symbol not in PERIODIC_TABLE:
                raise ValueError(f"unknown element symbol {symbol!r}")
            inner, pos = {symbol: 1}, pos + 1
        else:
            inner, pos = _parse_group(tokens, pos + 1, _CLOSING[opening])
        multiplier: float = 1
        if pos < len(tokens) and tokens[pos][3]:
            multiplier = float(tokens[pos][3]) if "." in tokens[pos][3] else int(tokens[pos][3])
            pos += 1
        for key, value in inner.items():
            counts[key] = counts.get(key, 0) + value * multiplier
    if closing is not None:
        raise ValueError(f"missing {closing!r} in formula")
    return counts, pos


def parse_formula(formula: str) -> dict:
    """Parse a condensed chemical formula into a ``{symbol: count}`` dict.

    Handles nested groups in ``()``, ``[]`` or ``{}`` (``"Ca(OH)2"``,
    ``"K4[Fe(CN)6]"``), non-integer counts for non-stoichiometric solids
    (``"Fe0.95O"``), and hydrate/adduct parts joined by ``*`` or ``·``,
    each with an optional leading multiplier (``"CuSO4*5H2O"``). A
    trailing ionic charge is stripped and ignored here -- see
    :func:`formula_charge`. Charge notation follows common usage: ``^``
    always marks the charge (``"SO4^2-"``, ``"Hg2^2+"``, ``"O2^-"``);
    without it, a digit before the sign is the charge of a *monatomic*
    ion (``"Fe3+"``, ``"O2-"`` for oxide) but a count in a polyatomic one
    (``"MnO4-"``, ``"NH4+"`` carry charge ±1). The electron, ``"e-"``,
    parses to an empty dict.

    Parameters
    ----------
    formula : str

    Returns
    -------
    dict of str to int or float
        Element counts, in order of first appearance.

    Raises
    ------
    ValueError
        On an unknown element symbol, unbalanced brackets, or any other
        malformed input.

    Examples
    --------
    >>> parse_formula("Ca(OH)2")
    {'Ca': 1, 'O': 2, 'H': 2}
    >>> parse_formula("CuSO4*5H2O")
    {'Cu': 1, 'S': 1, 'O': 9, 'H': 10}
    >>> parse_formula("K4[Fe(CN)6]")
    {'K': 4, 'Fe': 1, 'C': 6, 'N': 6}
    """
    body, _ = _split_charge(formula)
    if body == "e":
        return {}
    counts: dict = {}
    for part in re.split(r"[*·]", body):
        part = part.strip()
        if not part:
            raise ValueError(f"empty component in formula {formula!r}")
        lead = re.match(r"\d+", part)
        multiplier = int(lead.group()) if lead else 1
        part_counts, _ = _parse_group(_TOKEN.findall(part[lead.end() if lead else 0 :]), 0, None)
        if not part_counts:
            raise ValueError(f"no elements in formula {formula!r}")
        for key, value in part_counts.items():
            counts[key] = counts.get(key, 0) + value * multiplier
    return counts


def formula_charge(formula: str) -> int:
    """Return the ionic charge written at the end of a formula (0 if none).

    Uses the same charge notation as :func:`parse_formula`.

    Parameters
    ----------
    formula : str

    Returns
    -------
    int

    Examples
    --------
    >>> formula_charge("SO4^2-"), formula_charge("Fe3+"), formula_charge("MnO4-"), formula_charge("e-"), formula_charge("H2O")
    (-2, 3, -1, -1, 0)
    """
    return _split_charge(formula)[1]


def molar_mass(formula) -> float:
    """Molar mass of a formula string or ``{symbol: count}`` dict.

    A string is parsed with :func:`parse_formula`; any ionic charge is
    ignored (the electron-mass correction is below the 4-significant-figure
    precision of the tabulated atomic weights).

    Parameters
    ----------
    formula : str or dict of str to float
        E.g. ``"Ca(OH)2"``, or ``{"Na": 1, "Cl": 1}`` for NaCl.

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
    >>> round(molar_mass("Ca(OH)2"), 2)
    74.09
    """
    counts = parse_formula(formula) if isinstance(formula, str) else formula
    return sum(get_element(symbol).atomic_mass * count for symbol, count in counts.items())


#: Pauling-scale electronegativities, Z=1 (H) through Z=94 (Pu), as tabulated
#: by A. L. Allred, *J. Inorg. Nucl. Chem.* 17, 215 (1961) -- the standard
#: modern compilation of L. Pauling's original scale (*The Nature of the
#: Chemical Bond*, 3rd ed., 1960, Ch. 3) reproduced in most general-chemistry
#: textbooks and the CRC Handbook of Chemistry and Physics. The noble gases
#: He, Ne, and Ar have no conventionally tabulated Pauling value (they form
#: essentially no compounds on which to calibrate one) and are omitted, as
#: are Pm, Eu, Tb, Yb, Rn and everything past Pu, for which Allred gives no
#: value; :func:`electronegativity` raises `KeyError` for all of them.
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
    "Cs": 0.79,
    "Ba": 0.89,
    "La": 1.1,
    "Ce": 1.12,
    "Pr": 1.13,
    "Nd": 1.14,
    "Sm": 1.17,
    "Gd": 1.2,
    "Dy": 1.22,
    "Ho": 1.23,
    "Er": 1.24,
    "Tm": 1.25,
    "Lu": 1.27,
    "Hf": 1.3,
    "Ta": 1.5,
    "W": 2.36,
    "Re": 1.9,
    "Os": 2.2,
    "Ir": 2.2,
    "Pt": 2.28,
    "Au": 2.54,
    "Hg": 2.0,
    "Tl": 1.62,
    "Pb": 2.33,
    "Bi": 2.02,
    "Po": 2.0,
    "At": 2.2,
    "Fr": 0.7,
    "Ra": 0.9,
    "Ac": 1.1,
    "Th": 1.3,
    "Pa": 1.5,
    "U": 1.38,
    "Np": 1.36,
    "Pu": 1.28,
}

#: Number of valence (outermost-shell) electrons in the *neutral, free* atom,
#: following the main-group (representative-element) octet-rule counting
#: convention of G. N. Lewis (*J. Am. Chem. Soc.* 38, 762 (1916)) and I.
#: Langmuir (*J. Am. Chem. Soc.* 41, 868 (1919)): equal to the old-style
#: "A group" number (1, 2, then 3-8 for groups 13-18). Restricted to
#: main-group elements Z=1-88 (groups 1, 2, and 13-18) -- d- and f-block elements
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
    "Cs": 1,
    "Ba": 2,
    "Tl": 3,
    "Pb": 4,
    "Bi": 5,
    "Po": 6,
    "At": 7,
    "Rn": 8,
    "Fr": 1,
    "Ra": 2,
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
        light noble gas, or an actinide past Pu).

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
