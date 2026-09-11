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

__all__ = ["Element", "PERIODIC_TABLE", "get_element", "molar_mass"]


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
