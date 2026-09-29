r"""Reaction stoichiometry: equation balancing, limiting reagents, yields, empirical formulas.

Balancing is the linear-algebra problem of finding a positive integer
vector :math:`\nu` in the null space of the composition matrix
:math:`A_{e,s}` (the count of element `e` in species `s`, taken positive
for reactants and negative for products, plus a row for net charge when
ions are present): :math:`A\nu = 0` states that every element, and charge,
is conserved (W. R. Smith & R. W. Missen, *Chemical Reaction Equilibrium
Analysis*, Wiley, 1982, Ch. 2). The null space is computed exactly over the rationals
with :mod:`sympy`, so the integer coefficients carry no rounding error.
An equation whose null space is not one-dimensional is either impossible
to balance or a sum of independent reactions, and is rejected.

The limiting-reagent and yield helpers follow the usual general-chemistry
definitions via the extent of reaction :math:`\xi`: species `i` with
coefficient :math:`\nu_i` and initial amount :math:`n_i` supports at most
:math:`\xi_i = n_i/\nu_i`, and the reactant with the smallest
:math:`\xi_i` is limiting (Atkins & de Paula, *Physical Chemistry*, 11th
ed., Ch. 6.1 for :math:`\xi`).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

import sympy

from chemistrykit.periodic_table import formula_charge, get_element, molar_mass, parse_formula

__all__ = [
    "BalancedEquation",
    "balance_equation",
    "reaction_extent",
    "limiting_reagent",
    "theoretical_yield",
    "percent_yield",
    "mass_percent",
    "empirical_formula",
]

_ARROW = re.compile(r"\s*(?:<=>|<->|->|=>|=|→|⇌)\s*")


@dataclass(frozen=True)
class BalancedEquation:
    """A balanced chemical equation.

    Parameters
    ----------
    reactants : dict of str to int
        Reactant formula -> stoichiometric coefficient, in input order.
    products : dict of str to int
        Product formula -> stoichiometric coefficient, in input order.
    """

    reactants: dict
    products: dict

    @property
    def coefficients(self) -> dict:
        """dict of str to int: Signed (net) coefficients, negative for reactants.

        The same sign convention as a column of
        :attr:`chemistrykit.kinetics.StoichiometricNetwork.stoich_matrix`.
        """
        signed = {species: -nu for species, nu in self.reactants.items()}
        for species, nu in self.products.items():
            signed[species] = signed.get(species, 0) + nu
        return signed

    def __str__(self) -> str:
        def side(terms: dict) -> str:
            return " + ".join(f"{nu if nu != 1 else ''}{species}" for species, nu in terms.items())

        return f"{side(self.reactants)} -> {side(self.products)}"


def _split_side(side: str) -> list[str]:
    species = [re.sub(r"^\d+\s*(?=[A-Z(\[{e])", "", term.strip()) for term in re.split(r"\s+\+\s+", side.strip())]
    if not all(species):
        raise ValueError(f"empty species in {side!r}")
    return species


def balance_equation(equation: str) -> BalancedEquation:
    r"""Balance a chemical equation, including ionic and redox equations.

    Parameters
    ----------
    equation : str
        Reactants and products separated by ``->``, ``=``, ``<->``,
        ``→`` or ``⇌``; species on each side separated by ``" + "``
        (with spaces, since ``+`` also marks a charge). Formulas use the
        notation of :func:`chemistrykit.periodic_table.parse_formula`;
        ``"e-"`` is the electron. Any coefficients already present are
        ignored and recomputed.

    Returns
    -------
    BalancedEquation
        Smallest positive integer coefficients.

    Raises
    ------
    ValueError
        If the equation cannot be balanced with all species on the sides
        given, or has more than one independent balance (it combines
        several reactions whose proportions are arbitrary).

    Examples
    --------
    >>> str(balance_equation("Fe + O2 -> Fe2O3"))
    '4Fe + 3O2 -> 2Fe2O3'
    >>> str(balance_equation("C3H8 + O2 -> CO2 + H2O"))
    'C3H8 + 5O2 -> 3CO2 + 4H2O'

    Redox equations balance charge as well as mass, with no need to
    split them into half-reactions:

    >>> str(balance_equation("MnO4- + Fe2+ + H+ -> Mn2+ + Fe3+ + H2O"))
    'MnO4- + 5Fe2+ + 8H+ -> Mn2+ + 5Fe3+ + 4H2O'
    """
    sides = _ARROW.split(equation.strip())
    if len(sides) != 2:
        raise ValueError(f"expected exactly one reaction arrow in {equation!r}")
    reactants, products = _split_side(sides[0]), _split_side(sides[1])
    species = reactants + products
    signs = [1] * len(reactants) + [-1] * len(products)
    compositions = [parse_formula(s) for s in species]
    charges = [formula_charge(s) for s in species]
    elements = list(dict.fromkeys(e for comp in compositions for e in comp))
    rows = [[sign * sympy.nsimplify(comp.get(e, 0)) for comp, sign in zip(compositions, signs, strict=True)] for e in elements]
    if any(charges):
        rows.append([sign * q for q, sign in zip(charges, signs, strict=True)])
    null = sympy.Matrix(rows).nullspace()
    if len(null) != 1:
        reason = "cannot be balanced" if not null else f"has {len(null)} independent balances"
        raise ValueError(f"{equation!r} {reason}")
    vector = null[0]
    scale = sympy.ilcm(*[sympy.fraction(v)[1] for v in vector])
    ints = [int(v * scale) for v in vector]
    divisor = math.gcd(*ints)
    ints = [v // divisor for v in ints]
    if all(v < 0 for v in ints):
        ints = [-v for v in ints]
    if not all(v > 0 for v in ints):
        raise ValueError(f"{equation!r} cannot be balanced with every species on the side given")
    n_r = len(reactants)
    return BalancedEquation(reactants=dict(zip(reactants, ints[:n_r], strict=True)), products=dict(zip(products, ints[n_r:], strict=True)))


def _as_balanced(equation) -> BalancedEquation:
    return equation if isinstance(equation, BalancedEquation) else balance_equation(equation)


def reaction_extent(equation, moles: dict) -> float:
    r"""Maximum extent of reaction :math:`\xi_{max}=\min_i n_i/\nu_i` over the reactants supplied.

    Parameters
    ----------
    equation : str or BalancedEquation
    moles : dict of str to float
        Initial amount of each reactant, in mol. A reactant left out is
        treated as being in excess.

    Returns
    -------
    float
        Extent, in mol.

    Examples
    --------
    >>> reaction_extent("H2 + O2 -> H2O", {"H2": 4.0, "O2": 1.5})
    1.5
    """
    balanced = _as_balanced(equation)
    supplied = {s: n for s, n in moles.items() if s in balanced.reactants}
    if not supplied:
        raise ValueError("moles names none of the equation's reactants")
    return min(n / balanced.reactants[s] for s, n in supplied.items())


def limiting_reagent(equation, moles: dict) -> str:
    """The reactant that runs out first.

    Parameters
    ----------
    equation : str or BalancedEquation
    moles : dict of str to float
        Initial amount of each reactant, in mol.

    Returns
    -------
    str

    Examples
    --------
    Burning 4 mol H2 in 1.5 mol O2 (2H2 + O2 -> 2H2O): 4 mol H2 would need
    2 mol O2, so O2 is limiting:

    >>> limiting_reagent("H2 + O2 -> H2O", {"H2": 4.0, "O2": 1.5})
    'O2'
    """
    balanced = _as_balanced(equation)
    supplied = {s: n for s, n in moles.items() if s in balanced.reactants}
    if not supplied:
        raise ValueError("moles names none of the equation's reactants")
    return min(supplied, key=lambda s: supplied[s] / balanced.reactants[s])


def theoretical_yield(equation, moles: dict, product: str, in_grams: bool = False) -> float:
    """Amount of `product` formed if the limiting reagent is fully consumed.

    Parameters
    ----------
    equation : str or BalancedEquation
    moles : dict of str to float
        Initial amount of each reactant, in mol.
    product : str
        A product formula, exactly as written in the equation.
    in_grams : bool, default False
        Return a mass in g (via :func:`chemistrykit.periodic_table.molar_mass`)
        instead of an amount in mol.

    Returns
    -------
    float

    Examples
    --------
    >>> theoretical_yield("H2 + O2 -> H2O", {"H2": 4.0, "O2": 1.5}, "H2O")
    3.0
    >>> round(theoretical_yield("H2 + O2 -> H2O", {"H2": 4.0, "O2": 1.5}, "H2O", in_grams=True), 2)
    54.05
    """
    balanced = _as_balanced(equation)
    if product not in balanced.products:
        raise ValueError(f"{product!r} is not a product of {balanced}")
    n = reaction_extent(balanced, moles) * balanced.products[product]
    return n * molar_mass(product) if in_grams else n


def percent_yield(actual: float, theoretical: float) -> float:
    """Percent yield, ``100 * actual / theoretical`` (both in the same units).

    Examples
    --------
    >>> percent_yield(45.0, 54.05) > 83.0
    True
    """
    if theoretical <= 0:
        raise ValueError("theoretical yield must be positive")
    return 100.0 * actual / theoretical


def mass_percent(formula) -> dict:
    """Percent composition by mass of a formula.

    Parameters
    ----------
    formula : str or dict of str to float

    Returns
    -------
    dict of str to float
        Element -> mass percent; the values sum to 100.

    Examples
    --------
    >>> {e: round(p, 2) for e, p in mass_percent("H2O").items()}
    {'H': 11.19, 'O': 88.81}
    """
    counts = parse_formula(formula) if isinstance(formula, str) else formula
    total = molar_mass(counts)
    return {e: 100.0 * get_element(e).atomic_mass * n / total for e, n in counts.items()}


def empirical_formula(mass_percents: dict, max_multiplier: int = 12, tol: float = 0.05) -> dict:
    """Empirical formula (smallest whole-number ratio) from percent composition by mass.

    Converts each mass to moles, divides by the smallest, then multiplies
    by the smallest integer ``m <= max_multiplier`` that brings every ratio
    within `tol` of a whole number -- the textbook procedure, which also
    recovers ratios such as 1:1.5 (x2) or 1:1.33 (x3).

    Parameters
    ----------
    mass_percents : dict of str to float
        Element -> mass (or mass percent; only ratios matter).
    max_multiplier : int, default 12
    tol : float, default 0.05
        Largest accepted distance from a whole number.

    Returns
    -------
    dict of str to int

    Raises
    ------
    ValueError
        If no multiplier up to `max_multiplier` gives whole numbers (usually
        imprecise analytical data; raise `tol`).

    Examples
    --------
    Glucose's composition gives CH2O:

    >>> empirical_formula({"C": 40.00, "H": 6.71, "O": 53.29})
    {'C': 1, 'H': 2, 'O': 1}
    >>> empirical_formula(mass_percent("Fe2O3"))
    {'Fe': 2, 'O': 3}
    """
    moles = {e: m / get_element(e).atomic_mass for e, m in mass_percents.items()}
    smallest = min(moles.values())
    if smallest <= 0:
        raise ValueError("mass percents must be positive")
    ratios = {e: n / smallest for e, n in moles.items()}
    for m in range(1, max_multiplier + 1):
        scaled = {e: r * m for e, r in ratios.items()}
        if all(abs(v - round(v)) <= tol for v in scaled.values()):
            return {e: int(round(v)) for e, v in scaled.items()}
    raise ValueError("no whole-number ratio found; data may be too imprecise (try a larger tol)")
