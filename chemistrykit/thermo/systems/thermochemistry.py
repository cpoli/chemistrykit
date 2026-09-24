r"""Thermochemistry: Hess's law of constant heat summation.

Because enthalpy is a state function, the enthalpy change of a reaction
is the same whether it happens in one step or several (G. H. Hess,
"Thermochemische Untersuchungen," *Ann. Phys. Chem.* 50, 385-404
(1840); Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 2C).
Two consequences are implemented here: a target reaction's enthalpy is
the weighted sum of the enthalpies of any steps that add up to it, and a
reaction's standard enthalpy follows from tabulated standard enthalpies
of formation.
"""

from __future__ import annotations

import numpy as np

__all__ = ["hess_law_enthalpy", "reaction_enthalpy_from_formation"]


def hess_law_enthalpy(multipliers, step_enthalpies) -> float:
    r"""Enthalpy of a reaction built as a linear combination of steps, :math:`\Delta H = \sum_k c_k\,\Delta H_k`.

    Parameters
    ----------
    multipliers : array-like of float
        Coefficient :math:`c_k` applied to each step reaction (negative to
        reverse a step, 2 to double it, and so on).
    step_enthalpies : array-like of float
        Enthalpy change of each step as written, in J/mol (or any
        consistent unit).

    Returns
    -------
    float

    Examples
    --------
    Burning graphite to CO (-110.5 kJ/mol) and then CO to CO2 (-283.0
    kJ/mol) releases the same heat as burning graphite straight to CO2:

    >>> hess_law_enthalpy([1.0, 1.0], [-110.5, -283.0])
    -393.5
    """
    c = np.asarray(multipliers, dtype=np.float64)
    dH = np.asarray(step_enthalpies, dtype=np.float64)
    if c.shape != dH.shape:
        raise ValueError("multipliers and step_enthalpies must have the same shape")
    return float(np.sum(c * dH))


def reaction_enthalpy_from_formation(stoich_coeffs, formation_enthalpies) -> float:
    r"""Standard reaction enthalpy from enthalpies of formation, :math:`\Delta_rH^\circ = \sum_i \nu_i\,\Delta_fH_i^\circ`.

    This is Hess's law applied to the cycle "decompose reactants into
    their elements, then build the products from them".

    Parameters
    ----------
    stoich_coeffs : array-like of float
        Signed stoichiometric coefficients (negative for reactants).
    formation_enthalpies : array-like of float
        Standard enthalpy of formation of each species, same order, in
        J/mol or kJ/mol (0 for elements in their reference state).

    Returns
    -------
    float

    Examples
    --------
    Combustion of methane, :math:`CH_4 + 2O_2 \to CO_2 + 2H_2O(l)`, with
    formation enthalpies -74.8, 0, -393.5 and -285.8 kJ/mol:

    >>> round(reaction_enthalpy_from_formation([-1, -2, 1, 2], [-74.8, 0.0, -393.5, -285.8]), 1)
    -890.3
    """
    nu = np.asarray(stoich_coeffs, dtype=np.float64)
    dHf = np.asarray(formation_enthalpies, dtype=np.float64)
    if nu.shape != dHf.shape:
        raise ValueError("stoich_coeffs and formation_enthalpies must have the same shape")
    return float(np.sum(nu * dHf))
