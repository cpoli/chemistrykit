r"""Thermodynamics of the fuel cell: reversible voltage and efficiency limits.

See Larminie & Dicks, *Fuel Cell Systems Explained*, 2nd ed., Ch. 2, or
Atkins & de Paula, *Physical Chemistry*, 11th ed., Topic 6D, for
:math:`\Delta G = -nFE` and its temperature dependence. The original
gas battery: W. R. Grove, *Phil. Mag.* 21, 417 (1842).
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import FARADAY

__all__ = [
    "reversible_cell_voltage",
    "reversible_cell_voltage_at_temperature",
    "fuel_cell_efficiency_limit",
]


def reversible_cell_voltage(delta_G: float, n: int, F: float = FARADAY) -> float:
    r"""Reversible (open-circuit) voltage from the reaction Gibbs energy: :math:`E = -\Delta G/(nF)`.

    Parameters
    ----------
    delta_G : float
        Gibbs energy of the overall cell reaction, in J/mol.
    n : int
        Electrons transferred per mole of reaction.
    F : float, default :data:`chemistrykit.constants.FARADAY`

    Returns
    -------
    float
        Reversible cell voltage, in V.

    Examples
    --------
    Hydrogen-oxygen cell, :math:`H_2 + \tfrac{1}{2}O_2 \to H_2O(l)`,
    :math:`\Delta G^\circ = -237.13` kJ/mol, n=2:

    >>> round(reversible_cell_voltage(-237.13e3, n=2), 3)
    1.229
    """
    return -delta_G / (n * F)


def reversible_cell_voltage_at_temperature(delta_H: float, delta_S: float, n: int, T, F: float = FARADAY):
    r"""Reversible voltage at temperature `T`: :math:`E(T) = -(\Delta H - T\Delta S)/(nF)`.

    Assumes :math:`\Delta H` and :math:`\Delta S` independent of
    temperature (a standard first approximation over modest ranges).

    Parameters
    ----------
    delta_H : float
        Reaction enthalpy, in J/mol.
    delta_S : float
        Reaction entropy, in J/(mol K).
    n : int
        Electrons transferred per mole of reaction.
    T : float or array-like of float
        Temperature, in K.
    F : float, default :data:`chemistrykit.constants.FARADAY`

    Returns
    -------
    float or ndarray
        Reversible cell voltage, in V.

    Examples
    --------
    The voltage falls with temperature when :math:`\Delta S<0`, with slope
    exactly :math:`\Delta S/(nF)`:

    >>> e1 = reversible_cell_voltage_at_temperature(-285.83e3, -163.3, 2, 300.0)
    >>> e2 = reversible_cell_voltage_at_temperature(-285.83e3, -163.3, 2, 301.0)
    >>> round((e2 - e1) * 2 * 96485.33212 / -163.3, 6)
    1.0
    """
    T = np.asarray(T, dtype=np.float64)
    result = -(delta_H - T * delta_S) / (n * F)
    return float(result) if result.ndim == 0 else result


def fuel_cell_efficiency_limit(delta_G: float, delta_H: float) -> float:
    r"""Maximum (thermodynamic) efficiency of a fuel cell: :math:`\eta_{max} = \Delta G/\Delta H`.

    The fraction of the fuel's heat of combustion that a reversible cell
    can deliver as electrical work -- not bounded by the Carnot factor of
    a heat engine (Larminie & Dicks, *Fuel Cell Systems Explained*, 2nd
    ed., Ch. 2.3).

    Parameters
    ----------
    delta_G, delta_H : float
        Reaction Gibbs energy and enthalpy (same units, same sign).

    Returns
    -------
    float

    Examples
    --------
    >>> round(fuel_cell_efficiency_limit(-237.13e3, -285.83e3), 3)
    0.83
    """
    return delta_G / delta_H
