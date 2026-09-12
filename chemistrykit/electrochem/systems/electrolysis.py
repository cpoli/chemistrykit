r"""Galvanic vs. electrolytic cells, and Faraday's laws of electrolysis.

See Bard & Faulkner, *Electrochemical Methods: Fundamentals and
Applications*, 2nd ed., Ch. 1.3.3, or Atkins & de Paula, *Physical
Chemistry*, 11th ed., Ch. 6.9-6.10, throughout. The original: M. Faraday,
*Experimental Researches in Electricity*, Series VII (1834).
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import FARADAY

__all__ = [
    "charge_from_current",
    "moles_from_charge",
    "mass_from_charge",
    "faradays_law_mass",
    "minimum_applied_voltage_electrolytic",
]


def charge_from_current(current, time):
    r"""Charge passed at constant current: :math:`Q = It`.

    Parameters
    ----------
    current : float or array-like of float
        Current, in A.
    time : float or array-like of float
        Time, in s.

    Returns
    -------
    float or ndarray
        Charge, in C.

    Examples
    --------
    >>> round(float(charge_from_current(current=2.0, time=3600.0)), 2)
    7200.0
    """
    result = np.asarray(current, dtype=np.float64) * np.asarray(time, dtype=np.float64)
    return float(result) if result.ndim == 0 else result


def moles_from_charge(charge, n: int, F: float = FARADAY):
    r"""Faraday's first law: moles of substance transformed :math:`= Q/(nF)`.

    One mole of electrons (:math:`F` coulombs) transforms exactly
    :math:`1/n` mole of a species undergoing an `n`-electron
    electrode reaction -- charge passed is directly proportional to the
    amount of chemical change, the historically first-discovered
    electrochemical stoichiometry relationship (Faraday, 1834; Bard &
    Faulkner, *Electrochemical Methods*, 2nd ed., Ch. 1.3.3).

    Parameters
    ----------
    charge : float or array-like of float
        Charge passed, in C.
    n : int
        Electrons transferred per formula unit of the species.
    F : float, default :data:`chemistrykit.constants.FARADAY`
        Faraday constant, in C/mol.

    Returns
    -------
    float or ndarray
        Amount transformed, in mol.

    Examples
    --------
    Doubling the charge exactly doubles the moles transformed (Faraday's
    first law's defining linearity):

    >>> moles_from_charge(2.0 * 96485.332, n=1) == 2.0 * moles_from_charge(96485.332, n=1)
    True
    """
    result = np.asarray(charge, dtype=np.float64) / (n * F)
    return float(result) if result.ndim == 0 else result


def mass_from_charge(charge, molar_mass: float, n: int, F: float = FARADAY):
    r"""Mass of substance deposited/consumed: :math:`m = QM/(nF)`.

    Parameters
    ----------
    charge : float or array-like of float
        Charge passed, in C.
    molar_mass : float
        Molar mass of the species, in g/mol.
    n : int
        Electrons transferred per formula unit.
    F : float, default :data:`chemistrykit.constants.FARADAY`
        Faraday constant, in C/mol.

    Returns
    -------
    float or ndarray
        Mass, in g.

    Examples
    --------
    Silver electroplating, :math:`Ag^+ + e^- \to Ag` (n=1, M=107.87 g/mol):
    depositing 1 mol of electrons' worth of charge deposits exactly one
    mole (107.87 g) of silver:

    >>> round(mass_from_charge(charge=96485.332, molar_mass=107.87, n=1), 2)
    107.87
    """
    return moles_from_charge(charge, n, F) * molar_mass


def faradays_law_mass(current, time, molar_mass: float, n: int, efficiency: float = 1.0, F: float = FARADAY):
    r"""Mass deposited/consumed by electrolysis at constant current: :math:`m = ItM\eta/(nF)`.

    Combines :func:`charge_from_current` and :func:`mass_from_charge`
    (Faraday's laws of electrolysis together, Bard & Faulkner,
    *Electrochemical Methods*, 2nd ed., Ch. 1.3.3), with an optional
    current (faradaic) efficiency :math:`\eta\in(0,1]` for the common
    real-world case where a side reaction (e.g. competing gas evolution)
    consumes some of the current without depositing the desired product.

    Parameters
    ----------
    current : float
        Applied current, in A.
    time : float
        Duration, in s.
    molar_mass : float
        Molar mass of the deposited/consumed species, in g/mol.
    n : int
        Electrons transferred per formula unit.
    efficiency : float, default 1.0
        Current (faradaic) efficiency, in :math:`(0, 1]`; 1.0 is the
        ideal, no-side-reaction case.
    F : float, default :data:`chemistrykit.constants.FARADAY`
        Faraday constant, in C/mol.

    Returns
    -------
    float
        Mass, in g.

    Raises
    ------
    ValueError
        If `efficiency` is not in :math:`(0, 1]`.

    Examples
    --------
    Copper electrorefining, :math:`Cu^{2+}+2e^-\to Cu` (M=63.546 g/mol,
    n=2), run at 10 A for 1 hour: mass scales linearly with both current
    and time, the direct experimental test of Faraday's first law.

    >>> m1 = faradays_law_mass(current=10.0, time=3600.0, molar_mass=63.546, n=2)
    >>> m2 = faradays_law_mass(current=20.0, time=3600.0, molar_mass=63.546, n=2)
    >>> round(m2 / m1, 6)
    2.0
    >>> round(m1, 4)
    11.8549

    Reduced current efficiency proportionally reduces the mass deposited:

    >>> m_ideal = faradays_law_mass(10.0, 3600.0, 63.546, n=2, efficiency=1.0)
    >>> m_90pct = faradays_law_mass(10.0, 3600.0, 63.546, n=2, efficiency=0.9)
    >>> round(m_90pct / m_ideal, 6)
    0.9
    """
    if not (0.0 < efficiency <= 1.0):
        raise ValueError("efficiency must be in (0, 1]")
    Q = charge_from_current(current, time)
    return mass_from_charge(Q, molar_mass, n, F) * efficiency


def minimum_applied_voltage_electrolytic(E_cell: float) -> float:
    r"""Minimum externally applied voltage to drive a non-spontaneous (electrolytic) cell reaction.

    A cell reaction with :math:`E_{cell}<0` is non-spontaneous as
    written and will not run galvanically; forcing it to run
    (electrolysis) requires an externally applied voltage of at least
    :math:`|E_{cell}|`, opposing the cell's natural (spontaneous
    reverse) direction (Atkins & de Paula, *Physical Chemistry*, 11th
    ed., Ch. 6.10). **Approximation flagged**: this is the
    thermodynamic minimum only -- a real electrolysis cell also needs
    to overcome activation overpotential at each electrode (see
    :mod:`chemistrykit.electrochem.systems.butler_volmer`) and ohmic
    (solution-resistance) losses, so the practically applied voltage is
    always somewhat larger than this minimum.

    Parameters
    ----------
    E_cell : float
        Cell potential of the reaction as written, in V (should be
        negative -- i.e. the reaction is non-spontaneous -- for
        electrolysis to be the relevant regime; see
        :func:`chemistrykit.electrochem.systems.standard_potentials.is_spontaneous`).

    Returns
    -------
    float
        Minimum applied voltage magnitude, in V.

    Examples
    --------
    Electrolyzing molten NaCl (reversing the spontaneous
    Na/Cl2 cell reaction, :math:`E_{cell}=-4.07` V as written for
    Na deposition):

    >>> round(minimum_applied_voltage_electrolytic(-4.07), 2)
    4.07
    """
    return abs(E_cell)
