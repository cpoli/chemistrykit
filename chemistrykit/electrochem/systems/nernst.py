r"""The Nernst equation: standard cells, concentration cells, and activity-corrected cell potentials.

See Bard & Faulkner, *Electrochemical Methods: Fundamentals and
Applications*, 2nd ed., Ch. 2.1, or Atkins & de Paula, *Physical
Chemistry*, 11th ed., Ch. 6.10, throughout. The original: W. Nernst,
*Z. Phys. Chem.* 4, 129 (1889).
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import FARADAY, STANDARD_TEMPERATURE, R
from chemistrykit.solutions.systems.activity import (
    activity_coefficient_debye_huckel_extended,
    activity_coefficient_debye_huckel_limiting,
    ionic_strength,
)

__all__ = [
    "nernst_potential",
    "concentration_cell_potential",
    "activity_corrected_reaction_quotient",
    "nernst_potential_with_activity",
]


def nernst_potential(E_standard: float, n: int, Q, T: float = STANDARD_TEMPERATURE, R_gas: float = R, F: float = FARADAY):
    r"""The Nernst equation: :math:`E = E^\circ - \frac{RT}{nF}\ln Q`.

    Relates a cell's potential under arbitrary (non-standard-state)
    conditions to its standard potential and the reaction quotient `Q`
    of the overall cell reaction (Bard & Faulkner, *Electrochemical
    Methods*, 2nd ed., Ch. 2.1, eq. 2.1.13).

    Parameters
    ----------
    E_standard : float
        Standard cell potential :math:`E^\circ`, in V.
    n : int
        Number of electrons transferred in the balanced overall cell
        reaction (see :func:`chemistrykit.electrochem.systems.standard_potentials.balance_redox_reaction`).
    Q : float or array-like of float
        Reaction quotient of the overall cell reaction (activities of
        products over reactants, each raised to its stoichiometric
        coefficient).
    T : float, default :data:`chemistrykit.constants.STANDARD_TEMPERATURE`
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.
    F : float, default :data:`chemistrykit.constants.FARADAY`
        Faraday constant, in C/mol.

    Returns
    -------
    float or ndarray

    Examples
    --------
    At :math:`Q = 1` (unit activity of every species -- the standard
    state by definition), the Nernst equation reduces exactly to
    :math:`E^\circ`:

    >>> round(float(nernst_potential(E_standard=1.10, n=2, Q=1.0)), 6)
    1.1

    Increasing the product-side activity (`Q > 1`) lowers the cell
    potential below standard, per Le Chatelier's principle applied to
    the driving reaction:

    >>> bool(nernst_potential(1.10, 2, Q=10.0) < 1.10)
    True
    """
    Q = np.asarray(Q, dtype=np.float64)
    result = E_standard - (R_gas * T) / (n * F) * np.log(Q)
    return float(result) if result.ndim == 0 else result


def concentration_cell_potential(n: int, C_cathode, C_anode, T: float = STANDARD_TEMPERATURE, R_gas: float = R, F: float = FARADAY):
    r"""Cell potential of a concentration cell: :math:`E = \frac{RT}{nF}\ln(C_{cathode}/C_{anode})`.

    A concentration cell pairs two half-cells built from the *same*
    electrode material and redox couple, differing only in the
    concentration of the dissolved species -- so :math:`E^\circ = 0`
    identically (both half-reactions have the same standard potential)
    and the entire driving force comes from the Nernst equation's
    activity term (Atkins & de Paula, *Physical Chemistry*, 11th ed.,
    Ch. 6.10). By convention the more concentrated compartment is the
    cathode (reduction is favored there, consuming ions and diluting
    that side, driving the system toward equalized concentrations).

    Parameters
    ----------
    n : int
        Number of electrons transferred per formula unit of the redox
        couple (e.g. 2 for :math:`M^{2+} + 2e^- \to M`).
    C_cathode, C_anode : float or array-like of float
        Ion concentrations in the cathode and anode compartments, in the
        same units (an activity-coefficient-corrected value is more
        accurate at high concentration -- see
        :func:`nernst_potential_with_activity`).
    T : float, default :data:`chemistrykit.constants.STANDARD_TEMPERATURE`
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.
    F : float, default :data:`chemistrykit.constants.FARADAY`
        Faraday constant, in C/mol.

    Returns
    -------
    float or ndarray
        Cell potential, in V. Positive when `C_cathode > C_anode`, zero
        at equal concentrations (no net driving force between two
        identical half-cells).

    Examples
    --------
    Equal concentrations give exactly zero cell potential:

    >>> round(float(concentration_cell_potential(n=2, C_cathode=0.10, C_anode=0.10)), 9)
    0.0

    A tenfold concentration difference in a 1-electron couple at 25 degC
    gives the textbook ~59 mV per decade:

    >>> round(float(concentration_cell_potential(n=1, C_cathode=0.10, C_anode=0.010)), 4)
    0.0592
    """
    return nernst_potential(E_standard=0.0, n=n, Q=np.asarray(C_anode, dtype=np.float64) / np.asarray(C_cathode, dtype=np.float64), T=T, R_gas=R_gas, F=F)


def activity_corrected_reaction_quotient(
    concentrations, charges, stoich_coeffs, T: float = STANDARD_TEMPERATURE, extended: bool = True, Ba: float = 1.0
) -> float:
    r"""Reaction quotient built from Debye-Huckel activity-corrected concentrations rather than raw concentrations.

    :math:`Q = \prod_i a_i^{\nu_i} = \prod_i (\gamma_i c_i / c^\circ)^{\nu_i}`
    (with :math:`c^\circ = 1` mol/L implicit throughout, matching the
    convention of :func:`chemistrykit.thermo.systems.equilibrium.reaction_quotient`),
    where every ion's activity coefficient :math:`\gamma_i` is evaluated
    at the *mixture's* ionic strength via
    :func:`chemistrykit.solutions.systems.activity.ionic_strength` and
    :func:`chemistrykit.solutions.systems.activity.activity_coefficient_debye_huckel_extended`
    -- reusing this package's existing Debye-Huckel machinery rather than
    assuming unit activity coefficients, per Bard & Faulkner,
    *Electrochemical Methods*, 2nd ed., Ch. 2.1.3. Only ionic species
    (:math:`z\neq0`) get a nontrivial activity coefficient; a neutral
    species or a pure solid/liquid passed with `charges=0` is treated as
    ideal (:math:`\gamma=1`), the usual convention for a species at unit
    activity by definition (e.g. the solid electrode itself).

    Parameters
    ----------
    concentrations : array-like of float
        Concentration of each species, in mol/L.
    charges : array-like of float
        Charge number of each species (0 for a neutral species / pure
        solid or liquid).
    stoich_coeffs : array-like of float
        Signed net stoichiometric coefficient of each species (positive
        for products, negative for reactants), same order as
        `concentrations`.
    T : float, default :data:`chemistrykit.constants.STANDARD_TEMPERATURE`
        Absolute temperature; only used insofar as
        :data:`chemistrykit.solutions.systems.activity.DEBYE_HUCKEL_A_25C`
        is itself a 25 degC value -- this function does not otherwise
        adjust `A` for `T`, an approximation reasonable near room
        temperature only.
    extended : bool, default True
        Use the extended Debye-Huckel law (valid to higher ionic
        strength) rather than the limiting law.
    Ba : float, default 1.0
        Ion-size parameter for the extended law; see
        :func:`chemistrykit.solutions.systems.activity.activity_coefficient_debye_huckel_extended`.

    Returns
    -------
    float

    Examples
    --------
    At vanishing ionic strength (dilute limit), activity coefficients
    are all 1 and this reduces to the raw-concentration reaction
    quotient:

    >>> c = [1e-9, 1e-12]
    >>> z = [1, -2]
    >>> nu = [-1.0, 1.0]
    >>> round(activity_corrected_reaction_quotient(c, z, nu), 6)
    0.001

    At higher ionic strength the two differently-charged species'
    activity coefficients diverge from each other and from 1, so the
    activity-corrected quotient deviates from the raw-concentration
    ratio (here :math:`c_1/c_0=0.1`):

    >>> Q_raw = 0.1
    >>> Q_corrected = activity_corrected_reaction_quotient([0.1, 0.1], [1, -2], nu)
    >>> bool(abs(Q_corrected - Q_raw) > 1e-3)
    True
    """
    c = np.asarray(concentrations, dtype=np.float64)
    z = np.asarray(charges, dtype=np.float64)
    nu = np.asarray(stoich_coeffs, dtype=np.float64)
    I = ionic_strength(c, z)
    gammas = np.array(
        [
            1.0 if zi == 0.0 else (activity_coefficient_debye_huckel_extended(zi, I, Ba=Ba) if extended else activity_coefficient_debye_huckel_limiting(zi, I))
            for zi in z
        ]
    )
    activities = gammas * c
    return float(np.prod(activities**nu))


def nernst_potential_with_activity(
    E_standard: float,
    n: int,
    concentrations,
    charges,
    stoich_coeffs,
    T: float = STANDARD_TEMPERATURE,
    extended: bool = True,
    Ba: float = 1.0,
) -> float:
    r"""The Nernst equation with `Q` built from Debye-Huckel activity-corrected concentrations.

    Combines :func:`activity_corrected_reaction_quotient` with
    :func:`nernst_potential` -- the "reuse `chemistrykit.solutions`'s
    activity-coefficient machinery" version of the Nernst equation, as
    opposed to :func:`nernst_potential`'s assumption of ideal (unit
    activity coefficient) behavior. Real electrochemical cells, like
    real solution equilibria, deviate from ideal Nernstian behavior at
    concentrations much above dilute (Bard & Faulkner, *Electrochemical
    Methods*, 2nd ed., Ch. 2.1.3).

    Parameters
    ----------
    E_standard : float
        Standard cell potential, in V.
    n : int
        Electrons transferred in the balanced overall cell reaction.
    concentrations, charges, stoich_coeffs : array-like of float
        As in :func:`activity_corrected_reaction_quotient`.
    T : float, default :data:`chemistrykit.constants.STANDARD_TEMPERATURE`
    extended : bool, default True
    Ba : float, default 1.0

    Returns
    -------
    float

    Examples
    --------
    In the dilute limit, this recovers the ideal Nernst potential to
    high accuracy:

    >>> c = [1e-6, 1e-12]
    >>> z = [2, 0]
    >>> nu = [-1.0, 1.0]
    >>> Q_raw = c[1] / c[0]
    >>> ideal = nernst_potential(0.34, n=2, Q=Q_raw)
    >>> corrected = nernst_potential_with_activity(0.34, n=2, concentrations=c, charges=z, stoich_coeffs=nu)
    >>> bool(abs(ideal - corrected) < 1e-4)
    True
    """
    Q = activity_corrected_reaction_quotient(concentrations, charges, stoich_coeffs, T=T, extended=extended, Ba=Ba)
    return nernst_potential(E_standard, n, Q, T=T)
