"""Chemistry constants shared across chemistrykit, in SI units.

Every subpackage in chemistrykit computes with the actual numeric value of
the physical constants it needs -- unlike physicskit, no chemistrykit
subpackage adopts a rescaled unit system (there is no chemistry analogue of
setting :math:`\\hbar = 1`), so this module is simply *the* place every
subpackage gets ``R``, ``NA``, etc. from, so that "R" means the same number
everywhere in chemistrykit rather than being independently retyped -- and
potentially drifting -- in each subpackage. Values are taken from
:mod:`scipy.constants` (2019 SI redefinition / CODATA), which is already a
chemistrykit dependency.

A handful of conversion helpers are also provided, but -- following
physicskit's own restraint here -- only for conversions that are
well-defined without an extra domain choice: energy <-> temperature (via
:math:`k_B`) and atm <-> Pa. Converting a concentration to, say, a mole
fraction or an activity requires picking a solvent/solution model (density,
molar volume, activity-coefficient model), so no generic helper is provided
for that -- see the relevant subpackage (:mod:`chemistrykit.solutions`,
:mod:`chemistrykit.thermo`) instead.

See also :mod:`chemistrykit.periodic_table` for the accompanying periodic-
table data (re-exported here for convenience).
"""

from __future__ import annotations

import scipy.constants as _sc

from chemistrykit.periodic_table import PERIODIC_TABLE, Element, get_element, molar_mass

__all__ = [
    "R",
    "NA",
    "K_B",
    "H",
    "HBAR",
    "C",
    "ELEMENTARY_CHARGE",
    "FARADAY",
    "ATOMIC_MASS_UNIT",
    "ELECTRON_MASS",
    "PROTON_MASS",
    "NEUTRON_MASS",
    "VACUUM_PERMITTIVITY",
    "STANDARD_TEMPERATURE",
    "STP_TEMPERATURE",
    "STANDARD_PRESSURE",
    "STANDARD_PRESSURE_BAR",
    "STANDARD_PRESSURE_ATM",
    "ELECTRONVOLT",
    "energy_to_temperature",
    "temperature_to_energy",
    "joules_to_ev",
    "ev_to_joules",
    "atm_to_pa",
    "pa_to_atm",
    "celsius_to_kelvin",
    "kelvin_to_celsius",
    "Element",
    "PERIODIC_TABLE",
    "get_element",
    "molar_mass",
]

#: Molar gas constant, in J mol^-1 K^-1.
R = _sc.R

#: Avogadro constant, in mol^-1 (exact by SI definition).
NA = _sc.N_A

#: Boltzmann constant, in J/K (exact by SI definition). ``R == NA * K_B``.
K_B = _sc.k

#: Planck constant, in J*s (exact by SI definition).
H = _sc.h

#: Reduced Planck constant, in J*s.
HBAR = _sc.hbar

#: Speed of light in vacuum, in m/s (exact by SI definition).
C = _sc.c

#: Elementary charge, in C (exact by SI definition).
ELEMENTARY_CHARGE = _sc.e

#: Faraday constant, in C/mol -- the charge of one mole of electrons,
#: ``NA * ELEMENTARY_CHARGE``. Central to :mod:`chemistrykit.electrochem`.
FARADAY = _sc.N_A * _sc.e

#: Unified atomic mass unit, in kg (1/12 the mass of carbon-12).
ATOMIC_MASS_UNIT = _sc.atomic_mass

#: Electron rest mass, in kg.
ELECTRON_MASS = _sc.m_e

#: Proton rest mass, in kg.
PROTON_MASS = _sc.m_p

#: Neutron rest mass, in kg.
NEUTRON_MASS = _sc.m_n

#: Vacuum electric permittivity, in F/m.
VACUUM_PERMITTIVITY = _sc.epsilon_0

#: One electronvolt, in J.
ELECTRONVOLT = _sc.eV

#: Standard thermodynamic-state temperature, 298.15 K (25 degC) -- the
#: reference temperature for tabulated standard-state values (e.g.
#: standard Gibbs energies, standard reduction potentials), per IUPAC.
STANDARD_TEMPERATURE = 298.15

#: "Standard temperature" in the STP (Standard Temperature and Pressure)
#: sense used for ideal-gas molar-volume calculations, 273.15 K (0 degC) --
#: a different conventional reference point from STANDARD_TEMPERATURE
#: above; the two are easy to conflate and chemistrykit keeps them as
#: separate named constants rather than one ambiguous "standard temperature".
STP_TEMPERATURE = 273.15

#: Standard pressure per the current (1982-) IUPAC convention, 1 bar =
#: 10^5 Pa exactly.
STANDARD_PRESSURE_BAR = 1.0e5

#: Standard pressure per the older (pre-1982) convention, still in common
#: use (e.g. US engineering practice): 1 atm, in Pa.
STANDARD_PRESSURE_ATM = _sc.atm

#: chemistrykit's default "standard pressure", matching current IUPAC
#: practice (1 bar). Functions that need the historical 1-atm convention
#: instead take it as an explicit argument or use STANDARD_PRESSURE_ATM.
STANDARD_PRESSURE = STANDARD_PRESSURE_BAR


def energy_to_temperature(energy_j):
    """Convert an energy (in J) to the temperature (in K) with :math:`E = k_B T`.

    Examples
    --------
    >>> round(float(energy_to_temperature(K_B)), 6)
    1.0
    """
    return energy_j / K_B


def temperature_to_energy(temperature_k):
    """Convert a temperature (in K) to an energy (in J) with :math:`E = k_B T`.

    Examples
    --------
    >>> temperature_to_energy(1.0) == K_B
    True
    """
    return temperature_k * K_B


def joules_to_ev(energy_j):
    """Convert an energy from joules to electronvolts."""
    return energy_j / ELECTRONVOLT


def ev_to_joules(energy_ev):
    """Convert an energy from electronvolts to joules."""
    return energy_ev * ELECTRONVOLT


def atm_to_pa(pressure_atm):
    """Convert a pressure from standard atmospheres to pascals.

    Examples
    --------
    >>> round(atm_to_pa(1.0), 2)
    101325.0
    """
    return pressure_atm * STANDARD_PRESSURE_ATM


def pa_to_atm(pressure_pa):
    """Convert a pressure from pascals to standard atmospheres.

    Examples
    --------
    >>> round(pa_to_atm(101325.0), 6)
    1.0
    """
    return pressure_pa / STANDARD_PRESSURE_ATM


def celsius_to_kelvin(temperature_c):
    """Convert a temperature from degrees Celsius to kelvin.

    Examples
    --------
    >>> celsius_to_kelvin(0.0)
    273.15
    """
    return temperature_c + STP_TEMPERATURE


def kelvin_to_celsius(temperature_k):
    """Convert a temperature from kelvin to degrees Celsius.

    Examples
    --------
    >>> kelvin_to_celsius(298.15)
    25.0
    """
    return temperature_k - STP_TEMPERATURE
