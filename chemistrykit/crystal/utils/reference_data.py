r"""A handful of reference values for the lattice-energy/XRD examples and tests: Born exponents and ionic radii.

Plain reference data, in the same spirit as :mod:`chemistrykit.periodic_table`
being plain data rather than a dependency: not exhaustive, just enough to
exercise :mod:`chemistrykit.crystal.systems.lattice_energy` on the
handful of textbook salts used in this domain's tests and examples.

Born exponents (`BORN_EXPONENTS`) are tabulated by electronic
configuration, following M. Born and A. Lande's original assignment
(reproduced in Atkins & de Paula, *Physical Chemistry*, 11th ed., Table
20A.1, and West, *Solid State Chemistry and its Applications*, 2nd ed.,
Table 1.5): each noble-gas-like ion is assigned an exponent by which
closed shell it has, and a Born-Lande calculation on a salt of two
different ion types averages the two ions' exponents.

Ionic radii (`SHANNON_IONIC_RADII_PM`) are 6-coordinate (octahedral)
Shannon effective ionic radii (R. D. Shannon, *Acta Cryst.* A32, 751
(1976)), in picometers, for the handful of common main-group ions used
here -- an approximation in the sense that Shannon radii are themselves
coordination-number- and spin-state-dependent (the 6-coordinate value is
used throughout this module regardless of a given salt's actual
coordination, which is exact for rock-salt-structure salts like NaCl/MgO
and only approximate otherwise).
"""

from __future__ import annotations

__all__ = ["BORN_EXPONENTS", "average_born_exponent", "SHANNON_IONIC_RADII_PM"]

#: dict: Born exponent `n` by isoelectronic noble-gas configuration
#: (Atkins & de Paula, *Physical Chemistry*, 11th ed., Table 20A.1).
BORN_EXPONENTS: dict = {
    "He": 5,
    "Ne": 7,
    "Ar": 9,
    "Kr": 10,
    "Xe": 12,
}


def average_born_exponent(configuration_1: str, configuration_2: str) -> float:
    """Average the Born exponents of two ions' noble-gas configurations.

    The standard prescription for a binary salt of two different ion
    types (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 20A.1):
    take the arithmetic mean of each ion's tabulated exponent.

    Parameters
    ----------
    configuration_1, configuration_2 : str
        Noble-gas configuration labels, keys of :data:`BORN_EXPONENTS`
        (e.g. ``"Ne"`` for Na+, ``"Ar"`` for Cl-).

    Returns
    -------
    float

    Examples
    --------
    NaCl: Na+ is Ne-like (n=7), Cl- is Ar-like (n=9), average 8:

    >>> average_born_exponent("Ne", "Ar")
    8.0
    """
    return (BORN_EXPONENTS[configuration_1] + BORN_EXPONENTS[configuration_2]) / 2.0


#: dict: 6-coordinate (octahedral) Shannon effective ionic radii, in
#: picometers (R. D. Shannon, *Acta Cryst.* A32, 751 (1976)), for the
#: common ions used in this domain's Kapustinskii-equation examples/tests.
SHANNON_IONIC_RADII_PM: dict = {
    "Li+": 76.0,
    "Na+": 102.0,
    "K+": 138.0,
    "Rb+": 152.0,
    "Cs+": 167.0,
    "Mg2+": 72.0,
    "Ca2+": 100.0,
    "F-": 133.0,
    "Cl-": 181.0,
    "Br-": 196.0,
    "I-": 220.0,
    "O2-": 140.0,
}
