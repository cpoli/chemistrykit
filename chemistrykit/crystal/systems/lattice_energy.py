r"""Ionic-crystal lattice energy: the Born-Lande and Kapustinskii equations.

Two independently-derived ways to estimate the same physical quantity --
the energy released forming 1 mol of an ionic solid from its gaseous ions
-- from different inputs, so they are usefully compared side by side on
the shared :class:`chemistrykit.crystal.core.base_system.LatticeEnergyModel`
interface. Both are approximations (ionic, point-charge, purely
electrostatic-plus-repulsion models of what is, in reality, a partly
covalent bonding picture for many real salts); see Atkins & de Paula,
*Physical Chemistry*, 11th ed., Ch. 20A, for the standard treatment of
both, and West, *Solid State Chemistry and its Applications*, 2nd ed.
(2014), Ch. 1.4, for the Born-Lande derivation in particular.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import ELEMENTARY_CHARGE, NA, VACUUM_PERMITTIVITY
from chemistrykit.crystal.core.base_system import LatticeEnergyModel
from chemistrykit.crystal.systems.madelung import madelung_constant_nacl

__all__ = ["BornLandeLatticeEnergy", "KapustinskiiLatticeEnergy"]


class BornLandeLatticeEnergy(LatticeEnergyModel):
    r"""The Born-Lande equation for ionic-crystal lattice energy.

    .. math::

        U = -\frac{N_A M |z_+z_-|e^2}{4\pi\varepsilon_0 r_0}\left(1-\frac1n\right)

    The Madelung constant `M` accounts for the electrostatic sum over the
    whole lattice; the Born-repulsion correction :math:`(1-1/n)` accounts
    for short-range electron-cloud repulsion at the equilibrium spacing
    `r0` (`n` the Born exponent, from
    :func:`chemistrykit.crystal.utils.reference_data.average_born_exponent`)
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., eq. 20A.5).

    Parameters
    ----------
    cation_charge : int
        Cation charge number (e.g. 1 for Na+, 2 for Mg2+); sign is
        ignored (absolute value used).
    anion_charge : int
        Anion charge number (e.g. 1 for Cl-, 2 for O2-); sign is ignored.
    r0 : float
        Equilibrium nearest-neighbor (cation-anion) distance, in meters.
    born_exponent : float
        Born exponent `n` (typically 5-12; see
        :data:`chemistrykit.crystal.utils.reference_data.BORN_EXPONENTS`).
    madelung_constant : float, default :func:`chemistrykit.crystal.systems.madelung.madelung_constant_nacl`
        Madelung constant of the crystal structure (the default is only
        appropriate for a rock-salt-structure salt; other structures
        (CsCl, fluorite, ...) have different Madelung constants not
        computed by this package, and must be supplied explicitly).

    Examples
    --------
    NaCl (:math:`r_0=282\,\text{pm}`, Born exponent 8 -- the Ne/Ar
    average): the Born-Lande estimate comes out close to, but (as
    expected for this level of approximation) somewhat less negative
    than, the experimental Born-Haber-cycle value of about
    :math:`-787\,\text{kJ/mol}`:

    >>> model = BornLandeLatticeEnergy(cation_charge=1, anion_charge=1, r0=282e-12, born_exponent=8.0)
    >>> round(model.lattice_energy() / 1000.0, 1)  # kJ/mol
    -753.4
    """

    def __init__(self, cation_charge: int, anion_charge: int, r0: float, born_exponent: float, madelung_constant: float | None = None):
        if r0 <= 0:
            raise ValueError("r0 must be positive")
        if born_exponent <= 1:
            raise ValueError("born_exponent must be > 1")
        self.cation_charge = abs(int(cation_charge))
        self.anion_charge = abs(int(anion_charge))
        self.r0 = float(r0)
        self.born_exponent = float(born_exponent)
        self.madelung_constant = float(madelung_constant) if madelung_constant is not None else madelung_constant_nacl()

    def lattice_energy(self) -> float:
        charge_product = self.cation_charge * self.anion_charge
        numerator = NA * self.madelung_constant * charge_product * ELEMENTARY_CHARGE**2
        prefactor = numerator / (4.0 * np.pi * VACUUM_PERMITTIVITY * self.r0)
        return -prefactor * (1.0 - 1.0 / self.born_exponent)


class KapustinskiiLatticeEnergy(LatticeEnergyModel):
    r"""The Kapustinskii equation: an ionic-radii-only estimate of lattice energy, no Madelung constant required.

    .. math::

        U = -\frac{\kappa\,\nu\,|z_+z_-|}{r_++r_-}\left(1-\frac{d}{r_++r_-}\right)

    with :math:`\kappa=1.2025\times10^5\,\text{kJ pm mol}^{-1}` and
    :math:`d=34.5\,\text{pm}` empirical constants fitted so this
    reproduces experimental/Born-Lande lattice energies without needing a
    structure-specific Madelung constant -- useful for salts (or
    hypothetical salts) whose crystal structure is unknown (A. F.
    Kapustinskii, *Q. Rev. Chem. Soc.* 10, 283 (1956); Atkins & de Paula,
    *Physical Chemistry*, 11th ed., Ch. 20A.2). This is itself an
    approximation, typically agreeing with more detailed calculations to
    within about 5%.

    Parameters
    ----------
    n_ions : int
        Number of ions per formula unit (e.g. 2 for NaCl, 3 for CaF2/Na2O).
    cation_charge : int
        Cation charge number; sign ignored.
    anion_charge : int
        Anion charge number; sign ignored.
    r_cation_pm : float
        Cation radius, in picometers (e.g. a Shannon effective ionic
        radius,
        :data:`chemistrykit.crystal.utils.reference_data.SHANNON_IONIC_RADII_PM`).
    r_anion_pm : float
        Anion radius, in picometers.

    Examples
    --------
    NaCl again, this time from ionic radii alone (no Madelung constant,
    no crystal structure assumed) -- close to, but not identical to, the
    Born-Lande estimate above, as expected of two different
    approximations to the same quantity:

    >>> model = KapustinskiiLatticeEnergy(n_ions=2, cation_charge=1, anion_charge=1, r_cation_pm=102.0, r_anion_pm=181.0)
    >>> round(model.lattice_energy() / 1000.0, 1)  # kJ/mol
    -746.2
    """

    _KAPPA_KJ_PM_PER_MOL = 1.2025e5
    _D_PM = 34.5

    def __init__(self, n_ions: int, cation_charge: int, anion_charge: int, r_cation_pm: float, r_anion_pm: float):
        if r_cation_pm <= 0 or r_anion_pm <= 0:
            raise ValueError("ionic radii must be positive")
        if n_ions < 2:
            raise ValueError("n_ions (ions per formula unit) must be at least 2")
        self.n_ions = int(n_ions)
        self.cation_charge = abs(int(cation_charge))
        self.anion_charge = abs(int(anion_charge))
        self.r_cation_pm = float(r_cation_pm)
        self.r_anion_pm = float(r_anion_pm)

    def lattice_energy(self) -> float:
        r0_pm = self.r_cation_pm + self.r_anion_pm
        u_kj_per_mol = -(self._KAPPA_KJ_PM_PER_MOL * self.n_ions * self.cation_charge * self.anion_charge / r0_pm) * (1.0 - self._D_PM / r0_pm)
        return u_kj_per_mol * 1000.0
