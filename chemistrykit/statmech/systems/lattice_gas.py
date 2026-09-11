r"""A canonical-ensemble lattice-gas model of adsorption.

The statistical-mechanical derivation of the Langmuir adsorption
isotherm: `M` identical, independent, non-interacting lattice sites, each
either empty or occupied by one adsorbed gas molecule with binding
(adsorption) energy :math:`\epsilon>0` released on binding. Treating the
sites as an open system in equilibrium with an ideal-gas reservoir at
pressure `P` and temperature `T` (a grand-canonical treatment of the
lattice, since molecules are exchanged with the gas phase) gives the
average fractional coverage exactly the Langmuir form -- see McQuarrie,
*Statistical Mechanics*, Ch. 8, or Hill, *An Introduction to Statistical
Thermodynamics* (1960), Ch. 7, for the full derivation. A fuller
empirical treatment of the Langmuir isotherm (fitting `K` directly to
adsorption data, plus Freundlich/BET isotherms) belongs in the future
``chemistrykit.surface``; this module is the microscopic derivation
*from* statistical mechanics, and :meth:`LatticeGasAdsorption.canonical_entropy`
is the purely combinatorial (fixed-`N`, canonical-ensemble) counting that
sits underneath it.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import K_B
from chemistrykit.statmech.utils.combinatorics import ln_binomial
from chemistrykit.statmech.utils.thermal_wavelength import thermal_de_broglie_wavelength

__all__ = ["LatticeGasAdsorption"]


class LatticeGasAdsorption:
    r"""Non-interacting lattice-gas model of adsorption on `M` independent sites.

    Each site's single-site grand partition function is
    :math:`\xi=1+\lambda e^{\beta\epsilon}` (empty, or occupied with
    Boltzmann factor :math:`e^{\beta\epsilon}` weighted by the gas-phase
    fugacity :math:`\lambda=e^{\beta\mu}`); since the `M` sites are
    independent, :math:`\Xi=\xi^M` and the average coverage is
    :math:`\theta=\langle N\rangle/M=1/(1+e^{-\beta(\mu+\epsilon)})`.
    Using the ideal-gas chemical potential
    :math:`\mu=k_BT\ln(P\Lambda^3/k_BT)` (with :math:`\Lambda` the
    thermal de Broglie wavelength,
    :func:`chemistrykit.statmech.utils.thermal_wavelength.thermal_de_broglie_wavelength`)
    turns this into exactly the Langmuir isotherm :math:`\theta=P/(P+P_0)`
    (see :meth:`coverage`, :meth:`p_half`).

    Parameters
    ----------
    adsorption_energy : float
        Binding energy per adsorbed molecule, in J (positive = favorable
        adsorption).
    mass : float
        Mass of the adsorbing gas molecule, in kg.
    T : float
        Absolute temperature, in K.
    """

    def __init__(self, adsorption_energy: float, mass: float, T: float):
        if mass <= 0 or T <= 0:
            raise ValueError("mass and T must be positive")
        self.adsorption_energy = float(adsorption_energy)
        self.mass = float(mass)
        self.T = float(T)

    def thermal_wavelength(self) -> float:
        """The gas molecule's thermal de Broglie wavelength at `T`. See :func:`chemistrykit.statmech.utils.thermal_wavelength.thermal_de_broglie_wavelength`."""
        return thermal_de_broglie_wavelength(self.mass, self.T)

    def p_half(self) -> float:
        r"""Pressure at half coverage, :math:`P_0=(k_BT/\Lambda^3)e^{-\epsilon/k_BT}`.

        Returns
        -------
        float
            Pressure, in Pa.
        """
        wavelength = self.thermal_wavelength()
        return (K_B * self.T / wavelength**3) * np.exp(-self.adsorption_energy / (K_B * self.T))

    def coverage(self, P):
        r"""Fractional coverage :math:`\theta(P) = P/(P+P_0)` -- the Langmuir isotherm.

        Parameters
        ----------
        P : float or array-like of float
            Gas pressure, in Pa.

        Returns
        -------
        float or ndarray

        Examples
        --------
        Coverage is exactly one-half at :math:`P=P_0` by construction,
        and approaches the physical limits of 0 (vacuum) and 1
        (saturation) at the extremes:

        >>> model = LatticeGasAdsorption(adsorption_energy=3.0e-20, mass=4.65e-26, T=300.0)
        >>> round(float(model.coverage(model.p_half())), 10)
        0.5
        >>> float(model.coverage(0.0))
        0.0
        >>> bool(model.coverage(1.0e15) > 0.9999)
        True
        """
        P = np.asarray(P, dtype=np.float64)
        P0 = self.p_half()
        return P / (P + P0)

    @staticmethod
    def canonical_entropy(N, M) -> float:
        r"""Canonical (fixed-`N`) combinatorial entropy :math:`S=k_B\ln\binom{M}{N}`.

        The exact configurational entropy of placing `N` indistinguishable
        adsorbed molecules on `M` distinguishable sites with no
        interactions, computed via the numerically exact
        :func:`chemistrykit.statmech.utils.combinatorics.ln_binomial`
        (rather than Stirling's approximation).

        Parameters
        ----------
        N : float
            Number of occupied sites (0 <= N <= M).
        M : float
            Total number of sites.

        Returns
        -------
        float
            Entropy, in J/K.

        Examples
        --------
        Zero at the two fully-ordered extremes (all sites empty, or all
        full -- only one microstate each), and maximal exactly at half
        filling, where it approaches the large-`M` Stirling estimate
        :math:`Mk_B\ln2` (McQuarrie, *Statistical Mechanics*, Ch. 8.1):

        >>> LatticeGasAdsorption.canonical_entropy(0, 100)
        0.0
        >>> LatticeGasAdsorption.canonical_entropy(100, 100)
        0.0
        >>> import numpy as np
        >>> M = 100_000
        >>> S_half = LatticeGasAdsorption.canonical_entropy(M // 2, M)
        >>> relative_error = abs(S_half - M * 1.380649e-23 * np.log(2)) / (M * 1.380649e-23 * np.log(2))
        >>> bool(relative_error < 1e-4)
        True
        """
        return float(K_B * ln_binomial(M, N))
