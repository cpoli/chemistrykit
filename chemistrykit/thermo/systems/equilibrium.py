r"""Chemical reaction equilibrium: Kp/Kc, Q, van't Hoff, and Gibbs-energy minimization.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 6 ("Chemical
equilibrium") for the reaction quotient, :math:`K_p`/:math:`K_c`, and the
van't Hoff equation, and Smith, Van Ness & Abbott, *Introduction to
Chemical Engineering Thermodynamics*, 7th ed., Ch. 13, for the
Gibbs-energy-minimization approach to equilibrium composition used by
:func:`solve_equilibrium_composition` below (the standard alternative to
solving ``Q = K`` directly, and the one that generalizes cleanly to
several simultaneous reaction equilibria).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from chemistrykit.constants import STANDARD_PRESSURE, R
from chemistrykit.thermo.utils.regression import linear_fit

__all__ = [
    "reaction_quotient",
    "kp_from_kc",
    "kc_from_kp",
    "van_t_hoff_equilibrium_constant",
    "VantHoffFit",
    "fit_van_t_hoff",
    "EquilibriumComposition",
    "gibbs_energy_of_mixture",
    "solve_equilibrium_composition",
]


def reaction_quotient(activities, stoich_coeffs):
    r"""The reaction quotient :math:`Q = \prod_i a_i^{\nu_i}`.

    `stoich_coeffs` are *signed net* stoichiometric coefficients (positive
    for products, negative for reactants), so a reactant with coefficient
    :math:`-\nu` contributes :math:`a^{-\nu} = 1/a^{\nu}` to the product,
    matching the usual "products over reactants" definition (Atkins & de
    Paula, *Physical Chemistry*, 11th ed., Ch. 6.1). `activities` may be
    partial pressures relative to the standard pressure (dimensionless,
    for :math:`Q_p`) or concentrations relative to a standard
    concentration (for :math:`Q_c`); which one determines whether the
    result should be compared against :math:`K_p` or :math:`K_c`.

    Parameters
    ----------
    activities : array-like of float
        Activity (or activity-like relative pressure/concentration) of
        each species.
    stoich_coeffs : array-like of float
        Signed net stoichiometric coefficient of each species, same
        order as `activities`.

    Returns
    -------
    float

    Examples
    --------
    For :math:`N_2O_4 \rightleftharpoons 2NO_2`, at the equilibrium mole
    fractions corresponding to :math:`Q = K = 4` (i.e. extent
    :math:`\xi = \sqrt{K/(4+K)}`, as in
    :func:`solve_equilibrium_composition`'s worked example below):

    >>> import numpy as np
    >>> xi = np.sqrt(4.0 / 8.0)
    >>> x_N2O4, x_NO2 = (1.0 - xi) / (1.0 + xi), 2.0 * xi / (1.0 + xi)
    >>> round(float(reaction_quotient([x_N2O4, x_NO2], [-1.0, 2.0])), 3)
    4.0
    """
    activities = np.asarray(activities, dtype=np.float64)
    stoich_coeffs = np.asarray(stoich_coeffs, dtype=np.float64)
    return float(np.prod(activities**stoich_coeffs))


def kp_from_kc(Kc: float, delta_n: float, T: float, R_gas: float = R) -> float:
    r"""Convert :math:`K_c` to :math:`K_p` via :math:`K_p = K_c(RT)^{\Delta n}`.

    Valid for reactions among ideal gases, where :math:`\Delta n` is the
    change in moles of gas (products minus reactants) (Atkins & de Paula,
    *Physical Chemistry*, 11th ed., Ch. 6.2). `Kc` must be expressed in
    concentration units consistent with `R_gas` (e.g. mol/m^3 with the SI
    `R_gas`, or mol/L with ``R_gas = 0.0831446 L bar / (mol K)``).

    Parameters
    ----------
    Kc : float
        Equilibrium constant in terms of concentration.
    delta_n : float
        Change in moles of gas, products minus reactants.
    T : float
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, consistent with the concentration units of `Kc`.

    Returns
    -------
    float

    Examples
    --------
    With ``delta_n = 0`` (no change in moles of gas), Kp = Kc:

    >>> round(kp_from_kc(Kc=2.5, delta_n=0.0, T=298.15), 6)
    2.5
    """
    return Kc * (R_gas * T) ** delta_n


def kc_from_kp(Kp: float, delta_n: float, T: float, R_gas: float = R) -> float:
    r"""Convert :math:`K_p` to :math:`K_c` via :math:`K_c = K_p/(RT)^{\Delta n}` (inverse of :func:`kp_from_kc`).

    Parameters
    ----------
    Kp : float
        Equilibrium constant in terms of pressure.
    delta_n : float
        Change in moles of gas, products minus reactants.
    T : float
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant.

    Returns
    -------
    float

    Examples
    --------
    Round-trips with :func:`kp_from_kc`:

    >>> Kp = kp_from_kc(Kc=1.8, delta_n=1.0, T=350.0)
    >>> round(kc_from_kp(Kp, delta_n=1.0, T=350.0), 6)
    1.8
    """
    return Kp / (R_gas * T) ** delta_n


def van_t_hoff_equilibrium_constant(T, T_ref: float, K_ref: float, delta_h: float, R_gas: float = R):
    r"""The (integrated) van't Hoff equation :math:`\ln(K/K_{ref}) = -\Delta H^\circ/R\,(1/T - 1/T_{ref})`.

    Assumes the standard reaction enthalpy :math:`\Delta H^\circ` is
    constant over the temperature range of interest (Atkins & de Paula,
    *Physical Chemistry*, 11th ed., Ch. 6.4) -- the same approximation, and
    the same functional form, as :class:`chemistrykit.thermo.systems.phase_equilibria.ClausiusClapeyron`
    (a coincidence of both following from :math:`d\ln(\cdot)/dT \propto
    1/T^2`, one for a phase boundary, the other for an equilibrium
    constant).

    Parameters
    ----------
    T : float or array-like of float
        Absolute temperature(s) at which to evaluate `K`, in K.
    T_ref : float
        Reference temperature at which `K_ref` is known, in K.
    K_ref : float
        Equilibrium constant at `T_ref`.
    delta_h : float
        Standard reaction enthalpy, in J/mol (positive for endothermic).
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.

    Returns
    -------
    float or ndarray

    Examples
    --------
    For an endothermic reaction, `K` increases with temperature (Le
    Chatelier's principle):

    >>> import numpy as np
    >>> T = np.array([280.0, 300.0, 320.0])
    >>> K = van_t_hoff_equilibrium_constant(T, T_ref=298.15, K_ref=1.0, delta_h=50e3)
    >>> bool(np.all(np.diff(K) > 0))
    True
    """
    T = np.asarray(T, dtype=np.float64)
    return K_ref * np.exp(-delta_h / R_gas * (1.0 / T - 1.0 / T_ref))


@dataclass
class VantHoffFit:
    r"""Result of fitting :math:`\ln K` vs. :math:`1/T` data to the van't Hoff equation.

    Since :math:`\ln K = -\Delta H^\circ/(RT) + \Delta S^\circ/R` (Atkins &
    de Paula, *Physical Chemistry*, 11th ed., Ch. 6.4, eq. 6.16), the
    slope of a plot of :math:`\ln K` against :math:`1/T` gives
    :math:`-\Delta H^\circ/R` and the intercept gives
    :math:`\Delta S^\circ/R`.
    """

    delta_h: float
    """float: Fitted standard reaction enthalpy, in J/mol."""

    delta_s: float
    """float: Fitted standard reaction entropy, in J mol^-1 K^-1."""

    r_squared: float
    """float: Coefficient of determination of the linear (ln K vs 1/T) fit."""

    R_gas: float = R
    """float: Gas constant used in the fit, in J mol^-1 K^-1."""

    def predict(self, T):
        """Evaluate the fitted van't Hoff equation at temperature(s) `T`.

        Parameters
        ----------
        T : float or array-like of float
            Absolute temperature(s), in K.

        Returns
        -------
        float or ndarray
        """
        T = np.asarray(T, dtype=np.float64)
        return np.exp(-self.delta_h / (self.R_gas * T) + self.delta_s / self.R_gas)


def fit_van_t_hoff(T, K, R_gas: float = R) -> VantHoffFit:
    r"""Fit equilibrium-constant vs. temperature data to the van't Hoff equation.

    Linearizes :math:`\ln K = -\Delta H^\circ/R \cdot (1/T) + \Delta
    S^\circ/R` and fits by ordinary least squares -- the "van't Hoff plot"
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 6.4),
    structurally identical to
    :func:`chemistrykit.kinetics.systems.arrhenius.fit_arrhenius`'s
    Arrhenius plot.

    Parameters
    ----------
    T : array-like of float
        Absolute temperatures, in K (at least 2 distinct values).
    K : array-like of float
        Equilibrium constants measured at each temperature in `T`.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.

    Returns
    -------
    VantHoffFit

    Examples
    --------
    Generate exact data from a known (ΔH°, ΔS°) and recover ΔH°:

    >>> import numpy as np
    >>> T = np.array([280.0, 300.0, 320.0, 340.0, 360.0])
    >>> K_ref, delta_h = 1.0, 45_000.0
    >>> K = van_t_hoff_equilibrium_constant(T, T_ref=300.0, K_ref=K_ref, delta_h=delta_h)
    >>> fit = fit_van_t_hoff(T, K)
    >>> round(fit.delta_h, 2)
    45000.0
    >>> round(fit.r_squared, 6)
    1.0
    """
    T = np.asarray(T, dtype=np.float64)
    K = np.asarray(K, dtype=np.float64)
    fit = linear_fit(1.0 / T, np.log(K))
    delta_h = -fit.slope * R_gas
    delta_s = fit.intercept * R_gas
    return VantHoffFit(delta_h=float(delta_h), delta_s=float(delta_s), r_squared=fit.r_squared, R_gas=R_gas)


@dataclass
class EquilibriumComposition:
    """Result of a :func:`solve_equilibrium_composition` call."""

    species: tuple
    """tuple of str: Species names, in the order used throughout."""

    n: np.ndarray
    """ndarray, shape (n_species,): Equilibrium mole amounts."""

    x: np.ndarray
    """ndarray, shape (n_species,): Equilibrium mole fractions."""

    extents: np.ndarray
    """ndarray, shape (n_reactions,): Equilibrium extent(s) of reaction :math:`\\xi_j`."""

    converged: bool
    """bool: Whether the underlying optimizer reported success."""

    def moles(self, name: str) -> float:
        """Return the equilibrium mole amount of a single named species.

        Parameters
        ----------
        name : str

        Returns
        -------
        float
        """
        return float(self.n[list(self.species).index(name)])


def gibbs_energy_of_mixture(n, gibbs_formation, T: float, P: float = STANDARD_PRESSURE, P_standard: float = STANDARD_PRESSURE, R_gas: float = R) -> float:
    r"""Total Gibbs energy of an ideal-gas mixture, :math:`G = \sum_i n_i \mu_i`.

    Using the ideal-gas chemical potential :math:`\mu_i = \Delta G_{f,i}^\circ
    + RT\ln(x_i P/P^\circ)` (Smith, Van Ness & Abbott, *Introduction to
    Chemical Engineering Thermodynamics*, 7th ed., Ch. 13.2), where `x_i`
    is species `i`'s mole fraction. This is the objective function
    minimized by :func:`solve_equilibrium_composition`.

    Parameters
    ----------
    n : array-like of float
        Mole amount of each species (must sum to a positive total).
    gibbs_formation : array-like of float
        Standard Gibbs energy of formation of each species, in J/mol,
        same order as `n`.
    T : float
        Absolute temperature, in K.
    P : float, default :data:`chemistrykit.constants.STANDARD_PRESSURE`
        Total pressure, in Pa.
    P_standard : float, default :data:`chemistrykit.constants.STANDARD_PRESSURE`
        Standard-state pressure, in Pa.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.

    Returns
    -------
    float
        Total Gibbs energy, in J.
    """
    n = np.asarray(n, dtype=np.float64)
    gibbs_formation = np.asarray(gibbs_formation, dtype=np.float64)
    total = n.sum()
    x = np.clip(n, 1e-300, None) / total
    mu = gibbs_formation + R_gas * T * np.log(x * P / P_standard)
    return float(np.sum(n * mu))


def solve_equilibrium_composition(
    species,
    stoich_matrix,
    n0,
    gibbs_formation,
    T: float,
    P: float = STANDARD_PRESSURE,
    R_gas: float = R,
) -> EquilibriumComposition:
    r"""Solve for the equilibrium composition of a reacting ideal-gas mixture.

    Rather than solving :math:`Q(\xi) = K` algebraically (which gets
    unwieldy for several simultaneous reactions), this parameterizes the
    feasible mole amounts by the extent(s) of reaction :math:`\xi_j`,
    :math:`n_i(\xi) = n_{i,0} + \sum_j \nu_{ij}\xi_j` -- which
    automatically satisfies atomic mass balance for any :math:`\xi` -- and
    numerically minimizes the total Gibbs energy
    :func:`gibbs_energy_of_mixture` over the feasible region
    :math:`n_i(\xi) \geq 0`, via SLSQP (Smith, Van Ness & Abbott,
    *Introduction to Chemical Engineering Thermodynamics*, 7th ed., Ch.
    13.7). At the unconstrained minimum, :math:`dG/d\xi_j = \sum_i
    \nu_{ij}\mu_i = \Delta G_{rxn,j} = 0`, which is exactly the classical
    :math:`Q_j = K_j` equilibrium condition for every reaction `j` --
    Gibbs-energy minimization and root-finding on the equilibrium
    constant are the same equilibrium condition, just reached by
    different numerical routes.

    Parameters
    ----------
    species : sequence of str
        Ordered species names.
    stoich_matrix : array-like, shape (n_species, n_reactions)
        Net stoichiometric coefficient of each species in each reaction.
    n0 : array-like, shape (n_species,)
        Initial (pre-reaction) mole amounts.
    gibbs_formation : array-like, shape (n_species,)
        Standard Gibbs energy of formation of each species, in J/mol.
    T : float
        Absolute temperature, in K.
    P : float, default :data:`chemistrykit.constants.STANDARD_PRESSURE`
        Total pressure, in Pa.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.

    Returns
    -------
    EquilibriumComposition

    Examples
    --------
    :math:`N_2O_4 \rightleftharpoons 2NO_2`, choosing standard Gibbs
    energies of formation so that :math:`K = \exp(-\Delta G^\circ_{rxn}/RT)
    = 4` exactly at 298.15 K (with :math:`\Delta G_f^\circ(N_2O_4) = 0` as
    the reference), starting from 1 mol of pure :math:`N_2O_4`. Since
    :math:`Q(\xi) = x_{NO_2}^2/x_{N_2O_4} = 4\xi^2/(1-\xi^2)` at :math:`P
    = P^\circ`, setting :math:`Q = K` gives the closed-form extent
    :math:`\xi = \sqrt{K/(4+K)} = \sqrt{4/8} \approx 0.7071`, which the
    Gibbs-minimization solver should reproduce:

    >>> import numpy as np
    >>> from chemistrykit.constants import R
    >>> T = 298.15
    >>> K_target = 4.0
    >>> delta_g_rxn = -R * T * np.log(K_target)
    >>> gibbs_formation = [0.0, delta_g_rxn / 2.0]  # [N2O4, NO2]
    >>> result = solve_equilibrium_composition(
    ...     species=("N2O4", "NO2"),
    ...     stoich_matrix=[[-1.0], [2.0]],
    ...     n0=[1.0, 0.0],
    ...     gibbs_formation=gibbs_formation,
    ...     T=T,
    ... )
    >>> round(float(result.extents[0]), 4)
    0.7071
    >>> result.converged
    True
    """
    stoich = np.asarray(stoich_matrix, dtype=np.float64)
    n0 = np.asarray(n0, dtype=np.float64)
    gibbs_formation = np.asarray(gibbs_formation, dtype=np.float64)
    n_reactions = stoich.shape[1]

    def n_of_xi(xi):
        return n0 + stoich @ xi

    def objective(xi):
        return gibbs_energy_of_mixture(n_of_xi(xi), gibbs_formation, T, P=P, R_gas=R_gas)

    constraints = [{"type": "ineq", "fun": lambda xi: n_of_xi(xi) - 1e-12}]
    res = minimize(objective, np.zeros(n_reactions), method="SLSQP", constraints=constraints)
    n_eq = n_of_xi(res.x)
    x_eq = n_eq / n_eq.sum()
    return EquilibriumComposition(species=tuple(species), n=n_eq, x=x_eq, extents=res.x, converged=bool(res.success))
