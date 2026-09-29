"""Regenerate the README's per-subpackage teaser figures.

    python docs/make_readme_subpackage_figures.py                 # all subpackages
    python docs/make_readme_subpackage_figures.py kinetics thermo  # just these

Writes docs/source/_static/images/readme_<subpackage>.png, three panels each,
at the same size as readme_hero.png (see make_readme_figure.py). README.md
embeds them by their raw.githubusercontent.com URLs so they also render on PyPI.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path(__file__).parent / "source" / "_static" / "images"


def _replace_with_3d(ax):
    """Swap a 2D axes for a 3D one in the same grid slot."""
    fig, spec = ax.figure, ax.get_subplotspec()
    ax.remove()
    return fig.add_subplot(spec, projection="3d")


def kinetics(axes):
    from chemistrykit.kinetics.systems.networks import consecutive_analytic
    from chemistrykit.kinetics.systems.oscillators import Brusselator, Oregonator
    from chemistrykit.kinetics.systems.stochastic import gillespie_ssa
    from chemistrykit.kinetics.visualizers.kinetics_plots import plot_phase_portrait

    ax1, ax2, ax3 = axes
    for x0, y0, color in ((1.0, 1.0, "steelblue"), (0.3, 4.5, "darkorange")):
        plot_phase_portrait(Brusselator(X0=x0, Y0=y0, A=1.0, B=3.0).integrate((0.0, 40.0), dt=1e-2, method="rk4"), "X", "Y", ax=ax1, color=color)
    ax1.plot(*Brusselator(X0=1.0, Y0=1.0, A=1.0, B=3.0).fixed_point(), "k*", ms=12, label="unstable fixed point")
    ax1.set_title("Brusselator: two starts, one limit cycle")
    ax1.legend(fontsize=8)

    result = Oregonator(f=1.0).integrate((0.0, 40.0), method="dopri5", rtol=1e-7, atol=1e-10, max_steps=2_000_000)
    for series, label in zip(result.y.T, ("[HBrO$_2$]", "[Br$^-$]", "[Ce(IV)]"), strict=True):
        ax2.semilogy(result.t, series, label=label)
    ax2.set_ylim(top=1e4)
    ax2.set_xlabel("t (scaled)")
    ax2.set_title("Oregonator: Belousov-Zhabotinsky oscillations")
    ax2.legend(fontsize=8, ncol=3, loc="upper center")

    n, t_grid = 20, np.linspace(0.0, 15.0, 301)
    rng = np.random.default_rng(1976)
    for _ in range(8):
        traj = gillespie_ssa([[-1, 0], [1, -1], [0, 1]], [1.0, 0.3], [[1, 0], [0, 1], [0, 0]], [n, 0, 0], t_max=15.0, species=("A", "B", "C"), seed=rng)
        ax3.step(traj.t, traj.count("B") / n, where="post", color="steelblue", alpha=0.5, lw=0.8)
    ax3.plot(t_grid, consecutive_analytic(1.0, 1.0, 0.3, t_grid)[1], "k--", label="rate equations")
    ax3.set_xlabel("t")
    ax3.set_ylabel("[B] / [A]$_0$")
    ax3.set_title(f"Gillespie: A → B → C with {n} molecules")
    ax3.legend(fontsize=8)


def thermo(axes):
    from chemistrykit.thermo import AntoineEquation, WilsonSolution
    from chemistrykit.thermo.systems.equations_of_state import VanDerWaals
    from chemistrykit.thermo.systems.phase_equilibria import ClausiusClapeyron

    ax1, ax2, ax3 = axes
    Tc, Pc = 304.13, 7.3773e6  # CO2
    vdw = VanDerWaals.from_critical_constants(Tc, Pc)
    Vc = 3.0 * vdw.b
    Vm = np.linspace(1.4 * vdw.b, 8.0 * Vc, 600)
    for Tr in (0.85, 0.9, 0.95, 1.0, 1.1, 1.2):
        ax1.plot(Vm / Vc, vdw.pressure(Vm, Tr * Tc) / Pc, "k-" if Tr == 1.0 else "-", label=rf"$T/T_c$ = {Tr}")
    ax1.plot([1.0], [1.0], "ro")
    ax1.set_xlim(0.4, 8.0)
    ax1.set_ylim(0.0, 2.0)
    ax1.set_xlabel(r"$V_m / V_c$")
    ax1.set_ylabel(r"$P / P_c$")
    ax1.set_title("van der Waals isotherms of CO$_2$")
    ax1.legend(fontsize=7)

    T = 343.15
    mix = WilsonSolution(
        Lambda12=0.1782,
        Lambda21=0.8703,
        P1_star=float(AntoineEquation(5.37229, 1670.409, -40.191).pressure(T)) * 100.0,
        P2_star=float(AntoineEquation(5.08354, 1663.125, -45.622).pressure(T)) * 100.0,
    )
    x = np.linspace(0.0, 1.0, 400)
    ax2.plot(x, mix.total_pressure(x), color="darkorange", label="bubble curve")
    ax2.plot(mix.vapor_composition(x), mix.total_pressure(x), color="steelblue", label="dew curve")
    x_az, P_az = mix.azeotrope()
    ax2.plot([x_az], [P_az], "ko", label=f"azeotrope, x = {x_az:.2f}")
    ax2.set_xlabel("mole fraction of ethanol")
    ax2.set_ylabel("P (kPa)")
    ax2.set_title("Ethanol-water azeotrope (Wilson model)")
    ax2.legend(fontsize=8)

    T_tp, P_tp = 273.16, 611.657  # water triple point
    T_vap, T_sub = np.linspace(T_tp, 300.0, 100), np.linspace(250.0, T_tp, 100)
    P_melt = np.logspace(np.log10(P_tp), 7, 100)
    ax3.plot(T_vap, ClausiusClapeyron(delta_h_vap=45050.0, T_ref=T_tp, P_ref=P_tp).pressure(T_vap), label="liquid-vapor")
    ax3.plot(T_sub, ClausiusClapeyron(delta_h_vap=51060.0, T_ref=T_tp, P_ref=P_tp).pressure(T_sub), label="solid-vapor")
    ax3.plot(T_tp + (P_melt - P_tp) * T_tp * -1.63e-6 / 6010.0, P_melt, label="solid-liquid")
    ax3.plot([T_tp], [P_tp], "ko", label="triple point")
    for tx, py, text in ((263.0, 1e4, "ice"), (288.0, 1e5, "liquid"), (288.0, 200.0, "vapor")):
        ax3.text(tx, py, text, ha="center")
    ax3.set_yscale("log")
    ax3.set_xlim(250.0, 300.0)
    ax3.set_ylim(50.0, 1e7)
    ax3.set_xlabel("T (K)")
    ax3.set_ylabel("P (Pa)")
    ax3.set_title("Water's phase diagram (Clausius-Clapeyron)")
    ax3.legend(fontsize=8, loc="upper left")


def solutions(axes):
    from chemistrykit.solutions.systems.acid_base import polyprotic_fractions
    from chemistrykit.solutions.systems.activity import activity_coefficient_davies, activity_coefficient_debye_huckel_limiting
    from chemistrykit.solutions.systems.titration import StrongAcidStrongBaseTitration, WeakAcidStrongBaseTitration

    ax1, ax2, ax3 = axes
    Vb = np.linspace(1e-6, 0.09, 2000)
    for titration, label in (
        (StrongAcidStrongBaseTitration(Ca=0.1, Va=0.05, Cb=0.1), "strong acid (HCl)"),
        (WeakAcidStrongBaseTitration(Ca=0.1, Va=0.05, Ka=1.8e-5, Cb=0.1), "weak acid (acetic)"),
    ):
        curve = titration.curve(Vb)
        ax1.plot(curve.Vb * 1000.0, curve.pH, label=label)
    ax1.axvline(50.0, color="gray", ls=":", lw=0.8)
    ax1.set_xlabel("NaOH added (mL)")
    ax1.set_ylabel("pH")
    ax1.set_title("Titration curves with 0.1 M NaOH")
    ax1.legend(fontsize=8)

    pH = np.linspace(0, 14, 500)
    pKas = (2.15, 7.20, 12.35)
    for row, label in zip(polyprotic_fractions(pH, [10.0**-pK for pK in pKas]), ("H$_3$PO$_4$", "H$_2$PO$_4^-$", "HPO$_4^{2-}$", "PO$_4^{3-}$"), strict=True):
        ax2.plot(pH, row, label=label)
    ax2.set_xlabel("pH")
    ax2.set_ylabel(r"fraction $\alpha_j$")
    ax2.set_title("Bjerrum speciation of phosphoric acid")
    ax2.legend(fontsize=8)

    ionic = np.linspace(1e-4, 0.8, 400)
    for z, color in ((1, "steelblue"), (2, "crimson")):
        ax3.plot(ionic, [activity_coefficient_debye_huckel_limiting(z, i) for i in ionic], ":", color=color, label=f"Debye-Hückel limiting, z = {z}")
        ax3.plot(ionic, [activity_coefficient_davies(z, i) for i in ionic], color=color, label=f"Davies, z = {z}")
    ax3.set_ylim(0.0, 1.05)
    ax3.set_xlabel("ionic strength I (mol/L)")
    ax3.set_ylabel(r"activity coefficient $\gamma$")
    ax3.set_title("Activity coefficients: Debye-Hückel vs. Davies")
    ax3.legend(fontsize=7)


def md(axes):
    from chemistrykit.md.systems.lj_fluid import LJFluid
    from chemistrykit.md.systems.thermostats import VelocityRescalingThermostat
    from chemistrykit.md.systems.transport import einstein_diffusion_coefficient, mean_squared_displacement
    from chemistrykit.md.visualizers.md_plots import plot_speed_distribution

    ax1, ax2, ax3 = axes
    T_liquid, rho_liquid = 94.4 / 120.0, 0.84  # Rahman's liquid argon, reduced units
    fluid = LJFluid.from_lattice(n_per_side=7, density=rho_liquid, temperature=T_liquid, cutoff=2.5, rng=0)
    fluid.run(dt=0.005, n_steps=2000, thermostat=VelocityRescalingThermostat(T_liquid, interval=10), sample_every=2000)
    result = fluid.run(dt=0.005, n_steps=1000, sample_every=100)
    r, _ = fluid.radial_distribution_function(n_bins=150)
    g = np.mean([fluid.radial_distribution_function(n_bins=150, positions=p)[1] for p in result.positions], axis=0)
    ax1.plot(r * 3.4, g, color="steelblue")
    ax1.axhline(1.0, color="gray", ls=":", lw=0.8)
    ax1.set_xlabel(r"r ($\mathrm{\AA}$)")
    ax1.set_ylabel("g(r)")
    ax1.set_title("Liquid argon: radial distribution function")

    gas = LJFluid.from_lattice(n_per_side=7, density=0.5, temperature=1.5, cutoff=2.5, rng=1)
    gas.run(dt=0.002, n_steps=1000, sample_every=1000)
    plot_speed_distribution(np.linalg.norm(gas.velocities, axis=1), mass=1.0, temperature=gas.temperature(), ax=ax2, n_bins=30)
    ax2.set_title("MD speeds vs. Maxwell-Boltzmann")

    dt, sample_every, max_lag = 0.005, 4, 500
    fluid = LJFluid.from_lattice(n_per_side=6, density=0.7, temperature=1.0, cutoff=2.5, rng=0)
    fluid.run(dt=dt, n_steps=2000, thermostat=VelocityRescalingThermostat(1.0, interval=10), sample_every=2000)
    result = fluid.run(dt=dt, n_steps=4000, sample_every=sample_every)
    msd = mean_squared_displacement(result.positions, box_length=result.box_length, max_lag=max_lag)
    lag_t = np.arange(max_lag + 1) * dt * sample_every
    D = einstein_diffusion_coefficient(lag_t, msd, fit_from=2.0)
    ax3.loglog(lag_t[1:], msd[1:], color="steelblue", label="simulated MSD")
    ax3.loglog(lag_t[1:40], 3.0 * lag_t[1:40] ** 2, color="gray", ls=":", label=r"ballistic $\langle v^2\rangle t^2$")
    ax3.loglog(lag_t[40:], 6.0 * D * lag_t[40:], "k--", label=f"diffusive 6Dt, D = {D:.3f}")
    ax3.set_xlabel("t (reduced units)")
    ax3.set_title("Einstein's diffusion law: ballistic, then diffusive")
    ax3.legend(fontsize=8)


def statmech(axes):
    from scipy.optimize import brentq

    from chemistrykit.constants import K_B, H, R
    from chemistrykit.statmech import DebyeSolid, Ising2DOnsager, MaxwellBoltzmannSpeedDistribution, VibrationalPartitionFunctionHarmonic

    ax1, ax2, ax3 = axes
    v = np.linspace(0.0, 2500.0, 400)
    for T, color in ((100.0, "steelblue"), (300.0, "seagreen"), (1000.0, "crimson")):
        dist = MaxwellBoltzmannSpeedDistribution(mass=6.63e-26, temperature=T)
        ax1.plot(v, dist.pdf(v), color=color, label=f"T = {T:.0f} K")
        ax1.axvline(dist.most_probable_speed(), color=color, ls=":", lw=0.8)
    ax1.set_xlabel("speed (m/s)")
    ax1.set_title("Maxwell-Boltzmann speeds of argon")
    ax1.legend(fontsize=8)

    J = 100.0 * K_B
    model = Ising2DOnsager(coupling=J)
    Tc = model.critical_temperature
    T = np.linspace(20.0, 1.6 * Tc, 500)
    T = T[np.abs(T - Tc) > 1e-6 * Tc]

    def mean_field(t):
        return 0.0 if t >= 4 * J / K_B else brentq(lambda m: m - np.tanh(4 * J * m / (K_B * t)), 1e-9, 1.0)

    ax2.plot(T / Tc, model.spontaneous_magnetization(T), label="Onsager-Yang (exact)")
    ax2.plot(T / Tc, [mean_field(t) for t in T], "--", color="gray", label="mean field")
    ax2.set_xlabel(r"$T / T_c$")
    ax2.set_ylabel("spontaneous magnetization")
    ax2.set_title("2D Ising: Onsager's exact solution")
    ax2.legend(fontsize=8)

    theta_D = 343.0  # copper
    T = np.linspace(1.0, 600.0, 400)
    einstein = VibrationalPartitionFunctionHarmonic(frequency=np.sqrt(3.0 / 5.0) * theta_D * K_B / H)
    ax3.plot(T / theta_D, DebyeSolid(debye_temperature=theta_D).heat_capacity_v(T) / (3 * R), label="Debye")
    ax3.plot(T / theta_D, 3.0 * einstein.heat_capacity_v(T) / (3 * R), "--", label="Einstein")
    ax3.axhline(1.0, color="gray", ls=":", lw=0.8, label="Dulong-Petit 3R")
    ax3.set_xlabel(r"$T / \Theta_D$")
    ax3.set_ylabel(r"$C_V / 3R$")
    ax3.set_title("Heat capacity of copper: Debye vs. Einstein")
    ax3.legend(fontsize=8)


def quantum(axes):
    from chemistrykit.quantum.systems.hartree_fock import H2PlusVariational
    from chemistrykit.quantum.systems.huckel import HuckelSystem
    from chemistrykit.quantum.systems.hydrogenlike import HydrogenLikeAtom
    from chemistrykit.quantum.visualizers.quantum_plots import plot_radial_distribution

    ax1, ax2, ax3 = axes
    atom = HydrogenLikeAtom(Z=1)
    for n, l, style in ((1, 0, "-"), (2, 0, "--"), (2, 1, "-."), (3, 2, ":")):
        plot_radial_distribution(atom, n, l, ax=ax1, r_max_bohr_radii=25.0, linestyle=style, label=f"n={n}, l={l}")
    ax1.set_title("Hydrogen radial distribution functions")
    ax1.legend(fontsize=8)

    bond = 106.0e-12
    h2plus = H2PlusVariational(bond_length=bond)
    alpha = h2plus.optimize_exponent(alpha_guess=1.0 / 5.29177e-11**2).optimized_alpha
    result = h2plus.solve(alpha)
    z = np.linspace(-3.0 * bond, 3.0 * bond, 801)
    chi_a, chi_b = ((2.0 * alpha / np.pi) ** 0.75 * np.exp(-alpha * (z - c) ** 2) for c in (-bond / 2.0, bond / 2.0))
    for k, label, color in ((0, "bonding", "steelblue"), (1, "antibonding", "crimson")):
        c_a, c_b = result.orbital_coefficients(k)
        ax2.plot(z * 1e12, (c_a * chi_a + c_b * chi_b) ** 2, color=color, label=label)
    for c in (-bond / 2.0, bond / 2.0):
        ax2.axvline(c * 1e12, color="black", ls="--", lw=0.6)
    ax2.set_xlabel("position along the H$_2^+$ axis (pm)")
    ax2.set_ylabel(r"$|\psi|^2$")
    ax2.set_yticks([])
    ax2.set_title("H$_2^+$ molecular orbitals (LCAO variational)")
    ax2.legend(fontsize=8)

    theta = np.linspace(0.0, 2.0 * np.pi, 200)
    for n, x0 in ((4, 0.0), (5, 5.0), (6, 10.0)):
        energies = np.sort(HuckelSystem.cyclic_polyene(n).solve().energies)
        angles = -np.pi / 2.0 + 2.0 * np.pi * np.arange(n) / n
        ax3.plot(x0 + 2.0 * np.cos(theta), 2.0 * np.sin(theta), color="lightgray")
        ax3.fill(x0 + 2.0 * np.cos(angles), 2.0 * np.sin(angles), facecolor="none", edgecolor="steelblue")
        ax3.plot(x0 + 2.0 * np.cos(angles), 2.0 * np.sin(angles), "o", color="crimson")
        ax3.hlines(energies, x0 - 2.2, x0 + 2.2, color="black", lw=0.5, ls=":")
        ax3.text(x0, -2.8, f"C$_{n}$H$_{n}$", ha="center")
    ax3.set_aspect("equal", adjustable="datalim")
    ax3.set_xticks([])
    ax3.set_ylabel(r"$(E - \alpha) / |\beta|$")
    ax3.set_title("Hückel π energies on Frost circles")


def spectro(axes):
    import scipy.constants as sc

    from chemistrykit.quantum.systems.rigid_rotor import RigidRotor
    from chemistrykit.spectro.systems.electronic import franck_condon_spectrum
    from chemistrykit.spectro.systems.nmr import first_order_multiplet
    from chemistrykit.spectro.systems.rotational import rotational_spectrum
    from chemistrykit.spectro.visualizers.spectro_plots import plot_broadened_spectrum, plot_stick_spectrum

    ax1, ax2, ax3 = axes
    ppm = np.linspace(0.8, 4.1, 6000)
    for shift, neighbors, label in ((1.2, 2, "CH$_3$ triplet"), (3.7, 3, "CH$_2$ quartet")):
        multiplet = first_order_multiplet(chemical_shift_ppm=shift, j_coupling_hz=7.0, n_neighbors=neighbors, spectrometer_frequency_mhz=60.0)
        plot_broadened_spectrum(multiplet, ppm, ax=ax1, shape="lorentzian", fwhm=0.005, label=label)
    ax1.invert_xaxis()
    ax1.set_xlabel("chemical shift (ppm)")
    ax1.set_title("Ethanol $^1$H NMR at 60 MHz: the n+1 rule")
    ax1.legend(fontsize=8)

    spectrum = franck_condon_spectrum(origin_wavenumber=22000.0, vibrational_wavenumber=1400.0, S=1.8, v_max=12)
    x = np.linspace(21000.0, 32000.0, 2000)
    envelope = sum(h * np.exp(-0.5 * ((x - p) / 350.0) ** 2) for p, h in zip(spectrum.positions, spectrum.intensities, strict=True))
    ax2.fill_between(x, envelope, color="indigo", alpha=0.2, label="broadened band")
    plot_stick_spectrum(spectrum, ax=ax2, color="darkviolet")
    ax2.set_xlim(21000.0, 32000.0)
    ax2.legend(fontsize=8)
    ax2.set_xlabel(r"wavenumber (cm$^{-1}$)")
    ax2.set_title("Franck-Condon vibronic progression")

    hcl = RigidRotor.from_diatomic(mass1=1.008 * sc.atomic_mass, mass2=34.97 * sc.atomic_mass, bond_length=127.5e-12)
    plot_stick_spectrum(rotational_spectrum(hcl, J_max=15, temperature=300.0), ax=ax3, color="steelblue")
    ax3.set_xlabel(r"wavenumber (cm$^{-1}$)")
    ax3.set_title("HCl rotational spectrum at 300 K")


def structure(axes):
    from chemistrykit.structure.systems.bonding import PAULING_C_C_CONSTANT, bond_order_from_length
    from chemistrykit.structure.systems.vsepr import build_vsepr_molecule
    from chemistrykit.structure.visualizers.structure_plots import plot_bond_order_correlation, plot_molecule_3d

    ax1, ax2, ax3 = axes
    ax1 = _replace_with_3d(ax1)
    plot_molecule_3d(build_vsepr_molecule(5, 1, bond_length=1.60, central_symbol="S", ligand_symbol="F"), ax=ax1)
    ax1.set_axis_off()
    ax1.set_title("VSEPR: SF$_4$'s seesaw geometry in 3D")

    for k in range(5):
        ax2.hlines(0.0, 0.1 + 0.16 * k, 0.24 + 0.16 * k, color="gray", lw=3)
    ax2.text(0.47, 0.06, "free ion", ha="center")
    for level, x0, count, color, label in ((0.6, 1.4, 2, "C0", r"$e_g$ ($+\frac{3}{5}\Delta_o$)"), (-0.4, 1.3, 3, "C3", r"$t_{2g}$ ($-\frac{2}{5}\Delta_o$)")):
        for k in range(count):
            ax2.hlines(level, x0 + 0.2 * k, x0 + 0.14 + 0.2 * k, color=color, lw=3)
        ax2.plot([0.9, x0], [0.0, level], color="lightgray", ls=":", lw=1)
        ax2.text(1.95, level, label, va="center")
    ax2.set_xlim(0, 2.6)
    ax2.set_ylim(-0.7, 0.9)
    ax2.set_xticks([])
    ax2.set_ylabel(r"energy / $\Delta_o$")
    ax2.set_title(r"Crystal-field splitting in an $O_h$ field")

    plot_bond_order_correlation(single_bond_length=1.54, bond_orders=np.linspace(0.8, 3.3, 200), c=PAULING_C_C_CONSTANT, ax=ax3, color="gray")
    for name, (n, length) in {"ethane": (1, 1.54), "ethylene": (2, 1.339), "acetylene": (3, 1.203)}.items():
        ax3.scatter([n], [length], zorder=3, label=name)
    ax3.scatter([bond_order_from_length(1.54, 1.397)], [1.397], marker="s", zorder=3, label="benzene (from length)")
    ax3.set_ylabel(r"C-C bond length ($\mathrm{\AA}$)")
    ax3.set_title("Pauling's bond-order/length correlation")
    ax3.legend(fontsize=8)


def electrochem(axes):
    from chemistrykit.electrochem.systems.butler_volmer import butler_volmer_current_density
    from chemistrykit.electrochem.systems.voltammetry import ilkovic_diffusion_current, polarographic_wave_current
    from chemistrykit.electrochem.visualizers.electrochem_plots import plot_tafel

    ax1, ax2, ax3 = axes
    i0, f = 1.0e-6, 96485.33212 / (8.314462618 * 298.15)
    eta = np.linspace(-0.2, 0.2, 400)
    ax1.plot(eta, i0 * np.exp(0.5 * f * eta) * 1e6, "--", label="anodic")
    ax1.plot(eta, -i0 * np.exp(-0.5 * f * eta) * 1e6, "--", label="cathodic")
    ax1.plot(eta, butler_volmer_current_density(i0, eta, alpha=0.5, n=1) * 1e6, "k", label="net")
    ax1.set_ylim(-50, 50)
    ax1.set_xlabel(r"overpotential $\eta$ (V)")
    ax1.set_ylabel(r"$i$ ($\mu$A cm$^{-2}$)")
    ax1.set_title("Butler-Volmer electrode kinetics")
    ax1.legend(fontsize=8)

    E = np.linspace(-0.3, -1.3, 400)
    total = sum(
        polarographic_wave_current(E, E_half, ilkovic_diffusion_current(2, D, 2.0, 4.0, C, average=True), n=2)
        for D, E_half, C in ((7.2e-6, -0.60, 0.5), (7.0e-6, -1.00, 1.0))
    )
    ax2.plot(E, total)
    ax2.annotate("Cd$^{2+}$", (-0.6, total[np.argmin(np.abs(E + 0.6))]), xytext=(-0.45, 3.0), arrowprops={"arrowstyle": "->"})
    ax2.annotate("Zn$^{2+}$", (-1.0, total[np.argmin(np.abs(E + 1.0))]), xytext=(-0.85, 8.0), arrowprops={"arrowstyle": "->"})
    ax2.invert_xaxis()
    ax2.set_xlabel("E (V vs. SCE)")
    ax2.set_ylabel(r"mean current ($\mu$A)")
    ax2.set_title("Heyrovský polarogram: one wave per ion")

    plot_tafel(i0=i0, eta_range=(-0.4, 0.4), alpha=0.5, n=1, ax=ax3)
    ax3.set_title("Tafel plot: log current linear in overpotential")


def photochem(axes):
    from chemistrykit.photochem import forster_efficiency, forster_radius
    from chemistrykit.photochem.systems.jablonski import jablonski_network
    from chemistrykit.photochem.systems.stern_volmer import dynamic_quenching_constant, fit_stern_volmer, stern_volmer_ratio
    from chemistrykit.photochem.visualizers.photochem_plots import plot_state_populations, plot_stern_volmer

    ax1, ax2, ax3 = axes
    plot_state_populations(jablonski_network(2.0, 1.0, 0.5, 0.3, 0.2, S1_0=1.0).integrate((0.0, 15.0), dt=1e-3, method="rk4"), ax=ax1)
    ax1.set_title("Jablonski kinetics: S$_1$ → T$_1$ → S$_0$")
    ax1.legend(fontsize=8)

    Q = np.array([0.0, 0.005, 0.010, 0.020, 0.040])
    ratio = stern_volmer_ratio(dynamic_quenching_constant(2.0e10, 5.0e-9), Q) * (1 + np.random.default_rng(0).normal(0.0, 0.01, Q.size))
    plot_stern_volmer(Q, ratio, fit=fit_stern_volmer(Q, ratio), ax=ax2)
    ax2.set_title("Stern-Volmer quenching")

    R0 = forster_radius(kappa2=2 / 3, n=1.4, quantum_yield_donor=0.5, overlap_J=2.0e15) / 10.0
    r = np.linspace(1.0, 12.0, 400)
    ax3.plot(r, forster_efficiency(r, R0))
    ax3.axvline(R0, color="k", ls=":", label=rf"$R_0$ = {R0:.1f} nm")
    ax3.axvspan(0.5 * R0, 1.5 * R0, color="0.92", label="useful ruler range")
    ax3.set_xlabel("donor-acceptor distance (nm)")
    ax3.set_ylabel("transfer efficiency")
    ax3.set_title("FRET: the inverse-sixth-power ruler")
    ax3.legend(fontsize=8)


def surface(axes):
    from chemistrykit.surface.systems.bet import BETIsotherm
    from chemistrykit.surface.systems.langmuir import LangmuirIsotherm, langmuir_coverage
    from chemistrykit.surface.systems.tpd import simulate_tpd

    ax1, ax2, ax3 = axes
    P = np.linspace(0.0, 8.0, 200)
    ax1.plot(P, BETIsotherm(Vm=5.0, C=80.0, P0=10.0).loading(P), label="BET (multilayer)")
    ax1.plot(P, LangmuirIsotherm(K=8.0, qmax=5.0).loading(P), "--", label="Langmuir (monolayer)")
    ax1.axhline(5.0, color="gray", ls=":", lw=0.8)
    ax1.set_xlabel("P")
    ax1.set_ylabel("amount adsorbed")
    ax1.set_title("Adsorption isotherms: BET vs. Langmuir")
    ax1.legend(fontsize=8)

    K = np.logspace(-3, 3, 4000)
    theta = langmuir_coverage(K, 1.0)
    ax2.semilogx(K, 4.0 * theta * (1.0 - theta))
    for name, k in {"M1": 0.005, "M2": 0.05, "M3": 0.4, "M4": 3.0, "M5": 30.0, "M6": 300.0}.items():
        th = langmuir_coverage(k, 1.0)
        ax2.plot(k, 4 * th * (1 - th), "o", color="crimson")
        ax2.annotate(name, (k, 4 * th * (1 - th)), textcoords="offset points", xytext=(4, 4))
    ax2.set_xlabel("adsorption constant K (binding strength)")
    ax2.set_ylabel("relative activity")
    ax2.set_title("Balandin's volcano: the Sabatier principle")

    for beta in (1.0, 10.0, 100.0):
        res = simulate_tpd(110.0e3, 1.0e13, beta, T_start=250.0, T_end=650.0)
        ax3.plot(res.T, res.desorption_rate / res.desorption_rate.max(), label=rf"$\beta$ = {beta:g} K/s")
    ax3.set_xlabel("T (K)")
    ax3.set_ylabel("desorption rate (normalized)")
    ax3.set_title("Temperature-programmed desorption")
    ax3.legend(fontsize=8)


def polymer(axes):
    from chemistrykit.polymer.systems.chain_statistics import freely_jointed_chain
    from chemistrykit.polymer.systems.flory_huggins import flory_huggins_critical_point, flory_huggins_spinodal_chi
    from chemistrykit.polymer.systems.molecular_weight_distribution import flory_schulz_number_fraction, flory_schulz_weight_fraction

    ax1, ax2, ax3 = axes
    walk = freely_jointed_chain(1000, 1.0, rng=7)[0]
    ax1.plot(walk[:, 0], walk[:, 1], lw=0.6)
    ax1.plot(*walk[0, :2], "go", label="start")
    ax1.plot(*walk[-1, :2], "rs", label="end")
    ax1.set_aspect("equal", adjustable="datalim")
    ax1.set_title("A 1000-segment freely jointed chain")
    ax1.legend(fontsize=8)

    x = np.arange(1, 150)
    for p, color in ((0.9, "steelblue"), (0.95, "crimson")):
        ax2.plot(x, flory_schulz_number_fraction(x, p), color=color, label=f"number fraction, p = {p}")
        ax2.plot(x, flory_schulz_weight_fraction(x, p), "--", color=color, label=f"weight fraction, p = {p}")
    ax2.set_xlabel("degree of polymerization")
    ax2.set_title("Flory-Schulz distribution (step growth)")
    ax2.legend(fontsize=7)

    phi = np.linspace(1e-3, 0.95, 800)
    for N in (10, 100, 1000):
        ax3.plot(phi, flory_huggins_spinodal_chi(phi, N), label=f"N = {N}")
        ax3.plot(*flory_huggins_critical_point(N), "ko", ms=4)
    ax3.axhline(0.5, color="gray", ls="--", lw=0.8, label=r"theta, $\chi$ = 1/2")
    ax3.set_ylim(0.3, 1.5)
    ax3.set_xlabel(r"polymer volume fraction $\phi$")
    ax3.set_ylabel(r"spinodal $\chi$")
    ax3.set_title("Flory-Huggins: spinodals and critical points")
    ax3.legend(fontsize=8)


def crystal(axes):
    from chemistrykit.crystal.systems.madelung import MADELUNG_CONSTANT_NACL_LITERATURE, madelung_constant_nacl
    from chemistrykit.crystal.systems.packing import BodyCenteredCubicPacking, FaceCenteredCubicPacking, HexagonalClosePacking, SimpleCubicPacking
    from chemistrykit.crystal.systems.xrd import powder_xrd_peaks
    from chemistrykit.crystal.visualizers.crystal_plots import plot_madelung_convergence, plot_packing_fractions

    ax1, ax2, ax3 = axes
    for peaks, label, color, offset in (
        (powder_xrd_peaks("FCC", a=361.5, wavelength=154.18, hkl_max=3), "Cu (FCC)", "C1", 0.0),
        (powder_xrd_peaks("BCC", a=286.65, wavelength=154.18, hkl_max=3), "Fe (BCC)", "C0", 1.2),
    ):
        intensities = np.array([p.relative_intensity * p.multiplicity for p in peaks], dtype=float)
        intensities /= intensities.max()
        for peak, intensity in zip(peaks, intensities, strict=True):
            ax1.vlines(peak.two_theta, offset, offset + intensity, color=color)
            ax1.text(peak.two_theta, offset + intensity + 0.02, "".join(map(str, peak.hkl)), ha="center", fontsize=6)
        ax1.text(22, offset + 0.8, label, color=color)
    ax1.set_yticks([])
    ax1.set_xlabel(r"$2\theta$ (degrees), Cu K$\alpha$")
    ax1.set_title("Powder XRD: systematic absences")

    shells = list(range(2, 16))
    plot_madelung_convergence(shells, [madelung_constant_nacl(n) for n in shells], literature_value=MADELUNG_CONSTANT_NACL_LITERATURE, ax=ax2)
    ax2.set_title("NaCl Madelung constant by Evjen's method")

    lattices = {"SC": SimpleCubicPacking(), "BCC": BodyCenteredCubicPacking(), "FCC": FaceCenteredCubicPacking(), "HCP": HexagonalClosePacking()}
    plot_packing_fractions(list(lattices), [lattice.packing_fraction() for lattice in lattices.values()], ax=ax3)
    ax3.axhline(np.pi / np.sqrt(18.0), color="k", ls="--", lw=1, label=r"Kepler bound $\pi/\sqrt{18}$")
    ax3.set_title("Hard-sphere packing fractions")
    ax3.legend(fontsize=8, loc="lower right")


def analytical(axes):
    from chemistrykit.analytical import EDTATitration, minimum_plate_height, optimum_flow_velocity, simulate_chromatogram, van_deemter_H

    ax1, ax2, ax3 = axes
    pigments = {"carotenes": (0.3, "orange"), "chlorophyll a": (1.5, "darkgreen"), "chlorophyll b": (2.4, "yellowgreen"), "xanthophylls": (3.6, "gold")}
    t_R = [20.0 * (1.0 + k) for k, _ in pigments.values()]
    t = np.linspace(0.0, 1.2 * max(t_R), 4000)
    trace = simulate_chromatogram(t, centers=t_R, N=400)
    ax1.plot(t, trace, color="black")
    for (name, (_k, color)), tr in zip(pigments.items(), t_R, strict=True):
        height = trace[np.argmin(np.abs(t - tr))]
        ax1.fill_between(t, trace, where=np.abs(t - tr) < 3.0 * tr / 20.0, color=color, alpha=0.6)
        ax1.text(tr, height + 0.01, name, ha="center", va="bottom", fontsize=8)
    ax1.set_ylim(0, 1.2 * trace.max())
    ax1.set_xlabel("time (min)")
    ax1.set_ylabel("detector signal")
    ax1.set_title("Tsvet's leaf pigments, separated")

    A, B, C = 1.5, 25.0, 0.05
    u = np.linspace(0.5, 60.0, 3000)
    ax2.plot(u, van_deemter_H(u, A, B, C), "k", lw=2, label="H = A + B/u + Cu")
    ax2.plot(u, np.full_like(u, A), "--", label="A (eddy diffusion)")
    ax2.plot(u, B / u, "--", label="B/u (longitudinal)")
    ax2.plot(u, C * u, "--", label="Cu (mass transfer)")
    ax2.plot([optimum_flow_velocity(B, C)], [minimum_plate_height(A, B, C)], "o", color="crimson")
    ax2.set_ylim(0, 12)
    ax2.set_xlabel("linear velocity u")
    ax2.set_ylabel("plate height H")
    ax2.set_title("The van Deemter equation")
    ax2.legend(fontsize=7)

    V = np.linspace(1e-6, 0.080, 4000)
    for logK in (6, 8, 10, 12):
        curve = EDTATitration(C_metal=0.01, V_metal=0.05, K_conditional=10.0**logK, C_edta=0.01).curve(V)
        ax3.plot(curve.V * 1000, curve.response, label=f"log $K_f'$ = {logK}")
    ax3.axvline(50.0, color="gray", ls="--", lw=0.8)
    ax3.set_xlabel("EDTA added (mL)")
    ax3.set_ylabel("pM")
    ax3.set_title("EDTA complexometric titration")
    ax3.legend(fontsize=8)


FIGURES = {
    f.__name__: f
    for f in (
        kinetics,
        thermo,
        solutions,
        md,
        statmech,
        quantum,
        spectro,
        structure,
        electrochem,
        photochem,
        surface,
        polymer,
        crystal,
        analytical,
    )
}


def main(names):
    for name in names or FIGURES:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
        FIGURES[name](axes)
        out = OUT_DIR / f"readme_{name}.png"
        fig.savefig(out, dpi=110)
        plt.close(fig)
        print(f"wrote {out}")


if __name__ == "__main__":
    main(sys.argv[1:])
