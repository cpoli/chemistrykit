r"""
Ideal gas, van der Waals, and Redlich-Kwong P-V isotherms
==========================================================

Comparing three equations of state for carbon dioxide at a fixed
temperature: the ideal gas law
(:class:`~chemistrykit.thermo.systems.equations_of_state.IdealGas`), van
der Waals (:class:`~chemistrykit.thermo.systems.equations_of_state.VanDerWaals`),
and Redlich-Kwong (:class:`~chemistrykit.thermo.systems.equations_of_state.RedlichKwong`).
Both real-gas equations are built directly from CO2's critical
temperature and pressure via :meth:`~chemistrykit.thermo.systems.equations_of_state.VanDerWaals.from_critical_constants`,
so no separately fitted `a`/`b` parameters are needed. The further below
the critical temperature, the more the real-gas isotherms deviate from
the ideal gas law's simple hyperbola.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.thermo.systems.equations_of_state import IdealGas, RedlichKwong, VanDerWaals

# CO2 critical constants (Atkins & de Paula, Physical Chemistry, Table 1.5).
Tc, Pc = 304.13, 7.3773e6

ideal = IdealGas()
vdw = VanDerWaals.from_critical_constants(Tc, Pc)
rk = RedlichKwong.from_critical_constants(Tc, Pc)

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

for ax, T, title in zip(axes, [1.5 * Tc, 0.95 * Tc], ["well above Tc", "just below Tc"]):
    Vm = np.linspace(1.5 * vdw.b, 6.0e-4, 400)
    ax.plot(Vm * 1000.0, ideal.pressure(Vm, T) / 1.0e6, label="ideal gas")
    ax.plot(Vm * 1000.0, vdw.pressure(Vm, T) / 1.0e6, label="van der Waals")
    ax.plot(Vm * 1000.0, rk.pressure(Vm, T) / 1.0e6, label="Redlich-Kwong")
    ax.axhline(Pc / 1.0e6, color="gray", linestyle=":", linewidth=0.8, label="Pc" if title == "well above Tc" else None)
    ax.set_xlabel("Vm (L/mol)")
    ax.set_ylabel("P (MPa)")
    ax.set_title(f"T = {T:.1f} K ({title})")
    ax.set_ylim(0, 3.0 * Pc / 1.0e6)
    ax.legend()

fig.suptitle("CO2 P-V isotherms: ideal gas vs. real-gas equations of state")
fig.tight_layout()

# %%
# Well above the critical temperature all three curves nearly coincide
# (CO2 behaves nearly ideally at low pressure/high temperature); close to
# the critical temperature, the real-gas equations bend over sharply
# where the ideal gas law cannot, since only a cubic (or higher) equation
# of state can produce the S-shaped isotherm needed to describe
# liquid-vapor coexistence.

Vm_c = 3.0 * vdw.b
print(f"van der Waals critical molar volume Vc = 3b = {Vm_c * 1000.0:.4f} L/mol")
print(f"van der Waals critical compressibility factor Zc = {Pc * Vm_c / (vdw.R * Tc):.4f} (exactly 3/8)")

plt.show()
