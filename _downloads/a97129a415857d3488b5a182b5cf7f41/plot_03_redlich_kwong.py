r"""
Redlich-Kwong isotherms of CO2 against van der Waals and the ideal gas
======================================================================

Redlich and Kwong divided van der Waals's attraction term by
:math:`\sqrt{T}`. This example draws CO2 P-V isotherms well above and
just below the critical temperature from
:class:`~chemistrykit.thermo.systems.equations_of_state.RedlichKwong`,
with :class:`~chemistrykit.thermo.systems.equations_of_state.VanDerWaals`
and :class:`~chemistrykit.thermo.systems.equations_of_state.IdealGas` for
comparison. Both real-gas equations are built from CO2's critical
temperature and pressure alone via
:meth:`~chemistrykit.thermo.systems.equations_of_state.RedlichKwong.from_critical_constants`.
Redlich-Kwong's universal critical compressibility factor is 1/3, closer
to real gases than van der Waals's 3/8.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.thermo.systems.equations_of_state import IdealGas, RedlichKwong, VanDerWaals

Tc, Pc = 304.13, 7.3773e6  # CO2
ideal = IdealGas()
vdw = VanDerWaals.from_critical_constants(Tc, Pc)
rk = RedlichKwong.from_critical_constants(Tc, Pc)

fig, axes = plt.subplots(1, 2, figsize=(11, 5))
for ax, T, title in zip(axes, [1.5 * Tc, 0.95 * Tc], ["well above Tc", "just below Tc"], strict=True):
    Vm = np.linspace(1.5 * vdw.b, 6.0e-4, 400)
    ax.plot(Vm * 1000.0, rk.pressure(Vm, T) / 1.0e6, color="crimson", linewidth=2, label="Redlich-Kwong")
    ax.plot(Vm * 1000.0, vdw.pressure(Vm, T) / 1.0e6, "--", label="van der Waals")
    ax.plot(Vm * 1000.0, ideal.pressure(Vm, T) / 1.0e6, ":", label="ideal gas")
    ax.axhline(Pc / 1.0e6, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel(r"$V_m$ (L/mol)")
    ax.set_ylabel("P (MPa)")
    ax.set_title(f"T = {T:.1f} K ({title})")
    ax.set_ylim(0, 3.0 * Pc / 1.0e6)
    ax.legend()
fig.suptitle("CO2 isotherms: Redlich-Kwong vs. van der Waals vs. ideal gas")
fig.tight_layout()

# %%
# The critical-point construction fixes Redlich-Kwong's critical molar
# volume at :math:`RT_c/(3P_c)`; the critical isotherm passes through it
# at :math:`P = P_c`.

Vc_rk = rk.R * Tc / (3.0 * Pc)
print(f"Redlich-Kwong: P(Vc, Tc)/Pc = {rk.pressure(Vc_rk, Tc) / Pc:.4f}, Zc = {Pc * Vc_rk / (rk.R * Tc):.4f}")
print(f"van der Waals: Zc = {Pc * 3.0 * vdw.b / (vdw.R * Tc):.4f}; measured CO2 Zc is about 0.274")

plt.show()
