r"""
Peng-Robinson: vapor pressures from the acentric factor
=======================================================

Peng and Robinson tuned the temperature dependence of their attraction
term, :math:`\alpha(T)`, to Pitzer's acentric factor :math:`\omega` so
that the equation reproduces pure-fluid vapor pressures. This example
computes the CO2 vapor-pressure curve predicted by
:class:`~chemistrykit.thermo.systems.equations_of_state.PengRobinson`
(:math:`\omega = 0.224`) with
:func:`~chemistrykit.thermo.systems.fugacity.saturation_pressure`, and
compares it with van der Waals and Redlich-Kwong built from the same
critical constants. Only Peng-Robinson passes through Pitzer's defining
point, :math:`\log_{10}(P_{sat}/P_c) = -1 - \omega` at
:math:`T/T_c = 0.7`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.thermo.systems.equations_of_state import PengRobinson, RedlichKwong, VanDerWaals
from chemistrykit.thermo.systems.fugacity import saturation_pressure

Tc, Pc, omega = 304.13, 7.3773e6, 0.224  # CO2
models = {
    "Peng-Robinson": PengRobinson(Tc, Pc, omega),
    "Redlich-Kwong": RedlichKwong.from_critical_constants(Tc, Pc),
    "van der Waals": VanDerWaals.from_critical_constants(Tc, Pc),
}

Tr = np.linspace(0.6, 0.98, 25)
fig, ax = plt.subplots(figsize=(7, 5))
for name, eos in models.items():
    Pr = [saturation_pressure(eos, t * Tc) / Pc for t in Tr]
    ax.plot(1.0 / Tr, np.log10(Pr), label=name)
ax.plot([1.0 / 0.7], [-1.0 - omega], "k*", markersize=12, label=r"acentric-factor point ($\omega$ = 0.224)")
ax.set_xlabel(r"$T_c / T$")
ax.set_ylabel(r"$\log_{10}(P_{sat}/P_c)$")
ax.set_title("CO2 vapor-pressure curve from three cubic equations of state")
ax.legend()
fig.tight_layout()

# %%
# The acentric factor implied by each equation at :math:`T/T_c = 0.7`:

for name, eos in models.items():
    implied = -1.0 - np.log10(saturation_pressure(eos, 0.7 * Tc) / Pc)
    print(f"{name:14s}: implied omega = {implied:+.3f}")
print(f"Peng-Robinson Zc = 0.3074 (critical isotherm check: P(Vc)/Pc = {models['Peng-Robinson'].pressure(0.30740 * 8.314462618 * Tc / Pc, Tc) / Pc:.4f})")

plt.show()
