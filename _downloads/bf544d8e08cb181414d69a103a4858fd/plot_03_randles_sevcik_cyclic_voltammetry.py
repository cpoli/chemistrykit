r"""
The Randles-Ševčík equation: peak current in cyclic voltammetry
=================================================================

Randles and Ševčík (1948) independently showed that the peak current of
a reversible linear-sweep voltammogram grows with the square root of the
scan rate, :math:`i_p = 0.4463\,nFAC\sqrt{nFvD/RT}`. This example
evaluates
:func:`~chemistrykit.electrochem.systems.voltammetry.randles_sevcik_peak_current`
over a range of scan rates, and recovers the diffusion coefficient from
the slope of the :math:`i_p` vs. :math:`\sqrt{v}` plot -- the standard
diagnostic for a diffusion-controlled electrode reaction.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.constants import FARADAY, STANDARD_TEMPERATURE, R
from chemistrykit.electrochem.systems.voltammetry import randles_sevcik_peak_current

n, A, C, D_true = 1, np.pi * (1.5e-3) ** 2, 1.0, 2.4e-9  # 3 mm disc, 1 mM ferrocene-like couple
v = np.array([0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0])
rng = np.random.default_rng(1948)
ip = randles_sevcik_peak_current(v, n, A, C, D_true) * (1 + rng.normal(0.0, 0.01, v.size))
for vi, ipi in zip(v, ip, strict=True):
    print(f"v = {vi * 1e3:6.0f} mV/s: i_p = {ipi * 1e6:6.2f} uA")

# %%
# Slope of i_p vs sqrt(v) through the origin gives D.
slope = np.sum(ip * np.sqrt(v)) / np.sum(v)
D_fit = (slope / (0.4463 * n * FARADAY * A * C)) ** 2 * R * STANDARD_TEMPERATURE / (n * FARADAY)
print(f"\nD true = {D_true:.2e} m^2/s, D from Randles-Sevcik plot = {D_fit:.2e} m^2/s")

x = np.linspace(0.0, 1.0, 50)
fig, ax = plt.subplots()
ax.plot(np.sqrt(v), ip * 1e6, "o", label="peak currents")
ax.plot(x, slope * x * 1e6, label=r"linear in $\sqrt{v}$")
ax.set_xlabel(r"$\sqrt{v}$ (V/s)$^{1/2}$")
ax.set_ylabel(r"$i_p$ ($\mu$A)")
ax.set_title("Randles-Ševčík plot")
ax.legend()
fig.tight_layout()
plt.show()
