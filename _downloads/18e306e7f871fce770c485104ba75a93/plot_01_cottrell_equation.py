r"""
The Cottrell equation: current decay after a potential step
=============================================================

Cottrell (1903) solved the diffusion problem for a planar electrode whose
potential is suddenly stepped so that the electroactive species is
consumed at the surface: the current decays as :math:`t^{-1/2}`,
:math:`i = nFAC\sqrt{D/(\pi t)}`. This example plots
:func:`~chemistrykit.electrochem.systems.voltammetry.cottrell_current`,
checks the constant :math:`i\sqrt{t}` signature, and recovers the
diffusion coefficient from a noisy "measured" transient via the Cottrell
plot of :math:`i` against :math:`t^{-1/2}`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.constants import FARADAY
from chemistrykit.electrochem.systems.voltammetry import cottrell_current

n, A, C, D_true = 1, 7.07e-6, 1.0, 7.0e-10  # 3 mm disc, 1 mM, D in m^2/s
t = np.linspace(0.05, 5.0, 80)
i_true = cottrell_current(t, n, A, C, D_true)
print(
    f"i*sqrt(t) at t = 0.1 s and 4 s: {cottrell_current(0.1, n, A, C, D_true) * np.sqrt(0.1):.4e}, {cottrell_current(4.0, n, A, C, D_true) * 2.0:.4e} A s^1/2"
)

# %%
# Cottrell plot: i vs t^-1/2 is a line through the origin with slope
# nFAC sqrt(D/pi).
rng = np.random.default_rng(1903)
i_meas = i_true * (1 + rng.normal(0.0, 0.01, t.size))
slope = np.sum(i_meas * t**-0.5) / np.sum(t**-1.0)
D_fit = np.pi * (slope / (n * FARADAY * A * C)) ** 2
print(f"D true = {D_true:.2e} m^2/s, D from Cottrell plot = {D_fit:.2e} m^2/s")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(t, i_meas * 1e6, ".", label="measured")
ax1.plot(t, i_true * 1e6, label="Cottrell")
ax1.set_xlabel("t (s)")
ax1.set_ylabel(r"i ($\mu$A)")
ax1.set_title(r"$t^{-1/2}$ decay after a potential step")
ax1.legend()
ax2.plot(t**-0.5, i_meas * 1e6, ".")
ax2.plot(t**-0.5, slope * t**-0.5 * 1e6)
ax2.set_xlabel(r"$t^{-1/2}$ (s$^{-1/2}$)")
ax2.set_ylabel(r"i ($\mu$A)")
ax2.set_title("Cottrell plot")
fig.tight_layout()
plt.show()
