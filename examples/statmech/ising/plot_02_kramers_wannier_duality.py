r"""
Kramers-Wannier duality and the exact Ising critical temperature
===================================================================

Kramers and Wannier (1941) showed that the square-lattice Ising model's
high-temperature series at coupling :math:`K=J/k_BT` equals its
low-temperature series at the dual coupling :math:`K^*`, with

.. math::

    \sinh 2K\,\sinh 2K^* = 1,\qquad K^*=-\tfrac12\ln\tanh K .

If there is a single phase transition, it must sit at the self-dual point
:math:`K_c=K_c^*`, which gives :math:`\sinh 2K_c=1` and
:math:`T_c=2J/[k_B\ln(1+\sqrt2)]\approx2.269\,J/k_B`. This was the first
exact critical temperature for a lattice model.
:func:`~chemistrykit.statmech.kramers_wannier_dual_coupling` and
:func:`~chemistrykit.statmech.ising_2d_critical_temperature` implement
both.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.constants import K_B
from chemistrykit.statmech import ising_2d_critical_temperature, kramers_wannier_dual_coupling

K = np.linspace(0.05, 1.5, 400)
K_dual = kramers_wannier_dual_coupling(K)
K_c = 0.5 * np.log(1.0 + np.sqrt(2.0))

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
axes[0].plot(K, K_dual, color="steelblue", label=r"dual coupling $K^*(K)$")
axes[0].plot(K, K, color="gray", linestyle="--", label=r"$K^* = K$")
axes[0].plot([K_c], [K_c], "o", color="crimson", label=rf"self-dual point $K_c$ = {K_c:.5f}")
# One pair of dual points: a high-T state (small K) maps to a low-T state (large K) and back.
K_hot = 0.2
K_cold = float(kramers_wannier_dual_coupling(K_hot))
axes[0].annotate("", xy=(K_hot, K_cold), xytext=(K_cold, K_hot), arrowprops={"arrowstyle": "<->", "color": "darkorange"})
axes[0].plot([K_hot, K_cold], [K_cold, K_hot], "s", color="darkorange", label=f"dual pair {K_hot} <-> {K_cold:.3f}")
axes[0].set_xlim(0.0, 1.5)
axes[0].set_ylim(0.0, 1.5)
axes[0].set_xlabel(r"$K = J / k_BT$")
axes[0].set_ylabel(r"$K^*$")
axes[0].set_title("High temperature maps onto low temperature")
axes[0].legend(fontsize=8)

axes[1].plot(K, np.sinh(2 * K) * np.sinh(2 * K_dual), color="steelblue")
axes[1].set_ylim(0.9, 1.1)
axes[1].axvline(K_c, color="crimson", linestyle=":", label=r"$\sinh 2K_c = 1$")
axes[1].set_xlabel(r"$K$")
axes[1].set_ylabel(r"$\sinh 2K \, \sinh 2K^*$")
axes[1].set_title("The duality relation holds identically")
axes[1].legend()
fig.tight_layout()

# %%
# The critical temperature for a coupling J/k_B = 100 K, and the checks
# that K_c is the fixed point with sinh(2 K_c) = 1:

J = 100.0 * K_B
print(f"T_c = {ising_2d_critical_temperature(J):.3f} K  (= 2.269 J/k_B)")
print(f"K_c = {K_c:.6f}, K*(K_c) = {float(kramers_wannier_dual_coupling(K_c)):.6f}, sinh(2 K_c) = {np.sinh(2 * K_c):.6f}")

plt.show()
