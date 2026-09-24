r"""
Ising's one-dimensional chain: no phase transition
=====================================================

Ernst Ising (1925) solved the chain of spins :math:`s_i=\pm1` with
nearest-neighbor coupling `J` and found no spontaneous magnetization at
any temperature above zero. :class:`~chemistrykit.statmech.Ising1D`
solves it with the 2x2 transfer matrix, whose eigenvalues
:math:`\lambda_\pm` give :math:`Z_N=\lambda_+^N+\lambda_-^N` exactly. The
plots show a checked partition function, a magnetization that is a smooth
function of field at every temperature (it only becomes a step at
:math:`T=0`), a correlation length that grows as :math:`e^{2J/k_BT}` but
stays finite, and a smooth heat-capacity maximum instead of a singularity.
"""

# %%
import itertools

import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.constants import K_B
from chemistrykit.statmech import Ising1D

J = 100.0 * K_B  # coupling, J/k_B = 100 K

# %%
# The transfer-matrix result equals brute-force enumeration of all 2^N states:

N, T_check, h_check = 10, 150.0, 20.0 * K_B
Z_exact = 0.0
for spins in itertools.product((-1, 1), repeat=N):
    s = np.array(spins)
    E = -J * np.sum(s * np.roll(s, 1)) - h_check * np.sum(s)
    Z_exact += np.exp(-E / (K_B * T_check))
print(f"brute force Z_10 = {Z_exact:.6e}, transfer matrix Z_10 = {Ising1D(J, h_check).partition_function(T_check, N):.6e}")

# %%
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
h = np.linspace(-50.0, 50.0, 400) * K_B
for T, color in [(25.0, "crimson"), (50.0, "darkorange"), (100.0, "seagreen"), (300.0, "steelblue")]:
    m = np.array([Ising1D(J, hi).magnetization(T) for hi in h])
    axes[0].plot(h / K_B, m, color=color, label=f"T = {T:.0f} K")
axes[0].set_xlabel(r"field $h/k_B$ (K)")
axes[0].set_ylabel("magnetization per spin")
axes[0].set_title("m(h) is smooth for every T > 0")
axes[0].legend()

T = np.linspace(20.0, 400.0, 300)
chain = Ising1D(J)
axes[1].semilogy(T, chain.correlation_length(T), color="steelblue", label=r"exact $\xi = -1/\ln\tanh K$")
axes[1].semilogy(T, 0.5 * np.exp(2 * J / (K_B * T)), color="gray", linestyle="--", label=r"$\frac{1}{2}e^{2J/k_BT}$")
axes[1].set_xlabel("T (K)")
axes[1].set_ylabel("correlation length (sites)")
axes[1].set_title("Finite correlation length")
axes[1].legend()

axes[2].plot(T, chain.heat_capacity_per_spin(T) / K_B, color="crimson")
axes[2].set_xlabel("T (K)")
axes[2].set_ylabel(r"$c / k_B$ per spin")
axes[2].set_title("Smooth heat capacity: no transition")
fig.tight_layout()

plt.show()
