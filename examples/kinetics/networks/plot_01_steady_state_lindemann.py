r"""
Bodenstein's steady-state approximation and Lindemann's unimolecular fall-off
================================================================================

Bodenstein (1913) set the net rate of change of a short-lived
intermediate to zero; Lindemann (1922) used exactly that trick on the
energized molecule :math:`A^*` of the mechanism

.. math::

    A + M \xrightarrow{k_1} A^* + M, \qquad
    A^* + M \xrightarrow{k_{-1}} A + M, \qquad
    A^* \xrightarrow{k_2} P

to explain why "unimolecular" gas reactions become second order at low
pressure. The first panel shows the steady-state approximation for a
chain :math:`A \to B \to C` becoming exact as the intermediate is made
short-lived; the second integrates Lindemann's full mechanism with
:class:`~chemistrykit.kinetics.systems.networks.StoichiometricNetwork` at
many bath-gas concentrations :math:`[M]` and compares the measured
first-order rate constant with the steady-state prediction

.. math::

    k_{uni} = \frac{k_1 k_2 [M]}{k_{-1}[M] + k_2}.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.systems.networks import (
    StoichiometricNetwork,
    consecutive_analytic,
    ssa_intermediate_concentration,
)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

# %%
# Steady-state approximation for A -> B -> C: as k2/k1 grows, the exact
# intermediate concentration is tracked ever more closely by
# :math:`[B]_{ssa} = (k_1/k_2)[A](t)`.

t = np.linspace(0.0, 8.0, 300)
for k2 in (2.0, 5.0, 25.0):
    _, B_exact, _ = consecutive_analytic(A0=1.0, k1=1.0, k2=k2, t=t)
    B_ssa = ssa_intermediate_concentration(A0=1.0, k1=1.0, k2=k2, t=t)
    (line,) = axes[0].plot(t, B_exact, label=f"[B] exact, k2/k1={k2:g}")
    axes[0].plot(t, B_ssa, "--", color=line.get_color(), label=f"[B] SSA, k2/k1={k2:g}")
axes[0].set_xlabel("t")
axes[0].set_ylabel("[B]")
axes[0].set_title("Bodenstein: steady-state approximation")
axes[0].legend(fontsize=8)

# %%
# Lindemann's mechanism with [M] held constant (the bath gas is in huge
# excess), so ``k1*[M]`` and ``k_-1*[M]`` act as pseudo-first-order rate
# constants. Species order (A, A*, P); reactions (activation,
# deactivation, decomposition).

k1, k_m1, k2 = 1.0, 100.0, 100.0
species = ("A", "A*", "P")
stoich = [[-1.0, 1.0, 0.0], [1.0, -1.0, -1.0], [0.0, 0.0, 1.0]]
orders = [[1.0, 0.0, 0.0], [0.0, 1.0, 1.0], [0.0, 0.0, 0.0]]

M_values = np.logspace(-3, 3, 9)
k_measured = []
for M in M_values:
    k_ssa = k1 * k2 * M / (k_m1 * M + k2)
    net = StoichiometricNetwork(species, stoich, [k1 * M, k_m1 * M, k2], orders, state0=[1.0, 0.0, 0.0])
    result = net.integrate((0.0, 3.0 / k_ssa), method="dopri5", rtol=1e-8, atol=1e-12)
    late = result.t > 0.5 / k_ssa
    k_measured.append(-np.polyfit(result.t[late], np.log(result.concentration("A")[late]), 1)[0])
k_measured = np.array(k_measured)

M_fine = np.logspace(-3, 3, 200)
axes[1].loglog(M_fine, k1 * k2 * M_fine / (k_m1 * M_fine + k2), color="darkorange", label="steady-state $k_{uni}$")
axes[1].loglog(M_values, k_measured, "o", color="steelblue", label="from full integration")
axes[1].loglog(M_fine, k1 * M_fine, ":", color="gray", label=r"low pressure: $k_1[M]$ (2nd order)")
axes[1].axhline(k1 * k2 / k_m1, color="gray", linestyle="--", linewidth=0.8, label=r"high pressure: $k_1k_2/k_{-1}$")
axes[1].set_ylim(5e-4, 3.0)
axes[1].set_xlabel("[M]")
axes[1].set_ylabel("effective first-order rate constant")
axes[1].set_title("Lindemann fall-off curve")
axes[1].legend(fontsize=8)

for M, k_meas in zip(M_values[::4], k_measured[::4]):
    print(f"[M]={M:9.3g}: measured k_uni={k_meas:.4f}, SSA k_uni={k1 * k2 * M / (k_m1 * M + k2):.4f}")

fig.tight_layout()
plt.show()
