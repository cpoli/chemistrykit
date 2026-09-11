r"""
Parallel, consecutive, and reversible reaction networks
============================================================

:class:`~chemistrykit.kinetics.systems.networks.StoichiometricNetwork`
integrates an arbitrary mass-action reaction network numerically via
:mod:`chemistrykit.integrators`. Three classic mechanisms -- built from
its named constructors -- are each checked here against a closed-form
solution, and the steady-state approximation for the consecutive chain
is shown converging as the second step is made much faster than the
first.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.systems.networks import (
    StoichiometricNetwork,
    consecutive_analytic,
    reversible_analytic,
    ssa_intermediate_concentration,
)
from chemistrykit.kinetics.visualizers.kinetics_plots import plot_concentration_vs_time

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# %%
# Parallel: A -> B (k1), A -> C (k2). The product ratio [B]/[C] is
# exactly k1/k2 at every instant.
k1, k2 = 2.0, 0.5
net = StoichiometricNetwork.parallel(k1=k1, k2=k2, A0=1.0)
result = net.integrate((0.0, 4.0), dt=1e-3, method="rk4")
plot_concentration_vs_time(result, ax=axes[0])
axes[0].set_title(f"Parallel (final [B]/[C] = {result.concentration('B')[-1] / result.concentration('C')[-1]:.3f}, expected {k1 / k2})")

# %%
# Consecutive: A -> B -> C, checked against the Bateman solution.
k1, k2, A0 = 1.0, 0.3, 1.0
net = StoichiometricNetwork.consecutive(k1=k1, k2=k2, A0=A0)
result = net.integrate((0.0, 15.0), dt=1e-3, method="rk4")
A_exact, B_exact, C_exact = consecutive_analytic(A0, k1, k2, result.t)
plot_concentration_vs_time(result, ax=axes[1])
axes[1].plot(result.t, B_exact, "k--", linewidth=1, label="[B] analytic")
axes[1].legend()
axes[1].set_title("Consecutive A -> B -> C")

# %%
# Reversible: A <-> B, relaxing to a Boltzmann-like equilibrium ratio kf/kr.
kf, kr, A0 = 2.0, 1.0, 1.0
net = StoichiometricNetwork.reversible(kf=kf, kr=kr, A0=A0)
result = net.integrate((0.0, 8.0), dt=1e-3, method="rk4")
A_exact, B_exact = reversible_analytic(A0, kf, kr, result.t)
plot_concentration_vs_time(result, ax=axes[2])
axes[2].axhline(kr * A0 / (kf + kr), color="gray", linestyle="--", linewidth=0.8, label="[A]_eq")
axes[2].legend()
axes[2].set_title("Reversible A <-> B")

fig.tight_layout()

# %%
# Steady-state approximation: as k2/k1 grows, the intermediate's exact
# concentration is tracked ever more closely by the SSA formula
# [B]_ssa = (k1/k2)*[A](t).
fig2, ax2 = plt.subplots(figsize=(7, 5))
t = np.linspace(0.0, 8.0, 300)
for k2 in (2.0, 5.0, 25.0):
    _, B_exact, _ = consecutive_analytic(A0=1.0, k1=1.0, k2=k2, t=t)
    B_ssa = ssa_intermediate_concentration(A0=1.0, k1=1.0, k2=k2, t=t)
    (line,) = ax2.plot(t, B_exact, label=f"[B] exact, k2={k2}")
    ax2.plot(t, B_ssa, "--", color=line.get_color(), label=f"[B] SSA, k2={k2}")
ax2.set_xlabel("t")
ax2.set_ylabel("[B]")
ax2.set_title("Steady-state approximation converging as k2/k1 grows")
ax2.legend(fontsize=8)
fig2.tight_layout()

plt.show()
