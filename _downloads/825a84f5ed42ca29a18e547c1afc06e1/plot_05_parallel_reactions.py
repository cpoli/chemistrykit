r"""
Parallel (competing) first-order reactions
============================================

For two competing channels :math:`A \xrightarrow{k_1} B` and
:math:`A \xrightarrow{k_2} C`, the reactant decays with the *sum*
:math:`k_1 + k_2`, while the product ratio :math:`[B]/[C] = k_1/k_2` is
fixed at every instant -- the kinetic basis of product selectivity.
Integrated with
:meth:`~chemistrykit.kinetics.systems.networks.StoichiometricNetwork.parallel`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.systems.networks import StoichiometricNetwork
from chemistrykit.kinetics.visualizers.kinetics_plots import plot_concentration_vs_time

k1, k2 = 2.0, 0.5
net = StoichiometricNetwork.parallel(k1=k1, k2=k2, A0=1.0)
result = net.integrate((0.0, 4.0), dt=1e-3, method="rk4")

fig, ax = plt.subplots(figsize=(7, 5))
plot_concentration_vs_time(result, ax=ax)
ax.plot(result.t, np.exp(-(k1 + k2) * result.t), "k--", linewidth=0.8, label="exp(-(k1+k2)t)")
ax.legend()
ratio = result.concentration("B")[-1] / result.concentration("C")[-1]
ax.set_title(f"Parallel channels: final [B]/[C] = {ratio:.3f} (k1/k2 = {k1 / k2:g})")
fig.tight_layout()
plt.show()
