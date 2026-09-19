r"""
The Carothers equation
=========================

Step-growth polymerization's number-average degree of polymerization
follows the Carothers equation :math:`\bar X_n=1/(1-p)`, which diverges
as the extent of reaction `p` approaches 1 -- reaching a useful high
molecular weight requires very high conversion. A stoichiometric
imbalance between the two functional groups caps :math:`\bar X_n` even
at complete conversion of the limiting group.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.polymer.systems.step_growth import (
    degree_of_polymerization,
    degree_of_polymerization_stoichiometric_imbalance,
    extent_of_reaction_for_DP,
)
from chemistrykit.polymer.visualizers.polymer_plots import plot_carothers_curve

for p in (0.5, 0.9, 0.95, 0.99, 0.999):
    print(f"p = {p:<6} -> Xn = {degree_of_polymerization(p):8.2f}")

# %%
# Reaching a modest Xn=100 requires 99% conversion, not 90%.
p_for_100 = extent_of_reaction_for_DP(100.0)
print(f"\nExtent of reaction needed for Xn=100: p = {p_for_100:.4f}")

# %%
# A stoichiometric imbalance (r < 1) caps the attainable Xn even as p -> 1.
p_scan = np.linspace(0.0, 0.999, 300)
Xn_balanced = degree_of_polymerization(p_scan)
Xn_imbalanced = degree_of_polymerization_stoichiometric_imbalance(p_scan, r=0.98)

Xn_balanced_at_p999 = degree_of_polymerization(0.999)
Xn_imbalanced_at_p999 = degree_of_polymerization_stoichiometric_imbalance(0.999, r=0.98)
print(f"\nAt p=0.999: Xn (balanced) = {Xn_balanced_at_p999:.1f}, Xn (r=0.98) = {Xn_imbalanced_at_p999:.1f}")

# %%
fig, ax = plt.subplots()
plot_carothers_curve(p_scan, Xn_balanced, ax=ax, label="stoichiometric balance (r=1)")
plot_carothers_curve(p_scan, Xn_imbalanced, ax=ax, linestyle="--", label="r=0.98 (stoichiometric imbalance)")
ax.set_ylim(0, 250)
ax.legend()
plt.tight_layout()
plt.show()
