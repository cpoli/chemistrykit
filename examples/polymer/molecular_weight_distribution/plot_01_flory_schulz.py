r"""
The Flory-Schulz molecular-weight distribution
=================================================

For an ideal step-growth polymerization at extent of reaction `p`, the
Flory-Schulz (most-probable) distribution gives the number and weight
fraction of chains of each length in exact closed form. Its
number-/weight-average degrees of polymerization,
:math:`\bar X_n=1/(1-p)` and :math:`\bar X_w=(1+p)/(1-p)`, are
cross-checked here against direct numerical summation of the
distribution, and the PDI :math:`=1+p` is shown approaching exactly 2 as
`p` approaches 1.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.polymer.systems.molecular_weight_distribution import (
    flory_schulz_number_average_DP,
    flory_schulz_number_fraction,
    flory_schulz_pdi,
    flory_schulz_weight_average_DP,
    flory_schulz_weight_fraction,
)
from chemistrykit.polymer.utils.moments import number_average, weight_average
from chemistrykit.polymer.visualizers.polymer_plots import plot_molecular_weight_distribution

p = 0.95
x = np.arange(1, 20000)
N_x = flory_schulz_number_fraction(x, p)
w_x = flory_schulz_weight_fraction(x, p)

print(f"Number fractions sum to {np.sum(N_x):.6f} (should be 1)")
print(f"Weight fractions sum to {np.sum(w_x):.6f} (should be 1)")

# %%
# Cross-check the closed-form averages against direct numerical
# summation of the distribution.
Xn_closed = flory_schulz_number_average_DP(p)
Xn_numeric = number_average(x, N_x)
Xw_closed = flory_schulz_weight_average_DP(p)
Xw_numeric = weight_average(x, N_x)
print(f"\nXn: closed form = {Xn_closed:.4f}, numerical sum = {Xn_numeric:.4f}")
print(f"Xw: closed form = {Xw_closed:.4f}, numerical sum = {Xw_numeric:.4f}")

# %%
# PDI approaches exactly 2 as p -> 1.
p_values = np.array([0.5, 0.9, 0.99, 0.999, 1.0])
for pv in p_values:
    print(f"p={pv:<7} PDI = {flory_schulz_pdi(pv):.6f}")

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
plot_molecular_weight_distribution(x[:150], N_x[:150], ax=axes[0], label="number fraction")
plot_molecular_weight_distribution(x[:150], w_x[:150], ax=axes[1], label="weight fraction", color="crimson")
plt.tight_layout()
plt.show()
