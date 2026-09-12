r"""
Ideal vs. real chain-size scaling
====================================

:class:`~chemistrykit.polymer.systems.chain_statistics.IdealChain` gives
the *exact* random-walk result :math:`\langle R^2\rangle=nb^2`.
:class:`~chemistrykit.polymer.systems.chain_statistics.RealChain` scales
instead as :math:`R\sim bn^\nu`, with the Flory exponent :math:`\nu`
set by solvent quality: a good solvent swells the chain
(:math:`\nu\approx0.6`), a poor solvent collapses it (:math:`\nu=1/3`),
and the theta solvent reproduces the ideal chain exactly
(:math:`\nu=1/2`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.polymer.systems.chain_statistics import IdealChain, RealChain, flory_exponent
from chemistrykit.polymer.visualizers.polymer_plots import plot_chain_scaling

b = 0.5  # segment length
n_values = np.logspace(1, 5, 30)

ideal = IdealChain()
theta = RealChain.theta_solvent()
good = RealChain.good_solvent()
poor = RealChain.poor_solvent()

print(f"Flory exponents: theta={flory_exponent('theta')}, good={flory_exponent('good'):.4f}, poor={flory_exponent('poor'):.4f}")

# %%
# The ideal chain's <R^2> scales exactly linearly with n.
n1, n2 = 1000.0, 4000.0
ratio = ideal.mean_square_end_to_end(n2, b) / ideal.mean_square_end_to_end(n1, b)
print(f"\n<R^2>(4000) / <R^2>(1000) for the ideal chain: {ratio:.6f} (exactly 4.0)")

# theta-solvent RealChain reproduces the ideal chain's end-to-end distance exactly.
n_test = 5000.0
match = np.isclose(theta.end_to_end_distance(n_test, b), ideal.end_to_end_distance(n_test, b))
print(f"theta-solvent RealChain matches IdealChain exactly: {bool(match)}")

# %%
# A good solvent swells the chain; a poor solvent collapses it, relative
# to the ideal chain of the same length.
n_test = 10000.0
print(f"\nAt n={n_test:.0f}, b={b}:")
print(f"  ideal chain R = {ideal.end_to_end_distance(n_test, b):.2f}")
print(f"  good-solvent R = {good.end_to_end_distance(n_test, b):.2f}  (swollen)")
print(f"  poor-solvent R = {poor.end_to_end_distance(n_test, b):.2f}  (collapsed)")

# %%
models = {"ideal (theta)": ideal, "good solvent": good, "poor solvent": poor}
ax = plot_chain_scaling(models, n_values, b)
plt.tight_layout()
plt.show()
