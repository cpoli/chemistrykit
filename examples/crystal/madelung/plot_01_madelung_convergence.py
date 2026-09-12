r"""
Converging the NaCl Madelung constant
========================================

The NaCl Madelung-constant lattice sum is only conditionally convergent:
a naive truncated sum over a growing cube of ions does not settle down as
the cutoff grows, it oscillates. :func:`~chemistrykit.crystal.systems.madelung.madelung_constant_nacl`
uses Evjen's method (fractional boundary-charge weighting) to build a
genuinely converging summation instead, and this example demonstrates
both behaviors side by side against the literature value.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.crystal.systems.madelung import MADELUNG_CONSTANT_NACL_LITERATURE, madelung_constant_nacl
from chemistrykit.crystal.visualizers.crystal_plots import plot_madelung_convergence


def naive_truncated_sum(n_shells: int) -> float:
    """A naive (unweighted) truncated lattice sum -- does NOT converge properly."""
    total = 0.0
    for i in range(-n_shells, n_shells + 1):
        for j in range(-n_shells, n_shells + 1):
            for k in range(-n_shells, n_shells + 1):
                if i == 0 and j == 0 and k == 0:
                    continue
                total += (-1) ** (i + j + k) / np.sqrt(i * i + j * j + k * k)
    return -total


n_shells = list(range(2, 16))
evjen_values = [madelung_constant_nacl(n) for n in n_shells]
naive_values = [naive_truncated_sum(n) for n in n_shells]

print(f"Literature value: {MADELUNG_CONSTANT_NACL_LITERATURE}")
print(f"\n{'n_shells':>8s}  {'Evjen (converges)':>18s}  {'naive (oscillates)':>18s}")
for n, e, nv in zip(n_shells, evjen_values, naive_values, strict=True):
    print(f"{n:8d}  {e:18.6f}  {nv:18.6f}")

# %%
# The Evjen-weighted sum settles to within 1e-4 of the literature value
# well before the naive sum's oscillation amplitude has shrunk at all:
print(f"\nEvjen std (n=10..15): {np.std(evjen_values[-6:]):.2e}")
print(f"Naive std (n=10..15): {np.std(naive_values[-6:]):.2e}")

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
plot_madelung_convergence(n_shells, evjen_values, literature_value=MADELUNG_CONSTANT_NACL_LITERATURE, ax=axes[0])
axes[0].set_title("Evjen's method (converges)")
plot_madelung_convergence(n_shells, naive_values, literature_value=MADELUNG_CONSTANT_NACL_LITERATURE, ax=axes[1])
axes[1].set_title("Naive truncation (oscillates)")
plt.tight_layout()
plt.show()
