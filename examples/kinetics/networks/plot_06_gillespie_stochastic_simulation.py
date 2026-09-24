r"""
Gillespie's stochastic simulation of reaction kinetics
=========================================================

Gillespie (1976-1977) showed how to sample *exact* trajectories of a
reaction network treated as random, discrete molecular events: draw an
exponential waiting time with rate equal to the total propensity, then
pick which reaction fires in proportion to its propensity.
:func:`~chemistrykit.kinetics.systems.stochastic.gillespie_ssa` implements
this direct method. With thousands of molecules the trajectories hug the
deterministic rate-equation solution; with a handful they are visibly
noisy, and the average over many runs recovers the deterministic curve.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.systems.networks import consecutive_analytic
from chemistrykit.kinetics.systems.stochastic import gillespie_ssa

# A -> B (c1), B -> C (c2): counts of molecules, not concentrations.
stoich = [[-1, 0], [1, -1], [0, 1]]
orders = [[1, 0], [0, 1], [0, 0]]
c1, c2 = 1.0, 0.3
species = ("A", "B", "C")
t_grid = np.linspace(0.0, 15.0, 301)
rng = np.random.default_rng(1976)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
for ax, N in zip(axes, (20, 2000)):
    _, B_det, _ = consecutive_analytic(1.0, c1, c2, t_grid)
    runs = []
    for i in range(100):
        traj = gillespie_ssa(stoich, [c1, c2], orders, [N, 0, 0], t_max=15.0, species=species, seed=rng)
        runs.append(traj.sample(t_grid)[:, 1] / N)
        if i < 5:
            ax.step(traj.t, traj.count("B") / N, where="post", color="steelblue", alpha=0.5, linewidth=0.8)
    ax.plot(t_grid, np.mean(runs, axis=0), color="crimson", label="mean of 100 SSA runs")
    ax.plot(t_grid, B_det, "k--", label="deterministic rate equations")
    ax.set_xlabel("t")
    ax.set_ylabel("[B] / N_A(0)")
    ax.set_title(f"A -> B -> C with N = {N} molecules (5 sample paths)")
    ax.legend()
    print(f"N={N}: max |mean SSA - deterministic| = {np.max(np.abs(np.mean(runs, axis=0) - B_det)):.4f}")

fig.tight_layout()
plt.show()
