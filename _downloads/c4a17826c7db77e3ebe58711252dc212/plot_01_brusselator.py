r"""
Prigogine and Lefever's Brusselator limit cycle
==================================================

Below the Hopf-bifurcation threshold (:math:`B < 1 + A^2`) the
Brusselator's single fixed point is a stable focus and any initial
concentration relaxes toward it. Above threshold, that same fixed point
becomes unstable and every trajectory is drawn instead onto a surrounding
stable limit cycle -- sustained chemical oscillation with no external
forcing, driven purely by the autocatalytic :math:`X^2 Y` term.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.kinetics.systems.oscillators import Brusselator
from chemistrykit.kinetics.visualizers.kinetics_plots import plot_concentration_vs_time, plot_phase_portrait

fig, axes = plt.subplots(2, 2, figsize=(11, 8))

# %%
# Below threshold: relaxation to the fixed point (X*, Y*) = (A, B/A).
below = Brusselator(X0=2.0, Y0=0.2, A=1.0, B=1.5)
result_below = below.integrate((0.0, 60.0), dt=1e-2, method="rk4")
plot_concentration_vs_time(result_below, ax=axes[0, 0])
axes[0, 0].set_title(f"Below threshold (B={below.B} < 1+A^2={1 + below.A**2}): relaxes to fixed point")
plot_phase_portrait(result_below, "X", "Y", ax=axes[0, 1])
fp_below = below.fixed_point()
axes[0, 1].plot(*fp_below, "k*", markersize=12, label="fixed point")
axes[0, 1].legend()

# %%
# Above threshold: a stable limit cycle around the now-unstable fixed point.
above = Brusselator(X0=1.0, Y0=1.0, A=1.0, B=3.0)
result_above = above.integrate((0.0, 60.0), dt=1e-2, method="rk4")
plot_concentration_vs_time(result_above, ax=axes[1, 0])
axes[1, 0].set_title(f"Above threshold (B={above.B} > 1+A^2={1 + above.A**2}): sustained oscillation")
plot_phase_portrait(result_above, "X", "Y", ax=axes[1, 1])
fp_above = above.fixed_point()
axes[1, 1].plot(*fp_above, "k*", markersize=12, label="unstable fixed point")
axes[1, 1].legend()

fig.tight_layout()
plt.show()
