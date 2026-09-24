r"""
Wilhelmy's sucrose inversion: the first measured first-order rate law
========================================================================

Wilhelmy (1850) followed the acid-catalyzed inversion of cane sugar by
polarimetry: the optical rotation :math:`\alpha(t)` drifts from the
dextrorotatory sucrose value towards the (levorotatory) value of the
glucose + fructose mixture. Because :math:`\alpha` is linear in the
concentrations, :math:`\alpha(t) - \alpha_\infty` is proportional to the
sucrose still present, and a straight line of
:math:`\ln[(\alpha - \alpha_\infty)/(\alpha_0 - \alpha_\infty)]` against
:math:`t` identifies :math:`-d[\text{sugar}]/dt = k[\text{sugar}]`.

The polarimeter readings below are synthetic (an illustrative rate
constant and rotation angles, plus measurement noise), generated from
:class:`~chemistrykit.kinetics.systems.rate_laws.FirstOrder`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.systems.rate_laws import FirstOrder

k_true = 0.012  # min^-1 (illustrative)
C0 = 1.0  # relative sucrose concentration
sucrose = FirstOrder(k=k_true, C0=C0)

# Rotation of the solution: +alpha_0 for pure sucrose, alpha_inf once fully
# inverted (fructose's strong levorotation wins over glucose).
alpha_0, alpha_inf = 13.0, -4.0  # degrees (illustrative)

rng = np.random.default_rng(1850)
t = np.linspace(0.0, 240.0, 17)  # minutes
fraction_left = sucrose.concentration(t) / C0
alpha = alpha_inf + (alpha_0 - alpha_inf) * fraction_left + rng.normal(0.0, 0.08, t.shape)

# %%
# Wilhelmy's test: the log of the remaining rotation is linear in time.

y = np.log((alpha - alpha_inf) / (alpha_0 - alpha_inf))
slope, intercept = np.polyfit(t, y, 1)
k_fit = -slope
print(f"true k = {k_true:.4f} /min, fitted k = {k_fit:.4f} /min")
print(f"half-life ln2/k = {np.log(2) / k_fit:.1f} min")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
t_fine = np.linspace(0.0, 240.0, 300)
axes[0].plot(t, alpha, "o", color="steelblue", label="polarimeter reading")
axes[0].plot(t_fine, alpha_inf + (alpha_0 - alpha_inf) * np.exp(-k_fit * t_fine), color="darkorange", label="first-order fit")
axes[0].axhline(0.0, color="gray", linewidth=0.8)
axes[0].axhline(alpha_inf, color="gray", linestyle="--", linewidth=0.8, label=r"$\alpha_\infty$ (invert sugar)")
axes[0].set_xlabel("t (min)")
axes[0].set_ylabel(r"optical rotation $\alpha$ (deg)")
axes[0].set_title("Rotation falls through zero as sucrose inverts")
axes[0].legend()

axes[1].plot(t, y, "o", color="steelblue")
axes[1].plot(t_fine, intercept + slope * t_fine, color="darkorange", label=f"slope = -k = {slope:.4f} /min")
axes[1].set_xlabel("t (min)")
axes[1].set_ylabel(r"$\ln[(\alpha-\alpha_\infty)/(\alpha_0-\alpha_\infty)]$")
axes[1].set_title("A straight line: first-order kinetics")
axes[1].legend()
fig.tight_layout()

# %%
# The half-life does not depend on how much sugar is left: every
# :math:`t_{1/2}` interval halves the remaining sucrose again.

for n in range(4):
    tn = n * sucrose.half_life()
    print(f"after {n} half-lives ({tn:6.1f} min): [sugar]/[sugar]_0 = {sucrose.concentration(tn) / C0:.4f}")

plt.show()
