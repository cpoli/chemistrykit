r"""
Michaelis-Menten kinetics, Lineweaver-Burk, and inhibition
==============================================================

Three related views of the same enzyme: the saturating
Michaelis-Menten rate-vs-substrate curve, its Lineweaver-Burk
double-reciprocal linearization (used to extract :math:`V_{max}` and
:math:`K_m` from noisy initial-rate data), and how competitive vs.
noncompetitive inhibition distort that curve differently.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.systems.enzyme import (
    competitive_inhibition_rate,
    fit_lineweaver_burk,
    michaelis_menten_rate,
    noncompetitive_inhibition_rate,
)
from chemistrykit.kinetics.visualizers.kinetics_plots import plot_lineweaver_burk

Vmax, Km = 10.0, 2.0
S = np.linspace(0.1, 20.0, 200)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

axes[0].plot(S, michaelis_menten_rate(S, Vmax, Km), label="uninhibited", color="steelblue")
axes[0].plot(S, competitive_inhibition_rate(S, I=3.0, Vmax=Vmax, Km=Km, Ki=1.0), label="competitive (I=3)", color="darkorange")
axes[0].plot(S, noncompetitive_inhibition_rate(S, I=3.0, Vmax=Vmax, Km=Km, Ki=1.0), label="noncompetitive (I=3)", color="seagreen")
axes[0].axhline(Vmax, color="gray", linestyle="--", linewidth=0.8)
axes[0].set_xlabel("[S]")
axes[0].set_ylabel("v")
axes[0].set_title("Competitive vs. noncompetitive inhibition")
axes[0].legend()

# %%
# Competitive inhibition can be "out-competed" by enough substrate (all
# three curves approach the same asymptote at high [S]); noncompetitive
# inhibition permanently lowers the ceiling.

# %%
# Lineweaver-Burk: fit (Vmax, Km) back out from noisy discrete samples.

rng = np.random.default_rng(1)
S_samples = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 16.0])
v_exact = michaelis_menten_rate(S_samples, Vmax, Km)
v_noisy = v_exact * (1.0 + rng.normal(0.0, 0.03, size=S_samples.shape))
fit = fit_lineweaver_burk(S_samples, v_noisy)
print(f"true:   Vmax={Vmax}, Km={Km}")
print(f"fitted: Vmax={fit.Vmax:.3f}, Km={fit.Km:.3f}, R^2={fit.r_squared:.5f}")

plot_lineweaver_burk(S_samples, v_noisy, fit=fit, ax=axes[1])
fig.tight_layout()
plt.show()
