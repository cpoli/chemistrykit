r"""
Michaelis and Menten's saturating enzyme rate law
====================================================

Michaelis and Menten (1913) showed that an enzyme-catalyzed rate
saturates with substrate, :math:`v = V_{max}[S]/(K_m + [S])`, reaching
exactly half of :math:`V_{max}` at :math:`[S] = K_m`. This example draws
that saturation curve with its two limiting regimes, integrates the
substrate-depletion progress curve an assay actually records, and shows
how competitive and noncompetitive inhibitors distort the curve.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.systems.enzyme import (
    MichaelisMentenProgress,
    competitive_inhibition_rate,
    michaelis_menten_rate,
    noncompetitive_inhibition_rate,
)
from chemistrykit.kinetics.visualizers.kinetics_plots import plot_concentration_vs_time

Vmax, Km = 10.0, 2.0
S = np.linspace(0.0, 20.0, 400)

fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

# %%
# The saturation curve: first order in [S] well below :math:`K_m`
# (:math:`v \approx (V_{max}/K_m)[S]`), zero order well above it
# (:math:`v \approx V_{max}`), and exactly :math:`V_{max}/2` at :math:`K_m`.

axes[0].plot(S, michaelis_menten_rate(S, Vmax, Km), color="steelblue", label="Michaelis-Menten")
axes[0].plot(S[S < 4], Vmax / Km * S[S < 4], ":", color="gray", label=r"$(V_{max}/K_m)[S]$")
axes[0].axhline(Vmax, color="gray", linestyle="--", linewidth=0.8, label=r"$V_{max}$")
axes[0].plot([Km, Km, 0], [0, Vmax / 2, Vmax / 2], color="crimson", linewidth=0.8)
axes[0].plot(Km, michaelis_menten_rate(Km, Vmax, Km), "o", color="crimson", label=r"$v(K_m) = V_{max}/2$")
axes[0].set_xlabel("[S]")
axes[0].set_ylabel("v")
axes[0].set_title("Saturation of the enzyme")
axes[0].legend()
print(f"v(Km) / Vmax = {michaelis_menten_rate(Km, Vmax, Km) / Vmax:.3f}")

# %%
# The progress curve: substrate first disappears at a nearly constant
# rate (enzyme saturated), then decays exponentially once [S] falls
# below :math:`K_m`.

progress = MichaelisMentenProgress(S0=20.0, Vmax=Vmax, Km=Km)
result = progress.integrate((0.0, 4.0), dt=1e-3, method="rk4")
plot_concentration_vs_time(result, ax=axes[1])
axes[1].axhline(Km, color="crimson", linestyle=":", linewidth=0.8, label="[S] = Km")
axes[1].legend()
axes[1].set_title("Substrate progress curve")

# %%
# Inhibition: competitive inhibitors raise the apparent :math:`K_m` but
# can be out-competed by enough substrate; noncompetitive inhibitors
# lower :math:`V_{max}` itself.

axes[2].plot(S, michaelis_menten_rate(S, Vmax, Km), label="uninhibited", color="steelblue")
axes[2].plot(S, competitive_inhibition_rate(S, I=3.0, Vmax=Vmax, Km=Km, Ki=1.0), label="competitive (I=3)", color="darkorange")
axes[2].plot(S, noncompetitive_inhibition_rate(S, I=3.0, Vmax=Vmax, Km=Km, Ki=1.0), label="noncompetitive (I=3)", color="seagreen")
axes[2].axhline(Vmax, color="gray", linestyle="--", linewidth=0.8)
axes[2].set_xlabel("[S]")
axes[2].set_ylabel("v")
axes[2].set_title("Competitive vs. noncompetitive inhibition")
axes[2].legend()

fig.tight_layout()
plt.show()
