r"""
The Langmuir adsorption isotherm
===================================

:class:`~chemistrykit.surface.systems.langmuir.LangmuirIsotherm` models
monolayer adsorption on a fixed pool of identical surface sites. Its
fractional coverage :math:`\theta=KP/(1+KP)` reaches exactly one-half at
the half-saturation pressure :math:`P=1/K`, and :func:`~chemistrykit.surface.systems.langmuir.fit_langmuir`
recovers `(K, qmax)` from noisy loading-vs-pressure data via the standard
1/q-vs-1/P linearization.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.surface.systems.langmuir import LangmuirIsotherm, fit_langmuir
from chemistrykit.surface.visualizers.surface_plots import plot_isotherm, plot_linearization

K_true, qmax_true = 2.5, 8.0
iso = LangmuirIsotherm(K=K_true, qmax=qmax_true)

P_half = iso.half_saturation_pressure()
theta_half = iso.fractional_coverage(P_half)
print(f"Half-saturation pressure P=1/K = {P_half:.4f}")
print(f"Coverage there: {theta_half:.6f} (should be exactly 0.5)")

# %%
# Generate synthetic noisy data and recover the parameters by fitting the
# linearized isotherm.
rng = np.random.default_rng(0)
P_data = np.array([0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0])
q_data = iso.loading(P_data) * (1.0 + rng.normal(scale=0.01, size=P_data.shape))

fit = fit_langmuir(P_data, q_data)
print(f"\nTrue (K, qmax)   = ({K_true}, {qmax_true})")
print(f"Fitted (K, qmax) = ({fit.K:.4f}, {fit.qmax:.4f})")
print(f"R^2 of linearized fit: {fit.r_squared:.6f}")

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
plot_isotherm(iso, P_max=20.0, P_data=P_data, q_data=q_data, ax=axes[0])
axes[0].axvline(P_half, color="gray", linestyle="--", linewidth=0.8)
axes[0].axhline(iso.qmax * theta_half, color="gray", linestyle="--", linewidth=0.8)

x_lin = 1.0 / P_data
y_lin = 1.0 / q_data
plot_linearization(x_lin, y_lin, fit=(1.0 / (fit.qmax * fit.K), 1.0 / fit.qmax), ax=axes[1], xlabel="1/P", ylabel="1/q")
plt.tight_layout()
plt.show()
