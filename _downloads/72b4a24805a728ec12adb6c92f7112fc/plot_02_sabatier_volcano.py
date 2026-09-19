r"""
The Sabatier principle and the catalytic volcano curve
==========================================================

Paul Sabatier's qualitative principle -- that the best catalyst binds a
reaction intermediate neither too weakly (nothing sticks long enough to
react) nor too strongly (the surface clogs with product and never frees
up for another turnover) -- is illustrated here with the standard modern
textbook toy model (see, e.g., Chorkendorff & Niemantsverdriet, *Concepts
of Modern Catalysis and Kinetics*, 2nd ed., Wiley-VCH, 2007, Ch. 8): a
two-step surface mechanism whose overall rate needs *both* an occupied
site (to hold the adsorbed intermediate) and an empty site (for the next
step to proceed), so the rate is proportional to
:math:`\theta(1-\theta)`, with :math:`\theta` the Langmuir coverage
(:func:`~chemistrykit.surface.systems.langmuir.langmuir_coverage`) set by
the adsorption equilibrium constant `K` -- a proxy for binding strength.
This is not a literal reproduction of Sabatier's own (qualitative) or
Balandin's (multiplet-theory) equations, just the simplest possible
model with the same volcano-shaped structure they established
conceptually and structurally, respectively.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.surface.systems.langmuir import LangmuirIsotherm, langmuir_coverage

# Fixed working pressure; K (a stand-in for binding strength/free energy of
# adsorption) is scanned over many orders of magnitude on both sides of the
# Langmuir half-saturation point K = 1/P.
P = 1.0
K_values = np.logspace(-3, 3, 4000)
theta = langmuir_coverage(K_values, P)
rate = theta * (1.0 - theta)

peak_index = int(np.argmax(rate))
K_peak = K_values[peak_index]
theta_peak = theta[peak_index]

print(f"Rate-maximizing K: {K_peak:.4f} (theory: 1/P = {1.0 / P:.4f})")
print(f"Coverage there: {theta_peak:.4f} (theory: exactly 0.5)")

# %%
# The peak sits exactly at Langmuir's own half-saturation condition
# K = 1/P, theta = 1/2 -- weak binding (K << 1/P, left leg) starves the
# surface of the adsorbed intermediate; strong binding (K >> 1/P, right
# leg) starves it of empty sites for the intermediate to react further.
# Both legs suppress the rate; only the intermediate binding strength in
# between maximizes it.
iso_at_peak = LangmuirIsotherm(K=K_peak, qmax=1.0)
print(f"\nLangmuirIsotherm(K=K_peak).half_saturation_pressure() = {iso_at_peak.half_saturation_pressure():.4f} (should equal P = {P})")

weak_binding_rate = rate[np.searchsorted(K_values, 0.01)]
strong_binding_rate = rate[np.searchsorted(K_values, 100.0)]
print(f"\nRate at K=0.01 (weak binding):   {weak_binding_rate:.4f}")
print(f"Rate at K_peak (intermediate):    {rate[peak_index]:.4f}")
print(f"Rate at K=100 (strong binding):   {strong_binding_rate:.4f}")

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

axes[0].semilogx(K_values, theta)
axes[0].axvline(K_peak, color="gray", linestyle="--", linewidth=0.8)
axes[0].axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
axes[0].set_xlabel("K (binding-strength proxy)")
axes[0].set_ylabel(r"coverage $\theta$")
axes[0].set_title("Langmuir coverage vs. binding strength")

axes[1].semilogx(K_values, rate)
axes[1].axvline(K_peak, color="crimson", linestyle="--", linewidth=0.8, label=f"peak at K={K_peak:.2f}")
axes[1].set_xlabel("K (binding-strength proxy)")
axes[1].set_ylabel(r"rate $\propto \theta(1-\theta)$")
axes[1].set_title("The Sabatier volcano curve")
axes[1].legend()
plt.tight_layout()
plt.show()
