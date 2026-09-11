r"""
Zero-, first-, and second-order integrated rate laws
======================================================

The three textbook elementary rate laws -- :class:`~chemistrykit.kinetics.systems.rate_laws.ZeroOrder`,
:class:`~chemistrykit.kinetics.systems.rate_laws.FirstOrder`, and
:class:`~chemistrykit.kinetics.systems.rate_laws.SecondOrder` -- all
started from the same initial concentration and with rate constants
chosen to give the *same* half-life, so the different curvature (linear,
exponential, hyperbolic) is the only thing distinguishing them. First
order is the only one whose half-life doesn't depend on the starting
concentration; the plot below marks all three half-lives to make that
concrete.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.systems.rate_laws import FirstOrder, SecondOrder, ZeroOrder

C0 = 1.0
t_half_target = 5.0

zero = ZeroOrder(k=C0 / (2.0 * t_half_target), C0=C0)
first = FirstOrder(k=np.log(2.0) / t_half_target, C0=C0)
second = SecondOrder(k=1.0 / (t_half_target * C0), C0=C0)

t = np.linspace(0.0, 20.0, 400)

fig, ax = plt.subplots(figsize=(7, 5))
for law, label, color in [(zero, "zero order", "steelblue"), (first, "first order", "darkorange"), (second, "second order", "seagreen")]:
    ax.plot(t, np.clip(law.concentration(t), 0.0, None), label=f"{label} (t_1/2={law.half_life():.2f})", color=color)
    ax.axvline(law.half_life(), color=color, linestyle=":", alpha=0.5)

ax.axhline(C0 / 2.0, color="gray", linestyle="--", linewidth=0.8, label="[A]_0 / 2")
ax.set_xlabel("t")
ax.set_ylabel("[A]")
ax.set_title("Same half-life, three different rate laws")
ax.legend()
fig.tight_layout()

# %%
# All three curves cross [A]_0/2 at their respective half-life by
# construction. The zero-order law is the only one that reaches exactly
# zero (and stays there -- the model doesn't allow negative
# concentration), while first- and second-order decay approach zero only
# asymptotically.

for law, label in [(zero, "zero"), (first, "first"), (second, "second")]:
    print(f"{label}-order: [A](t_1/2) = {law.concentration(law.half_life()):.6f} (expected {C0 / 2.0})")

plt.show()
