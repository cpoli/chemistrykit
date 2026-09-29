r"""
NRTL: a non-ideal liquid that splits into two phases
====================================================

Wilson's model cannot predict liquid-liquid immiscibility: its Gibbs
energy of mixing stays convex whatever the parameters. Renon and
Prausnitz's non-random two-liquid model
(:class:`~chemistrykit.thermo.NRTLSolution`) adds a non-randomness
parameter :math:`\alpha` and can. When

.. math::

    \frac{\Delta G_{mix}}{RT} = x_1\ln x_1 + x_2\ln x_2 + \frac{G^E}{RT}

develops two minima, a single liquid lowers its Gibbs energy by splitting
into two phases whose compositions share a common tangent. Below, NRTL
parameters typical of a partially miscible organic-water pair give such a
split, located by equating each component's activity
:math:`x_i\gamma_i` in the two phases.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from chemistrykit.constants import R
from chemistrykit.thermo import NRTLSolution, WilsonSolution

T = 298.15
nrtl = NRTLSolution(tau12=2.2, tau21=1.8, P1_star=1.0, P2_star=1.0, alpha=0.2)
wilson = WilsonSolution(Lambda12=0.05, Lambda21=0.05, P1_star=1.0, P2_star=1.0)

x = np.linspace(1e-4, 1 - 1e-4, 2000)


def g_mix(model, x1):
    return x1 * np.log(x1) + (1 - x1) * np.log(1 - x1) + model.excess_gibbs(x1, T) / (R * T)


def equal_activities(z):
    xa, xb = z
    g1a, g2a = nrtl.activity_coefficients(xa)
    g1b, g2b = nrtl.activity_coefficients(xb)
    return [np.log(xa * g1a) - np.log(xb * g1b), np.log((1 - xa) * g2a) - np.log((1 - xb) * g2b)]


xa, xb = fsolve(equal_activities, [0.02, 0.95])
assert xb - xa > 0.5, "expected two distinct liquid phases"
ga, gb = g_mix(nrtl, xa), g_mix(nrtl, xb)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
axes[0].plot(x, g_mix(nrtl, x), color="darkorange", label="NRTL")
axes[0].plot([xa, xb], [ga, gb], "k--", label="common tangent")
axes[0].plot([xa, xb], [ga, gb], "ko")
axes[0].set_xlabel("mole fraction of component 1")
axes[0].set_ylabel(r"$\Delta G_{mix}/RT$")
axes[0].set_title(f"NRTL: the liquid splits into x = {xa:.3f} and {xb:.3f}")
axes[0].legend()

axes[1].plot(x, g_mix(wilson, x), color="steelblue")
axes[1].set_xlabel("mole fraction of component 1")
axes[1].set_ylabel(r"$\Delta G_{mix}/RT$")
axes[1].set_title("Wilson, even strongly non-ideal: always one phase")
fig.tight_layout()

# %%
# A single phase is stable only where :math:`\Delta G_{mix}` is convex;
# NRTL's curvature changes sign, Wilson's never does.

print("NRTL curvature changes sign:", bool((np.diff(g_mix(nrtl, x), 2) < 0).any()))
print("Wilson curvature changes sign:", bool((np.diff(g_mix(wilson, x), 2) < 0).any()))

plt.show()
