r"""
Sabatier and Senderens: catalytic hydrogenation over finely divided nickel
============================================================================

Sabatier and Senderens showed from 1897 that finely divided nickel
hydrogenates ethylene (:math:`C_2H_4 + H_2 \to C_2H_6`) at temperatures
where the gas-phase reaction does not run at all. This example uses
illustrative Arrhenius parameters and
:func:`~chemistrykit.surface.systems.catalysis.compare_catalyzed_rate` to
compute the conversion reached in a one-second contact time as the
temperature rises. Over nickel the reaction lights off hundreds of
kelvin earlier. A second panel shows Sabatier's rule of thumb from
screening many metals: a metal that binds the reactant too weakly holds
almost none of it on the surface, and one that binds it too strongly
leaves no free sites. Only a moderate binder has both.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.surface.systems.catalysis import compare_catalyzed_rate
from chemistrykit.surface.systems.langmuir import langmuir_coverage

tau = 1.0  # contact time, s
T = np.linspace(250.0, 900.0, 400)
comparisons = [compare_catalyzed_rate(Ea_uncatalyzed=180e3, Ea_catalyzed=45e3, T=Ti, A_uncatalyzed=1e13, A_catalyzed=1e7) for Ti in T]
X_uncat = 1.0 - np.exp(-np.array([c.k_uncatalyzed for c in comparisons]) * tau)
X_ni = 1.0 - np.exp(-np.array([c.k_catalyzed for c in comparisons]) * tau)

T_half_ni = T[np.searchsorted(X_ni, 0.5)]
T_half_gas = T[np.searchsorted(X_uncat, 0.5)]
print(f"50% conversion over Ni at   ~{T_half_ni:.0f} K")
print(f"50% conversion in gas phase ~{T_half_gas:.0f} K")
print(f"Nickel lowers the working temperature by ~{T_half_gas - T_half_ni:.0f} K")

# %%
# Sabatier's rule: surface coverage of the reactant (Langmuir) for three
# hypothetical metals at the same pressure.
P = 1.0
metals = {"weak binder": 0.02, "nickel-like": 1.0, "strong binder": 50.0}
theta = {name: langmuir_coverage(K, P) for name, K in metals.items()}
for name, th in theta.items():
    print(f"{name:14s} reactant coverage {th:.3f}, free sites {1 - th:.3f}")

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(T, X_uncat, label="gas phase, no catalyst")
axes[0].plot(T, X_ni, label="over finely divided Ni")
axes[0].set_xlabel("temperature (K)")
axes[0].set_ylabel("conversion of ethylene in 1 s")
axes[0].set_title("Catalytic hydrogenation lights off far earlier")
axes[0].legend()

names = list(metals)
x = np.arange(len(names))
axes[1].bar(x - 0.2, [theta[n] for n in names], width=0.4, label="covered by reactant")
axes[1].bar(x + 0.2, [1 - theta[n] for n in names], width=0.4, label="free sites")
axes[1].set_xticks(x, names)
axes[1].set_ylabel("fraction of surface sites")
axes[1].set_title("Sabatier's rule: bind neither too weakly nor too strongly")
axes[1].legend()
plt.tight_layout()
plt.show()
