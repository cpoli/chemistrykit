r"""
Kováts retention indices from an n-alkane ladder
===================================================

Under isothermal gas chromatography, :math:`\log t'` of the n-alkanes
rises linearly with carbon number. Kováts placed every other compound on
that ladder: :func:`~chemistrykit.analytical.kovats_retention_index`
interpolates :math:`\log t'` between the bracketing alkanes, giving an
index (n-octane = 800, n-nonane = 900, ...) that is far more
transferable between instruments than a raw retention time.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.analytical import kovats_retention_index

t0 = 1.20  # dead time, min
carbons = np.arange(6, 12)
# Illustrative isothermal alkane retention times following log t' = a + b*n.
t_alkanes = t0 + 0.035 * 2.05**carbons

# Unknowns eluting between the alkanes (retention times in min, chosen to
# give indices close to literature values on a nonpolar phase):
unknowns = {"toluene": 9.57, "ethylbenzene": 17.40, "1-octanol": 71.8}

for name, t_x in unknowns.items():
    n = int(np.searchsorted(t_alkanes, t_x)) - 1 + carbons[0]
    i = n - carbons[0]
    I = kovats_retention_index(t_x, t_alkanes[i], t_alkanes[i + 1], n=n, dead_time=t0)
    print(f"{name:13s}: t_R = {t_x:5.2f} min, between C{n} and C{n + 1} -> I = {I:.0f}")

# %%
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(100 * carbons, np.log10(t_alkanes - t0), "o-", color="gray", label="n-alkanes (I = 100 n)")
for name, t_x in unknowns.items():
    n = int(np.searchsorted(t_alkanes, t_x)) - 1 + carbons[0]
    i = n - carbons[0]
    I = kovats_retention_index(t_x, t_alkanes[i], t_alkanes[i + 1], n=n, dead_time=t0)
    ax.plot([I], [np.log10(t_x - t0)], "s", markersize=8, label=f"{name}: I = {I:.0f}")
ax.set_xlabel("retention index I")
ax.set_ylabel(r"$\log_{10} t'_R$ (adjusted retention time, min)")
ax.set_title("Kováts index: interpolation on the alkane ladder")
ax.legend()
plt.tight_layout()
plt.show()
