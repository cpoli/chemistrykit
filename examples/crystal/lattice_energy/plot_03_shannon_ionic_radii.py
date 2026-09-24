r"""
Shannon's effective ionic radii: additivity and periodic trends
=================================================================

Shannon (1976) fitted one self-consistent radius per ion (for a given
coordination number) so that cation-anion distances in thousands of
crystals are reproduced by simple sums :math:`r_++r_-`.
``SHANNON_IONIC_RADII_PM`` holds the 6-coordinate values; here their sums
are checked against measured rock-salt nearest-neighbor distances, and
their periodic trends are plotted.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.crystal.utils.reference_data import SHANNON_IONIC_RADII_PM as R

measured_r0 = {  # rock-salt nearest-neighbor distances, pm
    ("Li+", "F-"): 201.0,
    ("Na+", "F-"): 232.0,
    ("Na+", "Cl-"): 282.0,
    ("K+", "Cl-"): 315.0,
    ("K+", "Br-"): 330.0,
    ("Rb+", "I-"): 367.0,
    ("Mg2+", "O2-"): 211.0,
    ("Ca2+", "O2-"): 240.0,
}
sums, measured, labels = [], [], []
for (cat, an), r0 in measured_r0.items():
    s = R[cat] + R[an]
    sums.append(s)
    measured.append(r0)
    labels.append(f"{cat.rstrip('+2')}{an.rstrip('-2')}")
    print(f"{cat:5s}+{an:4s} r+ + r- = {s:6.1f} pm   measured r0 = {r0:6.1f} pm   ({(s - r0) / r0:+.1%})")
    assert abs(s - r0) / r0 < 0.05

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
alkali = ["Li+", "Na+", "K+", "Rb+", "Cs+"]
halide = ["F-", "Cl-", "Br-", "I-"]
axes[0].plot(range(len(alkali)), [R[i] for i in alkali], "o-", label="alkali cations (M+)")
axes[0].plot(range(len(halide)), [R[i] for i in halide], "s-", label="halide anions (X-)")
axes[0].set_xticks(range(5), ["period 2", "3", "4", "5", "6"])
axes[0].set_ylabel("Shannon radius, CN 6 (pm)")
axes[0].set_title("Radii grow down a group")
axes[0].legend()
axes[1].scatter(measured, sums)
for label, x, y in zip(labels, measured, sums, strict=True):
    axes[1].annotate(label, (x, y), textcoords="offset points", xytext=(5, -10))
axes[1].plot([190, 380], [190, 380], "k--", lw=0.8)
axes[1].set_xlabel("measured r0 (pm)")
axes[1].set_ylabel(r"Shannon $r_+ + r_-$ (pm)")
axes[1].set_title("Additivity of Shannon radii")
plt.tight_layout()
plt.show()
