Examples
========

This gallery walks through every public feature of ``chemistrykit.spectro``:
the Beer-Lambert absorbance law and its stray-light deviation from
linearity; rigid-rotor rotational spectra with isotope shifts; harmonic
vs. Morse-potential vibrational band positions, plus a genuine triatomic
normal-mode calculation; Franck-Condon vibronic progressions for
electronic spectra; NMR from Larmor frequencies and chemical shifts to
J-coupling multiplets, the Karplus relation, and Fourier-transform NMR;
atomic line spectra (Fraunhofer, Kirchhoff-Bunsen, Balmer-Rydberg); and
spectral lineshapes.

Each script in this gallery is self-contained and can be run directly with
``python examples/spectro/<section>/<script>.py``. Every script also
carries an RST module docstring as its title/description and uses ``# %%``
markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source of
truth for what you see, not a copy of it.

Sections
--------

- **beer_lambert** -- the Beer-Lambert absorbance law and its stray-light
  deviation from linearity at high concentration.
- **atomic** -- Fraunhofer's dark solar lines, Kirchhoff and Bunsen's
  flame emission/absorption spectra, and the hydrogen series from the
  Balmer-Rydberg formula.
- **rotational** -- rigid-rotor rotational spectra, relative line
  intensities, and microwave isotope shifts.
- **vibrational** -- harmonic vs. Morse IR band positions and
  anharmonicity constants, Wilson GF-matrix normal modes for CO2 and H2O,
  Raman/IR mutual exclusion, and a linear molecule's degenerate bend.
- **electronic** -- Franck-Condon vibronic progressions in a UV-Vis
  absorption band.
- **nmr** -- Larmor frequencies and chemical shifts, first-order
  J-coupling multiplets, the Karplus relation, and Fourier-transform NMR.
- **lineshapes** -- Gaussian, Lorentzian, and Voigt line profiles.
