Examples
========

This gallery walks through every public feature of ``chemistrykit.spectro``:
the Beer-Lambert absorbance law and its stray-light deviation from
linearity; rigid-rotor rotational spectra with isotope shifts; harmonic
vs. Morse-potential vibrational band positions, plus a genuine triatomic
normal-mode calculation; Franck-Condon vibronic progressions for
electronic spectra; and a first-order NMR multiplet simulator.

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
- **rotational** -- rigid-rotor rotational spectra, relative line
  intensities, and the H/D isotope shift.
- **vibrational** -- harmonic vs. Morse IR band positions and
  anharmonicity constants, and a Wilson GF-matrix normal-mode calculation
  for CO2 and H2O.
- **electronic** -- Franck-Condon vibronic progressions in a UV-Vis
  absorption band.
- **nmr** -- first-order NMR multiplet simulation, from simple n+1
  multiplets to a genuine doublet of triplets.
