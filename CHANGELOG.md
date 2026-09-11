# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Initial repository scaffold: `pyproject.toml`, shared
  `chemistrykit.integrators` (RK4, leapfrog/velocity-Verlet, Yoshida4,
  adaptive Dormand-Prince, ported from physicskit), `chemistrykit.constants`
  (chemical constants from `scipy.constants` plus a built-in periodic-table
  data table), CI/pre-commit/readthedocs config, and root docs.
- `chemistrykit.kinetics` — the first complete domain subpackage, and the
  template the remaining 13 domains will follow: elementary integrated
  rate laws and half-lives (zero/first/second order); the Arrhenius
  equation and activation-energy fitting; Michaelis-Menten enzyme
  kinetics with Lineweaver-Burk linearization and competitive/
  noncompetitive inhibition; a general stoichiometry-matrix
  reaction-network engine (parallel, consecutive, reversible, and
  steady-state-approximation chains) integrated via
  `chemistrykit.integrators`; and the Brusselator oscillating reaction
  network. Includes a matching `examples/kinetics/` sphinx-gallery and
  `docs/source/api/kinetics.rst` / `docs/source/examples/kinetics.rst`.

[Unreleased]: https://github.com/chemistrykit/chemistrykit/compare/main...HEAD
