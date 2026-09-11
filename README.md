# chemistrykit

[![CI](https://github.com/chemistrykit/chemistrykit/actions/workflows/ci.yml/badge.svg)](https://github.com/chemistrykit/chemistrykit/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/chemistrykit/badge/?version=latest)](https://chemistrykit.readthedocs.io)
[![PyPI](https://img.shields.io/pypi/v/chemistrykit.svg)](https://pypi.org/project/chemistrykit/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/chemistrykit.svg)](pyproject.toml)

A unified numerical toolkit for computational chemistry: pure-numerical
models across reaction kinetics, thermodynamics, quantum chemistry,
electrochemistry, spectroscopy, and more -- no cheminformatics
dependencies (no RDKit/ASE/PySCF/OpenMM), sharing common ODE integrators,
chemical constants, and a consistent NumPy-based API. Conventionally
imported as `ck`. chemistrykit clones the architecture and engineering
conventions of the sibling project
[physicskit](https://github.com/physicskit/physicskit) (`pk`).

This is an early, in-progress build. Only `chemistrykit.kinetics` (reaction
kinetics) is implemented so far, as the reference domain the remaining 13
domains will follow -- see `chemistrykit-spec.md` for the full plan. It
covers elementary integrated rate laws, the Arrhenius equation,
Michaelis-Menten enzyme kinetics, a general stoichiometric
reaction-network engine, and the Brusselator oscillator.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Quick start

```python
import chemistrykit as ck
import numpy as np

# Closed-form first-order decay and its half-life
rate_law = ck.kinetics.FirstOrder(k=0.1, C0=1.0)
print(rate_law.half_life())  # ln(2) / k

# A general stoichiometric reaction network: A -> B -> C
network = ck.kinetics.StoichiometricNetwork.consecutive(k1=1.0, k2=0.3, A0=1.0)
result = network.integrate((0.0, 10.0), dt=1e-3, method="rk4")
print(result.concentration("B")[-1])
```

## Subpackages

Domain subpackages, each with runnable examples linked below:

- [`chemistrykit.kinetics`](https://chemistrykit.readthedocs.io/en/latest/examples/kinetics.html) -- reaction kinetics: integrated rate laws and half-lives, the Arrhenius equation and activation-energy fitting, Michaelis-Menten enzyme kinetics with Lineweaver-Burk linearization and inhibition, a general stoichiometry-matrix reaction-network engine (parallel/consecutive/reversible/steady-state-approximation chains), and the Brusselator oscillating reaction network.

Shared infrastructure, used across the subpackages above rather than
standalone toolkits:

- `chemistrykit.constants` -- chemical constants (R, NA, k_B, Faraday's
  constant, ...) from `scipy.constants`, plus a small built-in
  periodic-table data table and a few well-defined unit conversions.
- `chemistrykit.integrators` -- shared numerical ODE integrators (RK4,
  leapfrog, Yoshida4, adaptive Dormand-Prince) used across the other
  subpackages.

## Test

Tests live alongside each subpackage, at `chemistrykit/<name>/tests/`.

```bash
pytest                                              # everything
pytest chemistrykit/kinetics/tests                  # a single subpackage

# docstring examples, across every subpackage:
MPLBACKEND=Agg pytest --doctest-modules chemistrykit \
    --ignore-glob="*/tests/*"
```

Both commands, plus `ruff check`/`ruff format --check`, run in CI on
every PR (`.github/workflows/ci.yml`) across Python 3.9-3.12 on Linux and
macOS. See [CONTRIBUTING.md](CONTRIBUTING.md) before opening a PR.

## Docs

Built docs are hosted at <https://chemistrykit.readthedocs.io>. To build
locally:

```bash
pip install -e ".[docs]"
cd docs && make html
```

## Citation

If you use chemistrykit in your research, please cite it — see
[CITATION.cff](CITATION.cff).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Please note that this project
follows the [Contributor Covenant](CODE_OF_CONDUCT.md).

## License

MIT -- see [LICENSE](LICENSE).
