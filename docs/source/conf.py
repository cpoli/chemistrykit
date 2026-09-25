"""Sphinx configuration for chemistrykit."""

import os
import sys
from pathlib import Path

import matplotlib

sys.path.insert(0, os.path.abspath("../.."))


# Some example animations can embed just over Matplotlib's default 20MB
# HTML5-video limit. Sphinx-Gallery's built-in "matplotlib" reset_modules
# entry calls plt.rcdefaults() before *every* example script runs, which
# would immediately undo a plain module-level rcParams assignment here, so
# the raised limit has to be reapplied as a reset_modules callable of its
# own, ordered right after "matplotlib" (ported from physicskit's conf.py,
# which hit this with a large particle-scatter animation).
def _raise_animation_embed_limit(gallery_conf, fname):
    matplotlib.rcParams["animation.embed_limit"] = 50


project = "chemistrykit"
copyright = "2026, chemistrykit contributors"
author = "chemistrykit team"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # Supports NumPy-style docstrings
    "sphinx.ext.mathjax",
    "sphinx.ext.doctest",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
    "myst_parser",
    "sphinx_gallery.gen_gallery",
    "sphinx_design",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = False

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
autodoc_typehints = "description"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Single source of truth for the 14 chemistry subpackages. Every card grid
# (homepage, api/index, examples/index, history/index, and the per-subpackage
# hub pages) and the cross-link strip atop each api/examples/history page is
# generated from this table by _generate_subpackage_docs() below, instead of
# being hand-duplicated across five-plus RST files that used to drift out of
# sync with each other (different card ordering, stale blurbs, etc.) --
# ported structure-for-structure from physicskit's conf.py.
SUBPACKAGES = [
    {
        "name": "analytical",
        "category": "Analytical Methods",
        "blurb": "Redox/EDTA titration-curve simulation, chromatographic plate theory, and calibration curves.",
    },
    {
        "name": "crystal",
        "category": "Structure & Materials",
        "blurb": "Crystal systems and unit-cell volume, hard-sphere packing, ionic lattice energy, and powder XRD.",
    },
    {
        "name": "electrochem",
        "category": "Simulation & Electrochemistry",
        "blurb": "The Nernst equation, standard reduction potentials, Butler-Volmer kinetics, and battery discharge.",
    },
    {
        "name": "kinetics",
        "category": "Reaction Dynamics",
        "blurb": "Integrated rate laws, the Arrhenius equation, Michaelis-Menten enzyme kinetics, and reaction-network dynamics.",
    },
    {
        "name": "md",
        "category": "Simulation & Electrochemistry",
        "blurb": "The Lennard-Jones fluid, bonded potentials, and velocity-rescaling/Nose-Hoover thermostats.",
    },
    {
        "name": "photochem",
        "category": "Reaction Dynamics",
        "blurb": "Jablonski-diagram excited-state kinetics, quantum yields, Stern-Volmer quenching, and photostationary states.",
    },
    {
        "name": "polymer",
        "category": "Structure & Materials",
        "blurb": "Chain statistics and the Flory exponent, molecular-weight distributions, and step-/chain-growth kinetics.",
    },
    {
        "name": "quantum",
        "category": "Quantum & Spectroscopy",
        "blurb": "Particle-in-a-box, the harmonic oscillator and Morse potential, the rigid rotor, and Huckel theory.",
    },
    {
        "name": "solutions",
        "category": "Thermodynamics & Equilibrium",
        "blurb": "pH/pOH and acid-base equilibria, buffers, titration curves, and solubility/activity coefficients.",
    },
    {
        "name": "spectro",
        "category": "Quantum & Spectroscopy",
        "blurb": "Beer-Lambert absorbance, rotational/vibrational spectra, Franck-Condon progressions, and NMR multiplets.",
    },
    {
        "name": "statmech",
        "category": "Thermodynamics & Equilibrium",
        "blurb": "Partition functions and their thermodynamic functions, the Maxwell-Boltzmann distribution, and lattice-gas adsorption.",
    },
    {
        "name": "structure",
        "category": "Structure & Materials",
        "blurb": "VSEPR geometry prediction, point-group determination, bond order, and formal-charge assignment.",
    },
    {
        "name": "surface",
        "category": "Reaction Dynamics",
        "blurb": "Langmuir/Freundlich/BET adsorption isotherms and Langmuir-Hinshelwood surface-reaction kinetics.",
    },
    {
        "name": "thermo",
        "category": "Thermodynamics & Equilibrium",
        "blurb": "Equations of state, phase boundaries, reaction equilibrium, and Raoult's/Henry's law mixtures.",
    },
]

for _s in SUBPACKAGES:
    _s.setdefault("history_doc", f"{_s['name']}_breakthroughs")
    # Every subpackage's "Examples" link goes straight to its sphinx-gallery
    # index (:orphan: like every such index -- meant to be linked to
    # directly rather than placed in a toctree). Where a narrative tutorial
    # also exists, that gallery's examples/<name>/README.rst header links
    # out to it, so there's no separate examples/<name>.rst stub page whose
    # only job is forwarding to one or the other.
    _s["examples_doc"] = f"api/gallery/{_s['name']}/index"
del _s

_CATEGORY_ORDER = [
    "Reaction Dynamics",
    "Thermodynamics & Equilibrium",
    "Quantum & Spectroscopy",
    "Structure & Materials",
    "Simulation & Electrochemistry",
    "Analytical Methods",
]

# Subpackages with a sphinx-gallery-formatted examples/<name>/ directory.
# All 14 domains from chemistrykit-spec.md are now implemented.
_GALLERY_SUBPACKAGES = [s["name"] for s in SUBPACKAGES]

sphinx_gallery_conf = {
    "examples_dirs": [f"../../examples/{name}" for name in _GALLERY_SUBPACKAGES],
    "gallery_dirs": [f"api/gallery/{name}" for name in _GALLERY_SUBPACKAGES],
    "filename_pattern": r"/plot_",
    "download_all_examples": False,
    "within_subsection_order": "FileNameSortKey",
    "remove_config_comments": True,
    "matplotlib_animations": True,
    "reset_modules": ("matplotlib", "seaborn", _raise_animation_embed_limit),
    # Lets ".. minigallery::" (used throughout docs/source/history/) resolve
    # fully-qualified object names in addition to the file paths/globs it
    # already handles, and is required for it to not warn about falling
    # back to file-path resolution on every single invocation.
    "backreferences_dir": "gen_modules/backreferences",
    "doc_module": ("chemistrykit",),
}

# As more subpackages land, common attribute names (e.g. "dim", "species")
# become genuinely ambiguous across unrelated classes; autolinking to a
# single target for them isn't meaningful, so treat those as non-fatal
# rather than broken links (same rationale as physicskit's conf.py).
suppress_warnings = ["ref.python", "config.cache"]

html_theme = "pydata_sphinx_theme"
html_logo = "_static/images/chemistrykit_logo_transparent.png"
html_theme_options = {
    "github_url": "https://github.com/cpoli/chemistrykit",
    "icon_links": [
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/chemistrykit/",
            "icon": "fa-brands fa-python",
            "type": "fontawesome",
        },
    ],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "show_toc_level": 2,
    "navigation_with_keys": True,
    "navigation_depth": 2,
    "logo": {
        "alt_text": "chemistrykit logo",
    },
}
html_static_path = ["_static"]


def _card(link, blurb, link_title):
    """One sphinx-design grid-item-card, indented for direct concatenation into a ``.. grid::`` block."""
    return f"   .. grid-item-card:: {link_title}\n      :link: {link}\n      :link-type: doc\n\n      {blurb}\n\n"


def _grid(cards):
    return ".. grid:: 1 2 3 3\n   :gutter: 2\n\n" + "".join(cards)


def _grouped_grid(link_fn):
    """A ``.. grid::`` per category (in _CATEGORY_ORDER), each preceded by a
    rubric heading, with ``link_fn(subpackage)`` giving each card's target.
    Shared by the homepage hub grid and the api/examples/history grids so
    they show the same category grouping instead of a flat alphabetical
    list.
    """
    parts = []
    for category in _CATEGORY_ORDER:
        members = [s for s in SUBPACKAGES if s["category"] == category]
        if not members:
            continue
        parts.append(f".. rubric:: {category}\n\n")
        parts.append(_grid([_card(link_fn(s), s["blurb"], f"chemistrykit.{s['name']}") for s in members]))
        parts.append("\n")
    return "".join(parts)


def _write(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _generate_subpackage_docs(app):
    """Generate the card grids, per-subpackage hub pages, and cross-link
    strips derived from SUBPACKAGES.

    Runs at "builder-inited", the same point sphinx-gallery uses to write
    its own generated RST, so the output exists before Sphinx reads any
    source file that ``.. include::``/``.. toctree::``\\ s it. Everything
    lands under _generated/ (gitignored, like api/gallery/) rather than
    being committed, since it's entirely derived from SUBPACKAGES above --
    there's nothing hand-authored in it to keep in version control.
    """
    out = Path(app.srcdir) / "_generated"

    # RST substitutions for numbers/facts derived from SUBPACKAGES, so prose
    # elsewhere (e.g. the homepage's "N domains" pitch) doesn't hardcode a
    # count that silently goes stale the next time a subpackage is added.
    _write(out / "vars.rst", f".. |num_subpackages| replace:: {len(SUBPACKAGES)}\n")

    # Category-grouped grids reused by index.rst, api/index.rst,
    # examples/index.rst, and history/index.rst, so those five-plus listings
    # of the same 14 subpackages -- and their grouping -- come from one
    # place instead of being hand-copied (and drifting) independently.
    _write(out / "grid_api.rst", _grouped_grid(lambda s: f"/api/{s['name']}"))
    _write(out / "grid_examples.rst", _grouped_grid(lambda s: f"/{s['examples_doc']}"))
    _write(out / "grid_history.rst", _grouped_grid(lambda s: f"/history/{s['history_doc']}"))

    # Domain-first homepage grid: same grouping, linking to each
    # subpackage's own hub page rather than straight into a single
    # History/Examples/API silo.
    _write(out / "grid_hub.rst", _grouped_grid(lambda s: f"/_generated/subpackages/{s['name']}"))

    # Per-subpackage hub page: History, Examples, and API reference as
    # three equally-weighted doors into the same domain, so a reader who
    # wants "everything about kinetics" doesn't have to browse three
    # separate site-wide sections to find kinetics in each of them.
    hub_dir = out / "subpackages"
    for s in SUBPACKAGES:
        name = s["name"]
        title = f"chemistrykit.{name}"
        cards = _grid(
            [
                _card(
                    f"/history/{s['history_doc']}",
                    "The foundational breakthroughs behind this subpackage, linked to the implementation.",
                    "History",
                ),
                _card(
                    f"/{s['examples_doc']}",
                    "Runnable tutorials and the full example gallery.",
                    "Examples",
                ),
                _card(
                    f"/api/{name}",
                    "Every public class and function.",
                    "API reference",
                ),
            ]
        )
        _write(hub_dir / f"{name}.rst", f"{title}\n{'=' * len(title)}\n\n{s['blurb']}\n\n{cards}")

    # Cross-link strip included atop each hand-authored api/<name>.rst and
    # history/<name>_breakthroughs.rst page, so a reader who lands on just
    # one of them (e.g. from a search engine) can discover the others for
    # the same subpackage.
    nav_dir = out / "nav"
    for s in SUBPACKAGES:
        name = s["name"]
        links = [
            f":doc:`chemistrykit.{name} hub </_generated/subpackages/{name}>`",
            f":doc:`History </history/{s['history_doc']}>`",
            f":doc:`Examples </{s['examples_doc']}>`",
            f":doc:`API reference </api/{name}>`",
        ]
        _write(nav_dir / f"{name}.rst", f".. container:: subpkg-nav\n\n   {' · '.join(links)}\n")


def setup(app):
    app.connect("builder-inited", _generate_subpackage_docs)
