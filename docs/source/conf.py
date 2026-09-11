"""Sphinx configuration for chemistrykit."""

import os
import sys

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

# Subpackages with a sphinx-gallery-formatted examples/<name>/ directory.
# Add each remaining chemistry domain from chemistrykit-spec.md to this
# list as it lands.
_GALLERY_SUBPACKAGES: list[str] = [
    "kinetics",
    "thermo",
    "solutions",
    "md",
    "statmech",
]

sphinx_gallery_conf = {
    "examples_dirs": [f"../../examples/{name}" for name in _GALLERY_SUBPACKAGES],
    "gallery_dirs": [f"api/gallery/{name}" for name in _GALLERY_SUBPACKAGES],
    "filename_pattern": r"/plot_",
    "download_all_examples": False,
    "within_subsection_order": "FileNameSortKey",
    "remove_config_comments": True,
    "matplotlib_animations": True,
    "reset_modules": ("matplotlib", "seaborn", _raise_animation_embed_limit),
}

# As more subpackages land, common attribute names (e.g. "dim", "species")
# become genuinely ambiguous across unrelated classes; autolinking to a
# single target for them isn't meaningful, so treat those as non-fatal
# rather than broken links (same rationale as physicskit's conf.py).
suppress_warnings = ["ref.python"]

html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "github_url": "https://github.com/chemistrykit/chemistrykit",
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "show_toc_level": 2,
    "navigation_with_keys": True,
    "navigation_depth": 2,
    "logo": {
        "text": "chemistrykit",
    },
}
html_static_path = ["_static"]
