"""Sphinx configuration for smi-tiled documentation.

Build locally:
    pixi run docs-build
    pixi run docs-serve          # http://localhost:8765

Hosted at:
    https://smi-tiled.readthedocs.io     (via .readthedocs.yaml)
    https://EliotGann.github.io/smi-tiled (via .github/workflows/docs.yml)
"""
from __future__ import annotations

import os
import sys
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version as _pkg_version

# Make the source available for autodoc.
sys.path.insert(0, os.path.abspath("../../src"))

# ---------------------------------------------------------------------------
# Project information
# ---------------------------------------------------------------------------
project = "smi-tiled"
author = "NSLS-II SMI / Contributors"
copyright = f"{datetime.now().year}, {author}"

try:
    release = _pkg_version("smi-tiled")
except PackageNotFoundError:
    release = "0.0.0+unknown"
version = ".".join(release.split(".")[:2])

# ---------------------------------------------------------------------------
# Sphinx extensions
# ---------------------------------------------------------------------------
extensions = [
    # Auto-generate API docs from docstrings
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",            # NumPy/Google docstring styles
    "sphinx.ext.intersphinx",         # link to numpy, scipy, etc.
    "sphinx.ext.viewcode",            # "[source]" links
    "sphinx.ext.todo",
    # Rich content / better rendering
    "myst_parser",                    # Markdown alongside reStructuredText
    "sphinx_autodoc_typehints",       # nicer type signatures
    "sphinx_copybutton",              # copy button on code blocks
    "sphinx_design",                  # cards, tabs, grids
    "sphinxcontrib.mermaid",          # diagrams
]

# Auto-generate stub pages for each member listed in :autosummary:.
autosummary_generate = True
autosummary_imported_members = False

# autodoc: pull in __init__ docstrings, show inheritance.
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"      # show types in the body, not signature
autodoc_typehints_format = "short"

# Napoleon: parse both NumPy and Google docstrings; keep the docstrings
# we already have working.
napoleon_numpy_docstring = True
napoleon_google_docstring = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_examples = True
napoleon_use_rtype = False             # render Returns inline

# MyST: enable a sensible default set of Markdown extensions.
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 3
myst_url_schemes = ("http", "https", "mailto", "ftp")

# Source-file extensions Sphinx recognizes.
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Cross-reference into upstream docs.
intersphinx_mapping = {
    "python":  ("https://docs.python.org/3", None),
    "numpy":   ("https://numpy.org/doc/stable/", None),
    "scipy":   ("https://docs.scipy.org/doc/scipy/", None),
    "xarray":  ("https://docs.xarray.dev/en/stable/", None),
    "pandas":  ("https://pandas.pydata.org/docs/", None),
    "skimage": ("https://scikit-image.org/docs/stable/", None),
    "tiled":   ("https://blueskyproject.io/tiled/", None),
}

# Show :todo: directives in development; suppress in production.
todo_include_todos = bool(int(os.environ.get("SMI_TILED_DOCS_TODO", "0")))

templates_path = ["_templates"]

# Files/dirs to ignore during build.
exclude_patterns = [
    "_build", "Thumbs.db", ".DS_Store",
]

# Suppress noisy warnings that are cosmetic:
#  - autosectionlabel.*: duplicate section labels across pages
#  - ref.python:         missing intersphinx refs for third-party typing
#  - docutils:           napoleon-emitted asterisks in numpy "Attributes"
#                        sections occasionally trip docutils inline-emphasis
#                        detection; rendering is correct.
suppress_warnings = [
    "autosectionlabel.*",
    "ref.python",
    "docutils",
]

# ---------------------------------------------------------------------------
# HTML theme
# ---------------------------------------------------------------------------
html_theme = "pydata_sphinx_theme"
html_title = "smi-tiled"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sphinx = False

html_theme_options = {
    "github_url": "https://github.com/EliotGann/smi-tiled",
    "use_edit_page_button": True,
    "show_toc_level": 2,
    "navbar_align": "left",
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
    "secondary_sidebar_items": ["page-toc", "edit-this-page", "sourcelink"],
    "header_links_before_dropdown": 6,
    "icon_links": [
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/smi-tiled/",
            "icon": "fa-brands fa-python",
            "type": "fontawesome",
        },
    ],
}

html_context = {
    "github_user": "EliotGann",
    "github_repo": "smi-tiled",
    "github_version": "main",
    "doc_path": "docs/source",
}

# ---------------------------------------------------------------------------
# Code-block styling
# ---------------------------------------------------------------------------
pygments_style = "sphinx"
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

# Mermaid (for the mask-architecture diagram in user-guide/masks.md).
mermaid_version = "10.6.1"
