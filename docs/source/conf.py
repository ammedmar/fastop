"""Sphinx configuration for the generated fastop documentation."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fastop import __version__

project = "fastop"
author = "Federico Cantero-Morán and Anibal M. Medina-Mardones"
copyright = "2026, Federico Cantero-Morán and Anibal M. Medina-Mardones"
release = __version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]

exclude_patterns = ["_build"]
html_theme = "alabaster"
