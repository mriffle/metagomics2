"""Metagomics 2 - Metaproteomics annotation and aggregation tool."""

import os

# Version can be overridden by METAGOMICS_VERSION environment variable (set in Docker)
__version__ = os.getenv("METAGOMICS_VERSION", "0.1.0")
