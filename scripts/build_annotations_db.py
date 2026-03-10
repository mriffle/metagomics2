#!/usr/bin/env python3
"""Thin wrapper — delegates to the installed package module.

Usage:
    python scripts/build_annotations_db.py \
        --fasta uniprot_sprot.fasta.gz \
        --gaf goa_uniprot_all.gaf.gz \
        --output uniprot_sprot.annotations.db

Or via the CLI entry point:
    metagomics2-build-annotations \
        --fasta uniprot_sprot.fasta.gz \
        --gaf goa_uniprot_all.gaf.gz \
        --output uniprot_sprot.annotations.db
"""

import sys
from pathlib import Path

# Add src to path so we can import metagomics2 modules
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from metagomics2.scripts.build_annotations_db import main

if __name__ == "__main__":
    main()
