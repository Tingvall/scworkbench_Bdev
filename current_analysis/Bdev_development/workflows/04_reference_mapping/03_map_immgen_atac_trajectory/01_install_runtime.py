#!/usr/bin/env python3
"""Reuse the qualified Step 4 installer for this isolated ImmGen run library."""

from pathlib import Path
import runpy


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "02_qualify_bdev_peakvi_pseudocells"
    / "01_install_runtime.py"
)
runpy.run_path(str(SCRIPT), run_name="__main__")
