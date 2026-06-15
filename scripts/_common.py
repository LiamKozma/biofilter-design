"""Shared helpers for the driver scripts: config loading and path setup."""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)


def load_config(path: str) -> dict:
    cfg = yaml.safe_load(Path(path).read_text())
    return cfg


def operating_uL():
    """Return (u, length) for the bench column from the tidy operating file."""
    from biofilter import data

    op = data.load_operating()
    u = op["Q_m3_s"] / op["area_m2"]
    return u, op["bed_length_m"]
