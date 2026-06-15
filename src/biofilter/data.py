"""Loaders for the tidy bench-scale data produced by scripts/00_extract_data.py."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TIDY = ROOT / "data" / "tidy"


def load_profiles() -> pd.DataFrame:
    df = pd.read_csv(TIDY / "profiles.csv")
    df["run_date"] = df["run_date"].astype(str)
    return df


def load_operating() -> dict:
    row = pd.read_csv(TIDY / "operating_conditions.csv").iloc[0]
    return row.to_dict()


def compound_arrays(df: pd.DataFrame, compound: str):
    """Return (t, y, run_idx, n_runs, C0_obs) for one compound across runs.

    ``t`` residence time (min), ``y`` concentration (ppmv), ``run_idx`` integer
    run label, ``C0_obs`` the inlet concentration of each run.
    """
    sub = df[df["compound"] == compound].copy()
    runs = sorted(sub["run_date"].unique())
    run_to_idx = {r: i for i, r in enumerate(runs)}
    sub["ridx"] = sub["run_date"].map(run_to_idx)
    sub = sub.sort_values(["ridx", "position_m"])
    t = sub["residence_time_min"].to_numpy()
    y = sub["conc_ppmv"].to_numpy()
    ridx = sub["ridx"].to_numpy()
    C0_obs = np.array(
        [sub[(sub["ridx"] == i) & (sub["position_m"] == 0.0)]["conc_ppmv"].iloc[0]
         for i in range(len(runs))]
    )
    return t, y, ridx, len(runs), C0_obs, runs
