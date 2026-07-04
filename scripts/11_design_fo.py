#!/usr/bin/env python3
"""v2 Stage E: design under uncertainty from the first-order posterior.

For each posterior draw of the design compound's rate ``k1`` the required
industrial volume is analytic:  V = Q * EBRT,  EBRT = -ln(1-target)/k1.  Because
``k1`` is now tightly identified and bounded away from zero, the volume
distribution has no infinite tail and every draw is feasible -- the pathology of
the v1 root-finder (47% infeasible, 698,000 m^3 upper bound) is gone.

Reported as a decision-ready number: the median size plus a conservative
quantile (size on the slow-kinetics tail of ``k1``).

    python scripts/11_design_fo.py configs/first_order.yaml
"""
from __future__ import annotations

import json
import sys

import numpy as np

from _common import RESULTS, load_config
from biofilter import firstorder as FO


def main(cfg_path):
    cfg = load_config(cfg_path)
    name = cfg["run_name"]
    dcfg = cfg.get("design", {})
    target = dcfg.get("target_removal", FO.TARGET_REMOVAL)
    comp = dcfg.get("target_compound", "Hexanal")
    q = dcfg.get("conservative_quantile", 0.95)

    chain = np.load(RESULTS / f"calib_chain_{name}.npy")
    ci = cfg["compounds"].index(comp)
    k1 = np.exp(chain[:, ci])          # posterior draws of the design rate

    V = FO.required_volume(k1, target)          # vectorised over draws
    ebrt = FO.required_ebrt(k1, target)
    # Conservative design: slow-kinetics tail of k1 -> high-volume tail.
    k1_cons = np.percentile(k1, 100 * (1 - q))
    V_cons = float(FO.required_volume(k1_cons, target))
    L_cons, D_cons = FO.geometry_from_volume(V_cons)
    L_med, D_med = FO.geometry_from_volume(float(np.median(V)))

    out = {
        "run_name": name,
        "model": "pseudo_first_order",
        "target_compound": comp,
        "target_removal": target,
        "n_draws": int(k1.size),
        "frac_feasible": 1.0,                    # analytic: always feasible
        "k1_1_per_s": {
            "median": float(np.median(k1)),
            "ci95": [float(np.percentile(k1, 2.5)), float(np.percentile(k1, 97.5))],
        },
        "ebrt_s": {
            "median": float(np.median(ebrt)),
            "ci95": [float(np.percentile(ebrt, 2.5)), float(np.percentile(ebrt, 97.5))],
        },
        "volume_m3": {
            "median": float(np.median(V)),
            "ci95": [float(np.percentile(V, 2.5)), float(np.percentile(V, 97.5))],
            "ci50": [float(np.percentile(V, 25)), float(np.percentile(V, 75))],
        },
        "conservative_design": {
            "quantile": q,
            "volume_m3": V_cons,
            "length_m": float(L_cons),
            "diameter_m": float(D_cons),
        },
        "length_m_median": float(L_med),
        "diameter_m_median": float(D_med),
        "report_point_estimate_m3": 3500.0,
    }
    (RESULTS / f"design_uq_{name}.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out["volume_m3"], indent=2))
    print(f"conservative (P{int(q*100)}) volume: {V_cons:.0f} m3 "
          f"(L={L_cons:.1f} m, D={D_cons:.1f} m)")
    print(f"all draws feasible -> results/design_uq_{name}.json")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "configs/first_order.yaml")
