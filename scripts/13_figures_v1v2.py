#!/usr/bin/env python3
"""Comparison figures for the v1 (full Monod) vs v2 (identifiable first-order)
reanalysis.  Saves PNGs to figures/ (gitignored; regenerate locally).

    python scripts/13_figures_v1v2.py

Panels produced:
  fig_identifiability.png  v1 Rmax-Ks ridge (non-identifiable) vs v2 k1 (tight)
  fig_design_volume.png    required industrial volume: v1 500x CI vs v2 ~2.3x CI
  fig_fit.png              first-order model vs bench data (C/Cin collapse)
  fig_sbc.png              v2 SBC rank histograms (calibration check)
"""
from __future__ import annotations

import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _common import RESULTS, FIGURES, load_config
from biofilter import data, firstorder as FO

# Okabe-Ito colorblind-safe palette; fixed hue per compound, used everywhere.
C_3MB, C_2MB, C_HEX = "#0072B2", "#E69F00", "#D55E00"
COMP_COLOR = {"3-MB": C_3MB, "2-MB": C_2MB, "Hexanal": C_HEX}
C_V1, C_V2 = "#999999", "#009E73"          # old (recessive) vs new (headline)

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150,
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
    "legend.frameon": False,
})


def _v2_k1_draws(cfg):
    chain = np.load(RESULTS / "calib_chain_first_order.npy")
    return {c: np.exp(chain[:, i]) for i, c in enumerate(cfg["compounds"])}


def fig_identifiability(cfg):
    """v1 hexanal Rmax-Ks joint posterior (ridge) vs v2 k1 marginals (tight)."""
    v1 = np.load(RESULTS / "calib_chain_sapelo2_full.npy", mmap_mode="r")
    # hexanal: log_Rmax=col10, log_Ks=col11 (natural log). Convert to log10.
    rmax = v1[:, 10] / np.log(10.0)
    ks = v1[:, 11] / np.log(10.0)
    k1 = _v2_k1_draws(cfg)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.2, 4.0), constrained_layout=True)

    axL.hexbin(rmax, ks, gridsize=45, cmap="Greys", mincnt=1, linewidths=0)
    axL.set_title("v1 full Monod — Hexanal kinetics\n(collinear ridge: non-identifiable)")
    axL.set_xlabel(r"$\log_{10} R_{max}$")
    axL.set_ylabel(r"$\log_{10} K_s$")
    axL.text(0.05, 0.95, "each spans ~5 log units,\nbut only the ratio is fixed",
             transform=axL.transAxes, va="top", fontsize=8.5, color="0.3")
    axL.grid(False)

    xhi = max(0.15, float(np.percentile(np.concatenate(list(k1.values())), 99.8)) * 1.15)
    bins = np.linspace(0, xhi, 61)
    for c in cfg["compounds"]:
        axR.hist(k1[c], bins=bins, histtype="stepfilled", alpha=0.55,
                 color=COMP_COLOR[c], label=c, density=True)
        axR.hist(k1[c], bins=bins, histtype="step", color=COMP_COLOR[c], lw=1.4,
                 density=True)
    axR.set_xlim(0, xhi)
    axR.set_title("v2 first order — identified rate\n(tight, physically meaningful)")
    axR.set_xlabel(r"first-order rate $k_1$  (1/s)")
    axR.set_ylabel("posterior density")
    axR.legend(title=None)
    fig.savefig(FIGURES / "fig_identifiability.png")
    plt.close(fig)
    print("  fig_identifiability.png")


def fig_design_volume(cfg):
    """Required industrial volume: v1 500x CI vs v2 ~2.3x CI (log axis)."""
    k1_hex = _v2_k1_draws(cfg)["Hexanal"]
    V2 = FO.required_volume(k1_hex)                       # per-draw volumes
    v2med = float(np.median(V2))
    v2lo, v2hi = np.percentile(V2, [2.5, 97.5])
    v1 = json.loads((RESULTS / "design_uq_sapelo2_full.json").read_text())
    v1med = v1["volume_m3"]["median"]; v1lo, v1hi = v1["volume_m3"]["ci95"]
    dfo = json.loads((RESULTS / "design_uq_first_order.json").read_text())
    p95 = dfo["conservative_design"]["volume_m3"]

    fig, ax = plt.subplots(figsize=(8.6, 4.2), constrained_layout=True)
    ax.hist(V2, bins=np.logspace(np.log10(V2.min()*0.9), np.log10(V2.max()*1.1), 70),
            color=C_V2, alpha=0.75, label="v2 posterior (first order)")
    ax.axvspan(v2lo, v2hi, color=C_V2, alpha=0.12)

    # v1 95% CI as a wide error bar high on the axis (its distribution is not saved
    # and 47% of draws were infeasible, so we show its reported span honestly).
    y1 = ax.get_ylim()[1] * 0.86
    ax.plot([v1lo, v1hi], [y1, y1], color=C_V1, lw=3, solid_capstyle="round",
            label="v1 95% CI (full Monod)")
    ax.plot(v1med, y1, "o", color=C_V1, ms=7)
    ax.annotate("v1: 500× wide,\n47% infeasible", (v1hi, y1), color="0.4",
                fontsize=8.5, ha="right", va="bottom")

    ax.axvline(3500, color="0.2", ls="--", lw=1.2)
    ax.text(3500, ax.get_ylim()[1]*0.5, " original report\n 3500 m³", fontsize=8.5, color="0.2")
    ax.axvline(p95, color=C_HEX, ls=":", lw=1.6)
    ax.text(p95, ax.get_ylim()[1]*0.28, f" v2 conservative\n P95 {p95:.0f} m³",
            fontsize=8.5, color=C_HEX)

    ax.set_xscale("log")
    ax.set_xlabel("required industrial biofilter volume  (m³, log scale)")
    ax.set_ylabel("posterior draws")
    ax.set_title(f"Design under uncertainty — v2 median {v2med:.0f} m³, "
                 f"95% CI [{v2lo:.0f}, {v2hi:.0f}]")
    ax.legend(loc="upper left")
    fig.savefig(FIGURES / "fig_design_volume.png")
    plt.close(fig)
    print("  fig_design_volume.png")


def fig_fit(cfg):
    """First-order model vs bench data, C/Cin collapsed over runs."""
    df = data.load_profiles()
    fo = FO.build_fo_data(df, cfg["compounds"])
    k1 = _v2_k1_draws(cfg)
    tgrid = np.linspace(0, 60, 100)

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True,
                             constrained_layout=True)
    for ax, c in zip(axes, cfg["compounds"]):
        col = COMP_COLOR[c]
        for Cin, t_s, y in fo.per_compound[c]:
            ax.scatter(t_s, np.asarray(y) / Cin, s=34, color=col,
                       edgecolor="white", linewidth=0.6, zorder=3)
        klo, kmed, khi = np.percentile(k1[c], [2.5, 50, 97.5])
        ax.fill_between(tgrid, np.exp(-khi*tgrid), np.exp(-klo*tgrid),
                        color=col, alpha=0.18, zorder=1)
        ax.plot(tgrid, np.exp(-kmed*tgrid), color=col, lw=2, zorder=2)
        ax.set_title(f"{c}\n$k_1$={kmed:.3f} [{klo:.3f}, {khi:.3f}] 1/s", fontsize=10)
        ax.set_xlabel("residence time (s)")
    axes[0].set_ylabel(r"$C/C_{in}$")
    axes[0].set_ylim(-0.05, 1.05)
    fig.suptitle("First-order fit vs bench data (all three runs, inlet-normalised)",
                 fontsize=11)
    fig.savefig(FIGURES / "fig_fit.png")
    plt.close(fig)
    print("  fig_fit.png")


def fig_sbc(cfg):
    """v2 SBC rank histograms with the expected-uniform band."""
    red = json.loads((RESULTS / "sbc_reduce_first_order.json").read_text())
    pp = red["per_param"]
    n = red["n_datasets"]; nb = red["n_bins"]
    exp = n / nb
    band = 2.0 * np.sqrt(exp * (1 - 1/nb))    # ~95% Poisson-ish band around uniform

    fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.6), constrained_layout=True)
    for ax, (nm, d) in zip(axes.ravel(), pp.items()):
        counts = d["counts"]
        ok = d["p_uniform"] > 0.05
        color = C_V2 if ok else C_HEX
        ax.bar(np.arange(nb), counts, width=0.92, color=color, alpha=0.75)
        ax.axhline(exp, color="0.3", lw=1)
        ax.axhspan(exp - band, exp + band, color="0.5", alpha=0.15)
        ax.set_title(f"{nm}\np={d['p_uniform']:.3f} {'ok' if ok else 'CHECK'}",
                     fontsize=9.5, color="0.15" if ok else C_HEX)
        ax.set_xticks([])
        ax.set_ylim(0, max(counts) * 1.25)
    for ax in axes[:, 0]:
        ax.set_ylabel("count")
    fig.suptitle(f"v2 SBC rank uniformity — {n} datasets "
                 "(flat = calibrated)", fontsize=11)
    fig.savefig(FIGURES / "fig_sbc.png")
    plt.close(fig)
    print("  fig_sbc.png")


def main(cfg_path):
    cfg = load_config(cfg_path)
    print("writing figures to figures/ :")
    fig_identifiability(cfg)
    fig_design_volume(cfg)
    fig_fit(cfg)
    fig_sbc(cfg)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "configs/first_order.yaml")
