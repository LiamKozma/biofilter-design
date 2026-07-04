# Project status & to-do

_Last updated: 2026-07-04 (all jobs finished; pipeline COMPLETE)_

## ✅ Status: all six stages complete

All SLURM jobs from the 2026-07-03 batch finished cleanly and the final SBC
verdict has been produced. The three deliverable files all exist:

- `results/calib_summary_sapelo2_full.json`  posterior medians + 95% CIs
- `results/design_uq_sapelo2_full.json`      credible interval on industrial volume
- `results/sbc_reduce_sapelo2_full.json`     SBC rank-uniformity verdict — **all 12 params pass** (p_uniform > 0.05)

Nothing is left running. What remains is analysis / write-up and committing the
new summary JSONs (see "Next steps" below).

## Pipeline stages

| Stage | Script / job | Output | Status |
|-------|--------------|--------|--------|
| A — kinetics model select | `01_kinetics_modelselect.py` | `results/stageA_kinetics.json` | ✅ done (Jun 15) |
| D — Sobol sensitivity | `slurm/sobol.sub` | `results/sobol_sapelo2_full.json` | ✅ done (Jun 15) |
| C — Bayesian calibration | `slurm/calibrate.sub` | `calib_chain_sapelo2_full.{npy,h5}` + summary | ✅ done (job 46682671; 40000 steps, ESS≈3397, acc=0.26, ~14 h) |
| E — design under uncertainty | `slurm/design_uq.sub` | `design_uq_sapelo2_full.json` | ✅ done (job 46682753) |
| F — SBC (400 datasets) | `slurm/sbc_resubmit.sub` | `sbc_sapelo2_full_task*.json` | ✅ done (job 46682544; 400/400 on disk) |
| G — SBC reduce (verdict) | `06_sbc_reduce.py` | `sbc_reduce_sapelo2_full.json` | ✅ done (2026-07-04; all params rank-uniform) |

## Key results & caveat

- **SBC passed** for all 12 parameters (chi-square p_uniform between 0.13 and 0.96)
  → the inference machinery is statistically calibrated, no detectable bias.
- **Design-UQ is honest but very wide.** Median required industrial volume
  ≈ **2579 m³** (vs the original report's 3500 m³ point estimate), but the 95%
  credible interval is **[1393, 698131] m³** and only **53.4 %** of posterior draws
  yield a *feasible* design. This traces directly to weakly-identified kinetics: the
  posterior on `log_Rmax` and `log_Ks` spans ~5 log units each (the classic
  Monod Rmax/Ks trade-off — the bench data mainly constrain the *ratio*, not each
  separately). This is a legitimate UQ finding, not a bug, but it's the headline
  caveat for any write-up.

## Next steps

1. **Commit the three new summary JSONs** (tracked deliverables; not gitignored)
   and this STATUS update, then push. Do git ops on the cluster login node:
   ```bash
   cd /scratch/lmk04992/biofilter-design
   git add STATUS.md results/calib_summary_sapelo2_full.json \
           results/design_uq_sapelo2_full.json results/sbc_reduce_sapelo2_full.json
   git commit -m "Add final calibration/design/SBC summaries; pipeline complete"
   git push
   ```
2. **Write-up / figures.** Regenerate posterior + SBC-rank + design-volume figures
   into `figures/` for the report (they're gitignored, so produce locally).
3. **Optional follow-up** if the design interval is too wide to be useful:
   re-parametrise to the identifiable `Rmax/Ks` ratio, add a mildly informative
   prior on `Ks` from literature, or register the infeasible-draw fraction as a
   design-risk number rather than widening the CI.

## Git / GitHub

- Repo pushed to **https://github.com/LiamKozma/biofilter-design** (branch `master`).
- `.claude/` scrubbed from history (amended commit `f82802d`) and gitignored.
- Working tree clean as of the push; `diag.sh` / `run_diagnostics*.sh` were committed —
  remove them later if you don't want the ad-hoc diagnostics in the repo.

## What was fixed (2026-07-03)

Earlier runs timed out because `effectiveness_factor` (`src/biofilter/biofilm.py`)
thrashed the biofilm BVP in the diffusion-limited region **and** silently returned a
~10×-too-low `η` there via a wrong `tanh(φ)/φ` fallback. Fixed with a hybrid: clamp
`s ≥ 0`, exact deep-film closed form for `φ ≥ 12`, BVP in between (~15× faster on the
η-table, 0% error vs a high-res reference). Also added HDF5 checkpoint/resume to
calibration (`scripts/02_calibrate.py`) and the 44-ID `slurm/sbc_resubmit.sub`.
See README "Effectiveness factor, three regimes".
