# Project status & to-do

_Last updated: 2026-07-03_

## Pipeline stages

| Stage | Script / job | Output | Status |
|-------|--------------|--------|--------|
| A — kinetics model select | `01_kinetics_modelselect.py` | `results/stageA_kinetics.json` | ✅ done (Jun 15) |
| C — Bayesian calibration | `slurm/calibrate.sub` | `results/calib_chain_sapelo2_full.{npy,h5}` + summary | ⏳ running (job **46682671**) |
| D — Sobol sensitivity | `slurm/sobol.sub` | `results/sobol_sapelo2_full.json` | ✅ done (Jun 15) |
| E — design under uncertainty | `slurm/design_uq.sub` | `results/design_uq_sapelo2_full.json` | ⏳ queued (job **46682753**, waits on 46682671) |
| F — SBC (400 datasets) | `slurm/sbc_array.sub` → resubmit | `results/sbc_sapelo2_full_task*.json` | ⏳ 356/400 done; 44 rerunning (job **46682544**) |
| G — SBC reduce (verdict) | `06_sbc_reduce.py` | `results/sbc_reduce_sapelo2_full.json` | ⛔ TODO — run after F completes |

## To-do (Liam, on the cluster)

- [ ] **Wait for the SBC resubmit `46682544` to finish**, then run the aggregation
      (light — login node is fine):
      ```bash
      python scripts/06_sbc_reduce.py configs/sapelo2_full.yaml
      ```
- [ ] **Confirm calibration `46682671` finished cleanly** and produced
      `results/calib_chain_sapelo2_full.npy` + `calib_summary_sapelo2_full.json`.
      - If it hit the walltime, just `sbatch slurm/calibrate.sub` again — it
        **resumes** from `calib_chain_sapelo2_full.h5` (needs h5py, already installed).
- [ ] **Confirm design-UQ `46682753` ran** (it releases automatically once
      calibration succeeds) and wrote `results/design_uq_sapelo2_full.json`.
      - If calibration failed, the dependency cancels it — resubmit after fixing:
        `sbatch --dependency=afterok:<new_calib_id> slurm/design_uq.sub`
- [ ] **Final check**: all three deliverables present —
      - posterior: `calib_summary_sapelo2_full.json` (medians + 95% CIs)
      - design: `design_uq_sapelo2_full.json` (credible interval on industrial volume)
      - SBC verdict: `sbc_reduce_sapelo2_full.json` (rank uniformity)

## Handy commands

```bash
squeue --me                                   # what's running/queued (note: --me)
sacct -j 46682671 --format=JobID,State,Elapsed,ExitCode   # calibration outcome
ls -la results/ | grep -v sbc_.*_task         # non-SBC outputs
ls results/sbc_sapelo2_full_task*.json | wc -l  # should reach 400
```

## What was fixed (2026-07-03)

The earlier runs timed out because `effectiveness_factor` (`src/biofilter/biofilm.py`)
thrashed the biofilm BVP in the diffusion-limited region **and** silently returned
a ~10×-too-low `η` there via a wrong `tanh(φ)/φ` fallback. Fixed with a hybrid:
clamp `s ≥ 0`, use the exact deep-film closed form for `φ ≥ 12`, BVP in between
(~15× faster, 0% error vs a high-res reference). Also added HDF5 checkpoint/resume
to calibration and a 44-ID SBC resubmit script. See README "Effectiveness factor,
three regimes".
