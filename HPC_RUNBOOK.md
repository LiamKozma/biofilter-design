# Sapelo2 Runbook

Step-by-step to produce the real results on UGA's Sapelo2. Stages 0–1 are light
(run on a login node or interactively); stages C–F are the HPC workload.

## 0. Transfer the repo

From your laptop:

```bash
rsync -av --exclude .venv --exclude results --exclude figures \
  ~/Documents/biofilter-design/ \
  MyID@sapelo2.gacrc.uga.edu:/scratch/MyID/biofilter-design/
```

(Use `/scratch` for the run — it's the fast filesystem and where jobs should do I/O.)

## 1. One-time setup (login node)

```bash
cd /scratch/MyID/biofilter-design
bash slurm/setup_sapelo2.sh        # builds .venv incl. mpi4py against the MPI module
```

If module names differ, fix them with `module spider Python` / `module spider OpenMPI`
and edit the `module load` lines in `slurm/setup_sapelo2.sh` and the `*.sub` files.

## 2. Tidy data (login node, seconds)

```bash
source .venv/bin/activate
python scripts/00_extract_data.py
python scripts/01_kinetics_modelselect.py     # Stage A, light; writes results/stageA_kinetics.json
```

## 3. Submit the heavy stages

```bash
sbatch slurm/calibrate.sub      # Stage C  -> results/calib_chain_sapelo2_full.npy
sbatch slurm/sobol.sub          # Stage D  -> results/sobol_sapelo2_full.json
sbatch slurm/sbc_array.sub      # Stage F  -> results/sbc_sapelo2_full_task*.json
```

Watch them:

```bash
squeue --me
tail -f logs/calib_*.out
```

## 4. Post-process

```bash
# Stage E needs the calibration chain to be finished first:
python scripts/04_design_uq.py  configs/sapelo2_full.yaml   # -> results/design_uq_sapelo2_full.json
# Aggregate the SBC array:
python scripts/06_sbc_reduce.py configs/sapelo2_full.yaml   # -> results/sbc_reduce_sapelo2_full.json
```

## 5. Send results back

The whole `results/` directory is small (JSON + one `.npy` chain). Pull it home:

```bash
rsync -av MyID@sapelo2.gacrc.uga.edu:/scratch/MyID/biofilter-design/results/ \
  ~/Documents/biofilter-design/results/
```

Then hand me `results/` and I'll turn it into the figures and the website write-up.

## Expected wall-times (order of magnitude, full config)

| Stage | Job | Rough cost |
|-------|-----|-----------|
| C calibration | 1 node, 64 ranks | ~6–10 h |
| D Sobol | 1–2 nodes | ~10–30 min |
| F SBC | 400-task array, 50 concurrent | ~12–24 h wall |

If queue limits or allocation make 400 SBC datasets too much, drop `sbc.n_datasets`
(and the `--array` range in `sbc_array.sub`) to 100 — still a valid calibration check.

## What to send me when done

- `results/stageA_kinetics.json`
- `results/calib_chain_sapelo2_full.npy` + `results/calib_summary_sapelo2_full.json`
- `results/sobol_sapelo2_full.json`
- `results/design_uq_sapelo2_full.json`
- `results/sbc_reduce_sapelo2_full.json`
- the `logs/*.out` files (so I can sanity-check acceptance fractions / ESS / timings)
