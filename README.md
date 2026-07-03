# Two-Scale Bayesian Biofilter Design

A graduate-level re-analysis of a bench-scale biofiltration experiment, rebuilt
as a reproducible computational pipeline. It takes the original BCHE 3420 design
project (raw concentration profiles from a packed-bed biofilter) and replaces its
back-of-the-envelope kinetics and deterministic scale-up with:

1. **Bayesian reaction-order selection** on the real data (hierarchical model,
   PSIS-LOO comparison) instead of "pick the highest R²".
2. A **mechanistic two-scale forward model** — a biofilm reaction–diffusion BVP
   coupled to a column advection–dispersion–reaction BVP — that resolves the
   diffusion limitation a lumped rate constant hides.
3. **PDE-in-the-loop Bayesian calibration** of that model to the bench data.
4. **Global (Sobol) sensitivity analysis** to find which parameters control
   removal.
5. **Posterior-predictive design under uncertainty**: a credible interval on the
   required industrial biofilter volume, not a single number.
6. **Simulation-Based Calibration (SBC)** to prove the inference is statistically
   calibrated.

Stages 3–6 embed a coupled BVP solve inside every likelihood/sample evaluation,
so the workload is genuinely HPC-scale and is designed to run on UGA's **Sapelo2**
cluster (SLURM scripts in `slurm/`).

> **Provenance.** The bench-scale measurements and the industrial duty spec are
> from the original undergraduate project. Everything in the modelling/inference
> pipeline here is a subsequent, more rigorous re-analysis — not a claim about
> what the original report contained.

## Layout

```
data/raw/        original spreadsheet + report PDF
data/tidy/       extracted, documented CSVs (scripts/00)
src/biofilter/   biofilm BVP, column BVP, simulator, kinetics, likelihood, design
scripts/         00 extract · 01 kinetics · 02 calibrate · 03 sobol
                 04 design-UQ · 05 SBC · 06 SBC-reduce
configs/         local_smoke.yaml (laptop) · sapelo2_full.yaml (cluster)
slurm/           Sapelo2 batch scripts: calibrate · sobol · design_uq ·
                 sbc_array (full) · sbc_resubmit (only the timed-out tasks)
results/         JSON summaries + posterior chains
```

## Quick start (laptop smoke test)

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt            # mpi4py optional locally
python scripts/00_extract_data.py
python scripts/01_kinetics_modelselect.py
python scripts/02_calibrate.py  configs/local_smoke.yaml
python scripts/03_sobol.py      configs/local_smoke.yaml
python scripts/04_design_uq.py  configs/local_smoke.yaml
python scripts/05_sbc.py        configs/local_smoke.yaml
python scripts/06_sbc_reduce.py configs/local_smoke.yaml
```

## Production run (Sapelo2)

```bash
# build an env that includes mpi4py against the cluster MPI module
module load Python/3.11 OpenMPI
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt     # includes h5py (resume) + tqdm (progress)

# Stage C (MPI). Writes results/calib_chain_sapelo2_full.h5 as it runs, so a
# walltime kill can be RESUMED by simply resubmitting this same script.
JC=$(sbatch --parsable slurm/calibrate.sub)

# Stage E depends on the calibration chain -- chain it so it auto-starts on
# clean completion instead of erroring on a missing .npy.
sbatch --dependency=afterok:$JC slurm/design_uq.sub

# Stage D (Sobol) and Stage F (SBC) are independent of C/E and can run anytime.
sbatch slurm/sobol.sub              # Stage D, MPI
sbatch slurm/sbc_array.sub          # Stage F, full 400-task array

# After the SBC array finishes, aggregate the rank files (light; login node OK):
python scripts/06_sbc_reduce.py configs/sapelo2_full.yaml
```

**Which stage runs where.** `calibrate`, `design_uq`, `sobol`, `sbc_array`
(and `sbc_resubmit`) all do forward solves and/or use MPI, so they are `sbatch`
batch jobs — never run them as bare `python` on a login node. Only
`06_sbc_reduce.py` is light enough (it just tallies rank files) for the login
node.

**Resuming a partial SBC.** Each SBC task is fully determined by its seed
(`default_rng(seed+1000+task)`), so if some array tasks time out you can rerun
*only* those without redoing the rest -- put their IDs in
`slurm/sbc_resubmit.sub`'s `--array=` line and `sbatch` it. `06_sbc_reduce.py`
globs whatever task files are present, so it picks up the backfilled ones
automatically.

## The model in one paragraph

Inside the biofilm, substrate `S(z)` obeys `Df S'' = Rmax S/(Ks+S)` with no-flux
at the support and bulk concentration at the interface; non-dimensionalising
gives the Thiele modulus `φ = Lf·√(Rmax/(Df·Ks))` and an effectiveness factor
`η ∈ (0,1]` (the fraction of biomass actually fed). At column scale the gas
concentration obeys `Dax Cg'' − u Cg' − a·Lf·η(Cg/H)·Rmax(Cg/H)/(Ks+Cg/H) = 0`
with Danckwerts boundary conditions. One forward solve = the η-table BVP family
plus the column BVP; calibration runs ~10⁶ of them, SBC ~10⁹.

**Effectiveness factor, three regimes.** `effectiveness_factor` avoids solving
the stiff biofilm BVP where it is unnecessary or unreliable: the first-order
form `tanh(φ)/φ` at tiny `sb`; the exact deep-film closed form
`η = (1+sb)/sb·√(2(sb−ln(1+sb)))/φ` for `φ ≥ 12` (diffusion-limited, validated
<0.1% against the solved BVP and correct at large `sb` where `tanh(φ)/φ`
under-predicts several-fold); and the solved BVP in between, with the RHS clamped
to the physical branch `s ≥ 0` so the solver cannot wander into the `1+s` pole.
This is both faster (~15× on the η-table) and more accurate than solving the BVP
everywhere — the naive version silently returned an order-of-magnitude-low `η`
in the diffusion-limited corner where the solver failed to converge.
