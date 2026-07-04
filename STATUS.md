# Project status & to-do

_Last updated: 2026-07-04 (v2 reanalysis + comparison figures complete)_

## ✅ Where things stand — nothing running, resume-ready

Both v1 (full two-scale Monod, HPC) and v2 (identifiable first-order, local) are
complete. v2 resolves the v1 design-UQ pathology (500× credible interval, 47%
infeasible → ~2.3×, all feasible). Details in the two sections below.

**To finish: commit the v2 deliverables + figures (git ops run cluster-side).**

```bash
cd /scratch/lmk04992/biofilter-design
git add README.md STATUS.md .gitignore configs/first_order.yaml \
        src/biofilter/firstorder.py \
        scripts/10_calibrate_fo.py scripts/11_design_fo.py \
        scripts/12_sbc_fo.py scripts/13_figures_v1v2.py \
        results/calib_summary_first_order.json results/design_uq_first_order.json \
        results/sbc_reduce_first_order.json \
        figures/fig_identifiability.png figures/fig_design_volume.png \
        figures/fig_fit.png figures/fig_sbc.png \
        deepresearch_prompt.md deepresearchresults.txt
git commit -m "v2: identifiable first-order reanalysis + comparison figures"
git push
```

**Optional (not blocking), for the write-up:**
- Longer-inner-chain SBC (8–10k steps) to clean up the mild `sigma` non-uniformity
  (nuisance only — all design-relevant `k1` are already calibrated; see SBC note).
- The thesis/report prose itself (the numbers, figures, and diagnosis are ready).
- Recommended lab experiments to separate `Rmax`/`Ks` (higher inlet load / shorter
  EBCT / transient step tests) — Discussion-section material, not compute.

**How to regenerate v2 end-to-end** (laptop, ~a few min; chain `.npy` is gitignored):
```bash
python scripts/10_calibrate_fo.py configs/first_order.yaml   # k1 posteriors
python scripts/11_design_fo.py     configs/first_order.yaml   # design volume + CI
python scripts/12_sbc_fo.py        configs/first_order.yaml   # SBC (~45 min, 400 sets)
python scripts/13_figures_v1v2.py  configs/first_order.yaml   # the 4 figures
```

---

## 🔬 v2 reanalysis — identifiable pseudo-first-order model

Deep Research (`deepresearchresults.txt`) confirmed the v1 diagnosis: the wide
design interval is **structural non-identifiability**, not a bug (SBC was right).
The bench inlet concentrations sit far below `Ks`, so the column runs entirely in
the **pseudo-first-order regime** where only the lumped ratio `Rmax/Ks` — an
effective first-order rate — is identifiable, never `Rmax` and `Ks` separately.
Free `H` and the `Lf`↔`a` product added more confounding; the "over-performing
bed" (near-complete removal in the first third) left the slow-kinetics tail
unbounded → the 47%-infeasible, 698,000 m³ pathology.

**v2 fix** (`src/biofilter/firstorder.py`, `scripts/10–12`, `configs/first_order.yaml`):
infer one identifiable first-order rate `k1` (1/s) + noise `sigma` per compound;
fix `H`, `Df` to literature; drop `Dax` (plug flow, Pe≫100); analytic forward
model; quantile-based design. The forward model is analytic, so **v2 runs locally
in minutes — no cluster** (the HPC machinery is what you'd need only for the full
two-scale model, *if* you had data rich enough to identify it).

| v2 stage | Script | Output | Result |
|-------|--------|--------|--------|
| C — calibration | `10_calibrate_fo.py` | `calib_summary_first_order.json` | ✅ acc=0.50, ESS≈4900. k₁: 3-MB 0.068 [.058,.078], 2-MB 0.064 [.056,.071], Hexanal 0.058 [.041,.092] 1/s |
| E — design UQ | `11_design_fo.py` | `design_uq_first_order.json` | ✅ V median **2421 m³**, 95% CI **[1527, 3493]** (~2.3×, was ~500×); all draws feasible; conservative P95 = 3277 m³ (L=33.5, D=11.2 m) |
| F+G — SBC | `12_sbc_fo.py` | `sbc_reduce_first_order.json` | ✅ all three `k1` calibrated; `sigma` (nuisance) mildly off. See SBC note. |

**SBC note (400 datasets, thinned).** The first pass had a clip-at-0 in the SBC
*simulator* (`firstorder.simulate`) that didn't match the plain-Gaussian
likelihood — a bug in the *test*, not the inference; fixed (simulate un-clipped).
Final verdict (400 datasets, chains thinned ×12):

- **All three rate constants `k1` are rank-uniform** (p = 0.81 / 0.11 / 0.55),
  including the design-critical `k1[Hexanal]` (p=0.55). The design rests entirely
  on these, so **the design UQ is calibrated.**
- The **nuisance noise scales** `sigma[3-MB]`/`sigma[2-MB]` show mild residual
  non-uniformity (p ≈ 0.003 / 0.006); `sigma[Hexanal]` is fine (p=0.46). It
  survived 400 datasets + thinning, so it is a small *real* effect, not MC noise
  or autocorrelation. Verified **not a bug** (likelihood `-N logσ`, log-space
  prior, and simulator are mutually consistent). Leading cause: the alcohols'
  steeper decay gives a more correlated `k1`–`sigma` posterior that 2500 inner
  steps slightly under-explore. It does **not** affect any reported quantity.
- Optional refinement if a spotless SBC is wanted for the write-up: longer inner
  chains (e.g. 8–10k steps) or an independent-draw (integer-rank) SBC for the
  scale parameters. Not blocking.

The conservative P95 design (3277 m³, L=33.5 m, D=11.2 m) essentially reproduces
the original report's point estimate (~3500 m³, L=34.23 m, D=11.41 m) — but now
with a principled credible interval instead of a bare number. The median (2421 m³)
matches a conventional EBRT sanity check (~50 s EBRT → ~2360 m³).

**Only way to ever separate `Rmax` from `Ks`** = new experiments that leave the
starvation regime: higher inlet load (approach/exceed `Ks`), shorter EBCT, or
transient step-loading. That's a Discussion recommendation, not more compute.

---

## ✅ v1 status: all six stages complete (retained for provenance)

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

## Key v1 results & caveat (superseded by v2)

- **SBC passed** for all 12 parameters (chi-square p_uniform between 0.13 and 0.96)
  → the inference machinery is statistically calibrated, no detectable bias.
- **Design-UQ is honest but very wide.** Median required industrial volume
  ≈ **2579 m³** (vs the original report's 3500 m³ point estimate), but the 95%
  credible interval is **[1393, 698131] m³** and only **53.4 %** of posterior draws
  yield a *feasible* design. This traces directly to weakly-identified kinetics: the
  posterior on `log_Rmax` and `log_Ks` spans ~5 log units each (the classic
  Monod Rmax/Ks trade-off — the bench data mainly constrain the *ratio*, not each
  separately). **v2 (top of file) resolves this** by inferring only the identifiable
  first-order rate; the v1 wide interval and the diagnosis are kept here as the
  motivating result.

_(v1 summary JSONs were already committed on 2026-07-04. The remaining commit is
the v2 deliverables + figures — see the resume checklist at the top.)_

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
