# Deep Research prompt — Bayesian biofilter calibration & design UQ

> Paste everything below the line into Gemini Deep Research. It is self-contained.
> When it returns, save the report as `deepresearch.txt` in this directory.

---

## Role and goal

You are helping a graduate researcher diagnose and improve a Bayesian
calibration + design-under-uncertainty study of a **gas-phase packed-bed
biofilter** treating odorous VOCs. I have a completed computational pipeline that
is statistically self-consistent (it passes Simulation-Based Calibration) but
produces an engineering answer with an implausibly wide uncertainty. I need a
literature-grounded diagnosis of **why**, and concrete, cited options to fix it.
Prioritize peer-reviewed biofiltration / bioreactor engineering and Bayesian
inverse-problem literature. Give numbers with units and citations, not
generalities.

## System being modeled

- **Process**: aerobic gas-phase biofiltration in a packed bed (media: compost
  and a proprietary biosorbent "BIOSORBENS", bulk density ~1447 kg/m³, BET
  surface area ~47.6 m²/g; compost ~730 kg/m³, ~3.6 m²/g).
- **Target compounds** (odorants): **3-methyl-1-butanol** (isoamyl alcohol),
  **2-methyl-1-butanol**, and **hexanal**. Bench inlet concentrations ~12–40 ppmv
  for the alcohols and ~1–7 ppmv for hexanal.
- **Bench column**: 0.1 m diameter, 0.5 m packed depth, air flow 4 L/min,
  empty-bed contact time ~59 s. A single steady-state run (Dec 2006) with gas
  concentration measured at the inlet plus ~4 axial ports. Removal is nearly
  complete: the alcohols and hexanal fall from inlet levels to ≈0 ppmv within the
  first 0.15–0.25 m of bed.

## Forward model I am calibrating

A **two-scale mechanistic model**:

1. **Biofilm scale** (reaction–diffusion BVP): substrate `S(z)` inside a biofilm
   of thickness `Lf` obeys `Df·S'' = Rmax·S/(Ks+S)`, no-flux at the support,
   bulk concentration at the gas–biofilm interface. This yields a Thiele modulus
   `φ = Lf·√(Rmax/(Df·Ks))` and an effectiveness factor `η ∈ (0,1]`.
2. **Column scale** (advection–dispersion–reaction BVP with Danckwerts BCs):
   `Dax·Cg'' − u·Cg' − a·Lf·η·(Cg/H)·Rmax·(Cg/H)/(Ks + Cg/H) = 0`, where `a` is
   specific biofilm interfacial area per bed volume and `H` a dimensionless
   (gas/liquid) Henry constant.

**Parameters** (all sampled in log space; shared across compounds unless noted):
`Df` (effective diffusivity, m²/s), `Lf` (film thickness, m), `a` (1/m),
`Dax` (axial dispersion, m²/s), `H` (dimensionless Henry), observation noise
`σ` (ppmv), and **per-compound** Monod `Rmax` (volumetric max rate) and `Ks`
(half-saturation).

My current **weakly-informative log-normal priors** (median, and ±1σ multiplicative
factor `e^σ`): `Df` 1e-9 m²/s (×2.7), `Lf` 1e-4 m (×2.7), `a` 2e3 1/m (×4.5),
`Dax` 1e-4 m²/s (×4.5), `H` 20 (×2.7), `Rmax` 2 (×7.4), `Ks` 5 (×7.4).

## What the calibration produced (the problem)

- Posteriors are **very wide and the kinetics look weakly identified**: per-compound
  `log Rmax` and `log Ks` each have 95% credible intervals spanning ~5 natural-log
  units (~150× range). `log H` also spans a ~28× range.
- **SBC passes** for all 12 parameters (rank-uniformity chi-square p = 0.13–0.96),
  so the inference is not *buggy* — the data genuinely fail to constrain the
  parameters individually.
- **Design step**: I size an industrial biofilter for **95% removal of hexanal at
  2 ppmv inlet**, air flow 2831.68 m³/min, under a length ≤ 3×diameter structural
  rule, by root-finding the bed length over 1–200 m. Result: median required
  volume ≈ 2580 m³, but the **95% credible interval is [1393, 698000] m³**, and
  only **53%** of posterior draws achieve the target at all (the rest cannot reach
  95% removal even in a 200 m bed). A global (Sobol) sensitivity analysis says
  removal is driven most by `a` and `H`, then `Rmax` and `Ks`, with `Dax` and
  `Lf` negligible.

## Research questions

### A. Root cause — identifiability
1. Is the **joint non-identifiability of Monod `Rmax` and `Ks` from steady-state
   packed-bed concentration-vs-depth profiles** a documented result? Specifically,
   when the local substrate concentration is well below `Ks` (pseudo-first-order
   regime), is it established that only the ratio `Rmax/Ks` (or an equivalent
   first-order rate) is identifiable, not the two separately? Cite structural /
   practical identifiability analyses for Monod kinetics in biofilters or
   biofilms.
2. Does **near-complete removal early in the bed** (data hitting ≈0 within the
   first third of the column) further destroy kinetic information — i.e., do the
   downstream zero-concentration points only impose a *lower bound* on activity?
   How do practitioners avoid this "over-performing bed" problem when designing
   kinetic experiments?
3. Which parameters are, in principle, **separately identifiable** from this kind
   of single-run steady-state profile data, and which are confounded (e.g.,
   `a·Lf`, `Rmax/Ks`, `η·Rmax`)?

### B. Literature values for informative priors
4. Published **Monod / Michaelis–Menten kinetic parameters** (`Rmax` or `Vmax` or
   elimination capacity, and `Ks`) for **gas-phase biofiltration of C4–C6
   alcohols and aldehydes** — ideally isoamyl alcohol / 2-/3-methylbutanol and
   hexanal, otherwise the closest analogues on compost/biosorbent media. Give
   values with units and phase (gas- vs liquid-phase concentration basis) and note
   typical **elimination capacities (g m⁻³ h⁻¹)** and critical loads.
5. Literature ranges for the **shared physical parameters** in VOC biofilters:
   effective diffusivity `Df` of these VOCs in wet biofilm (m²/s), biofilm
   thickness `Lf`, specific interfacial area `a` (m²/m³) for compost media, axial
   dispersion `Dax` / bed Peclet number, and the **dimensionless air/water Henry
   constant `H`** for 3-methyl-1-butanol, 2-methyl-1-butanol, and hexanal at
   ~20–25 °C.
6. **Henry's constant is tabulated physical chemistry** — should `H` be *fixed*
   from literature rather than inferred? What are the accepted dimensionless
   air/water Henry values (and their temperature dependence) for these three
   compounds, and how much would fixing `H` be expected to tighten the design?

### C. Model structure and effectiveness factor
7. Is the **two-scale (biofilm reaction–diffusion + column advection–dispersion)
   structure with Danckwerts BCs** the standard mechanistic framework for
   steady-state biofilter modeling (e.g., Ottengraf-type, Deshusses-type models)?
   What are its known limitations for these operating conditions?
8. Are there **validated effectiveness-factor / Thiele-modulus correlations for
   Monod (not purely first-order) kinetics** in a slab biofilm that I can
   benchmark my numerical `η(φ, S_bulk)` against, including the diffusion-limited
   (large-φ) asymptote?

### D. Design-under-uncertainty methodology
9. Is a design credible interval spanning **~3 orders of magnitude, with ~50% of
   posterior draws infeasible**, a recognized pathology when propagating a
   weakly-identified posterior to a sizing decision? How is this handled in the
   Bayesian engineering-design / robust-design literature?
10. What is **best practice for reporting and deciding** the size in this
    situation — e.g., quantile-based / percentile sizing, decision-theoretic
    cost-optimal sizing, equivalence to a conventional safety factor, or reporting
    **required EBRT and elimination capacity** rather than an absolute volume?
11. Is sizing on the **single hardest compound (hexanal, 2 ppmv, 95% removal)**
    the right target, versus sizing on total VOC mass loading or the
    mass-transfer-limited compound? What does standard biofilter practice use?
12. **Sanity-check the industrial duty** (2831.68 m³/min of air at ~2 ppmv
    hexanal, 95% removal): using conventional **EBRT-based design (typical EBRT
    15–60 s) and elimination-capacity charts**, what bed volume would a
    practitioner expect? Does that fall inside or outside my [1393, 698000] m³
    interval, and does it suggest the mechanistic root-find is mis-scaled?

### E. Experiments that would fix it
13. Ranked by expected reduction in design uncertainty, **what additional bench
    measurements** would most help: multiple inlet concentrations spanning `Ks`,
    transient step/loading tests, replicate steady-state runs, denser sampling
    near the inlet where the gradient is steep, or direct measurement of biofilm
    properties (`Lf`, `a`, biomass) and independent `H`/`Df`? Cite optimal
    experimental design work for biofilter kinetics if available.

## Requested output format

1. A short **diagnosis** (2–4 paragraphs) of the most likely root cause(s) of the
   wide design interval, tied to the specifics above.
2. A **parameter table**: for each of `Df, Lf, a, Dax, H, Rmax, Ks` (per relevant
   compound), give a literature central value + range + units + citation, and a
   recommended prior (or "fix at X").
3. A ranked list of **concrete next steps** (prior changes, reparametrization,
   design-metric changes, experiments), each with the expected effect and a
   citation.
4. Full **reference list** with DOIs/links.
