"""Two-scale Bayesian biofilter model.

Subpackages:
    biofilm   -- reaction-diffusion BVP -> effectiveness factor eta(Sb)
    column    -- advection-dispersion-reaction BVP for the packed bed
    simulator -- couples biofilm + column into a forward predictor
    kinetics  -- integral-method kinetic models (Stage A model selection)
    likelihood-- priors + likelihood for PDE-in-the-loop Bayesian calibration
    design    -- scale-up and design-under-uncertainty
"""
__all__ = ["biofilm", "column", "simulator", "kinetics", "likelihood", "design"]
