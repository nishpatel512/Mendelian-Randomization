import pandas as pd
import numpy as np
from scipy.stats import t,norm,chi2,stats
from scipy.optimize import minimize
import scipy.stats
import statsmodels.api as sm
from harmonise_data import *
from extract_instruments import *
from extract_outcome_data import *
from mr_weighted_median import *

def mr_penalised_weighted_median(b_exp, b_out, se_exp, se_out, parameters):
'''
  Description: This function performs Mendelian Randomization (MR) using the penalized weighted median method.
  
  Parameters:
  - b_exp (array-like): Beta values of the exposure variable.
  - b_out (array-like): Beta values of the outcome variable.
  - se_exp (array-like): Standard errors of the exposure variable.
  - se_out (array-like): Standard errors of the outcome variable.
  - parameters (dict): Dictionary containing additional parameters for the penalized weighted median method, including "penk" and "nboot".
  
  Returns:
  result (dict): Dictionary containing the MR estimates, standard errors, and p-values for the causal effect.
'''
    if np.sum(~np.isnan(b_exp) & ~np.isnan(b_out) & ~np.isnan(se_exp) & ~np.isnan(se_out)) < 3:
        return {
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nsnp": np.nan
        }
    
    betaIV = b_out / b_exp
    betaIVW = np.sum(b_out * b_exp * se_out**-2) / np.sum(b_exp**2 * se_out**-2)
    VBj = (se_out**2) / (b_exp**2) + (b_out**2) * ((se_exp**2) / (b_exp**4))
    weights = 1 / VBj
    bwm = mr_weighted_median(b_exp, b_out, se_exp, se_out, parameters)
    penalty = norm.sf(weights * (betaIV - bwm["b"])**2)
    pen_weights = weights * np.minimum(1, penalty * parameters["penk"])
    b = weighted_median(betaIV, pen_weights)
    se = weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, pen_weights, parameters["nboot"])
    pval = 2 * norm.sf(abs(b / se))
    
    return {
        "b": b,
        "se": se,
        "pval": pval,
        "nsnp": len(b_exp)
    }
