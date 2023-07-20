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

def mr_simple_median(b_exp, b_out, se_exp, se_out, parameters):
'''
  Description: This function calculates Mendelian Randomization (MR) estimates using the simple median method.
  
  Parameters:
  - b_exp (array-like): Beta values of the exposure variable.
  - b_out (array-like): Beta values of the outcome variable.
  - se_exp (array-like): Standard errors of the exposure variable.
  - se_out (array-like): Standard errors of the outcome variable.
  - parameters (dict): Dictionary containing additional parameters, including "nboot" for the number of bootstrap iterations.
  
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
    
    b_iv = b_out / b_exp
    b = weighted_median(b_iv, np.repeat(1 / len(b_exp), len(b_exp)))
    se = weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, np.repeat(1 / len(b_exp), len(b_exp)), parameters["nboot"])
    pval = 2 * norm.sf(abs(b / se))
    
    return {
        "b": b,
        "se": se,
        "pval": pval,
        "nsnp": len(b_exp)
    }
