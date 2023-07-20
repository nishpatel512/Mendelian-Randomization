import pandas as pd
import numpy as np
from scipy.stats import t,norm,chi2,stats
from scipy.optimize import minimize
import scipy.stats
import statsmodels.api as sm
from harmonise_data import *
from extract_instruments import *
from extract_outcome_data import *
from weighted_median import *
from weighted_median_bootstrap import *

def mr_weighted_median(b_exp, b_out, se_exp, se_out, parameters):
'''
  Description: This function calculates Mendelian Randomization (MR) estimates using the weighted median method.
  
  Parameters:
  - b_exp (array-like): Beta values of the exposure variable.
  - b_out (array-like): Beta values of the outcome variable.
  - se_exp (array-like): Standard errors of the exposure variable.
  - se_out (array-like): Standard errors of the outcome variable.
  - parameters (dict): Dictionary containing additional parameters for MR analysis.

  Returns:
  result (dict): Dictionary containing the MR estimates, standard errors, p-values, Q-statistic, degrees of freedom (Q_df),
                 and p-value of the Q-statistic (Q_pval), as well as the number of SNPs used for the analysis (nsnp).
'''
    if np.sum(~np.isnan(b_exp) & ~np.isnan(b_out) & ~np.isnan(se_exp) & ~np.isnan(se_out)) < 3:
        return {
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nsnp": np.nan
        }
    
    b_iv = b_out / b_exp
    VBj = (se_out**2) / (b_exp**2) + (b_out**2) * ((se_exp**2)) / (b_exp**4)
    b = weighted_median(b_iv, 1 / VBj)
    se = weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, 1 / VBj, parameters["nboot"])
    pval = 2 * norm.sf(abs(b / se))
    
    return {
        "b": b,
        "se": se,
        "pval": pval,
        "Q": np.nan,
        "Q_df": np.nan,
        "Q_pval": np.nan,
        "nsnp": len(b_exp)
    }
