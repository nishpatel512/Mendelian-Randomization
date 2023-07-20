import pandas as pd
import numpy as np
from scipy.stats import t,norm,chi2,stats
from scipy.optimize import minimize
import scipy.stats
import statsmodels.api as sm
from harmonise_data import *
from extract_instruments import *
from extract_outcome_data import *

def mr_raps(b_exp, b_out, se_exp, se_out, parameters):
'''
  Description: This function performs Mendelian Randomization (MR) using the robust adjusted profile score (MR-RAPS) method.
  
  Parameters:
  - b_exp (array-like): Beta values of the exposure variable.
  - b_out (array-like): Beta values of the outcome variable.
  - se_exp (array-like): Standard errors of the exposure variable.
  - se_out (array-like): Standard errors of the outcome variable.
  - parameters (dict): Dictionary containing additional parameters for the MR-RAPS method, including "over_dispersion", "loss_function",
   and "shrinkage".
  
  Returns:
  result (dict): Dictionary containing the MR estimates, robust standard errors, and p-values for the causal effect.
  '''
    data = pd.DataFrame({
        "beta.exposure": b_exp,
        "beta.outcome": b_out,
        "se.exposure": se_exp,
        "se.outcome": se_out
    })
    
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            out = mr_raps(data, diagnostics=False, over_dispersion=parameters["over_dispersion"], loss_function=parameters["loss_function"], shrinkage=parameters["shrinkage"])
    except Exception as e:
        return {
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nsnp": len(b_exp)
        }
    
    return {
        "b": out["beta.hat"],
        "se": out["beta.se"],
        "pval": stats.norm.sf(- np.abs(out["beta.hat"] / out["beta.se"])) * 2,
        "nsnp": len(b_exp)
    }
