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