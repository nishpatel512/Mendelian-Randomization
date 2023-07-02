import pandas as pd
import numpy as np
from scipy.stats import t,norm,chi2,stats
from scipy.optimize import minimize
import scipy.stats
import statsmodels.api as sm
from harmonise_data import *
from extract_instruments import *
from extract_outcome_data import *

def mr_sign(b_exp, b_out, se_exp, se_out, parameters):
    b_exp = np.copy(b_exp)
    b_out = np.copy(b_out)
    b_exp[b_exp == 0] = np.nan
    b_out[b_out == 0] = np.nan
    betaIV = b_out / b_exp
    
    return {
        "b": np.nanmedian(betaIV),
        "se": np.nanstd(betaIV),
        "pval": 2 * stats.norm.sf(np.abs(np.nanmedian(betaIV) / np.nanstd(betaIV))),
        "nsnp": np.sum(~np.isnan(b_exp) & ~np.isnan(b_out) & ~np.isnan(se_exp) & ~np.isnan(se_out))
    }