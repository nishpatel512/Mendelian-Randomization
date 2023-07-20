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

def weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, weights, nboot):
    '''
    Description: The weighted_median_bootstrap function calculates the standard error of the weighted median estimate of instrumental
    variable (IV) effect in Mendelian Randomization (MR) analysis using bootstrap resampling.
    
    Parameters:  
    -b_exp (array-like): Beta values of the exposure variable.
    -b_out (array-like): Beta values of the outcome variable.
    -se_exp (array-like): Standard errors of the exposure variable.
    -se_out (array-like): Standard errors of the outcome variable.
    -weights (array-like): Array of weights corresponding to each IV effect estimate.
    -nboot (int): Number of bootstrap iterations.
    
    Returns:
    se (float): Standard error of the weighted median estimate calculated using bootstrap resampling.
  '''
    res = np.zeros(nboot)
    
    for i in range(nboot):
        idx = np.random.choice(len(b_exp), size=len(b_exp), replace=True)
        b_iv = b_out[idx] / b_exp[idx]
        res[i] = weighted_median(b_iv, weights[idx])
    
    return np.std(res, ddof=0)
