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
    res = np.zeros(nboot)
    
    for i in range(nboot):
        idx = np.random.choice(len(b_exp), size=len(b_exp), replace=True)
        b_iv = b_out[idx] / b_exp[idx]
        res[i] = weighted_median(b_iv, weights[idx])
    
    return np.std(res, ddof=0)