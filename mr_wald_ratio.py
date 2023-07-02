import pandas as pd
import numpy as np
from scipy.stats import t,norm,chi2,stats
from scipy.optimize import minimize
import scipy.stats
import statsmodels.api as sm
from harmonise_data import *
from extract_instruments import *
from extract_outcome_data import *

def mr_wald_ratio(b_exp, b_out, se_exp, se_out, parameters):
    if len(b_exp) > 1:
        return {'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nsnp': np.nan}
    
    b = b_out / b_exp
    se = se_out / np.abs(b_exp)
    pval = norm.sf(np.abs(b) / se) * 2
    
    return {'b': b, 'se': se, 'pval': pval, 'nsnp': 1}

