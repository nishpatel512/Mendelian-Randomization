import pandas as pd
import numpy as np
from scipy.stats import t,norm,chi2,stats
from scipy.optimize import minimize
import scipy.stats
import statsmodels.api as sm
from harmonise_data import *
from extract_instruments import *
from extract_outcome_data import *

def mr_ivw(b_exp, b_out, se_exp, se_out, parameters):
    if np.sum(~np.isnan(b_exp) & ~np.isnan(b_out) & ~np.isnan(se_exp) & ~np.isnan(se_out)) < 2:
        return {
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nsnp": np.nan,
            "Q": np.nan,
            "Q_df": np.nan,
            "Q_pval": np.nan
        }
    
    ivw_res = sm.OLS(b_out, sm.add_constant(b_exp), weights=1/se_out**2).fit()
    b = ivw_res.params[1]
    se = ivw_res.bse[1] / min(1, ivw_res.mse_resid)
    pval = 2 * norm.sf(abs(b / se))
    Q_df = len(b_exp) - 1
    Q = ivw_res.mse_resid * Q_df
    Q_pval = chi2.sf(Q, Q_df)
    
    return {
        "b": b,
        "se": se,
        "pval": pval,
        "nsnp": len(b_exp),
        "Q": Q,
        "Q_df": Q_df,
        "Q_pval": Q_pval
    }
