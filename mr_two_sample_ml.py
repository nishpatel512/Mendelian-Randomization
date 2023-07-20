import pandas as pd
import numpy as np
from scipy.stats import t,norm,chi2,stats
from scipy.optimize import minimize
import scipy.stats
import statsmodels.api as sm
from harmonise_data import *
from extract_instruments import *
from extract_outcome_data import *

def mr_two_sample_ml(b_exp, b_out, se_exp, se_out, parameters):
'''
  Description: This function calculates Mendelian Randomization (MR) estimates using the two-sample instrumental
  variable method with maximum likelihood estimation.   
  
  Parameters:
  - b_exp (array-like): Beta values of the exposure variable.
  - b_out (array-like): Beta values of the outcome variable.
  - se_exp (array-like): Standard errors of the exposure variable.
  - se_out (array-like): Standard errors of the outcome variable.
  - parameters (dict): Dictionary containing additional parameters for MR analysis.
  
  Returns:
  result (dict): Dictionary containing the MR estimates, standard errors, p-values, and other statistics.
'''
    valid_indices = np.where(~np.isnan(b_exp) & ~np.isnan(b_out) & ~np.isnan(se_exp) & ~np.isnan(se_out))
    if len(valid_indices[0]) < 2:
        return {"b": np.nan, "se": np.nan, "pval": np.nan, "nsnp": np.nan, "Q": np.nan, "Q_df": np.nan, "Q_pval": np.nan}

    def loglikelihood(param):
        b_exp_valid = b_exp[valid_indices]
        b_out_valid = b_out[valid_indices]
        se_exp_valid = se_exp[valid_indices]
        se_out_valid = se_out[valid_indices]
        return 0.5 * np.sum(((b_exp_valid - param[:len(b_exp_valid)]) / se_exp_valid)**2) + 0.5 * np.sum(((b_out_valid - param[len(b_exp_valid)] * param[:len(b_exp_valid)]) / se_out_valid)**2)

    try:
        opt = minimize(loglikelihood, np.concatenate((b_exp[valid_indices], [np.sum(b_exp[valid_indices] * b_out[valid_indices] / se_out[valid_indices]**2) / np.sum(b_exp[valid_indices]**2 / se_out[valid_indices]**2)])), hess=True, options={"maxiter": 25000})
        if not opt.success:
            print("mr_two_sample_ml failed to converge")
            return {"b": np.nan, "se": np.nan, "pval": np.nan, "nsnp": np.nan, "Q": np.nan, "Q_df": np.nan, "Q_pval": np.nan}

        b = opt.x[len(b_exp[valid_indices])]
        se = np.sqrt(np.linalg.inv(opt.hess_inv)[len(b_exp[valid_indices]), len(b_exp[valid_indices])])
        pval = 2 * (1 - stats.norm.cdf(np.abs(b) / se))

        Q = 2 * opt.fun
        Q_df = len(b_exp[valid_indices]) - 1
        Q_pval = 1 - chi2.cdf(Q, Q_df)

        return {"b": b, "se": se, "pval": pval, "nsnp": len(valid_indices[0]), "Q": Q, "Q_df": Q_df, "Q_pval": Q_pval}

    except:
        print("mr_two_sample_ml failed to converge")
        return {"b": np.nan, "se": np.nan, "pval": np.nan, "nsnp": np.nan, "Q": np.nan, "Q_df": np.nan, "Q_pval": np.nan}
