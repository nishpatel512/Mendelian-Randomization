import pandas as pd
import numpy as np
from scipy.stats import t,norm,chi2,stats
from scipy.optimize import minimize
import scipy.stats
import statsmodels.api as sm
from harmonise_data import *
from extract_instruments import *
from extract_outcome_data import *

def linreg(x, y, w=None):
    if w is None:
        w = np.ones_like(x)
    
    xp = w * x
    yp = w * y
    bhat = np.cov(xp, yp, ddof=0)[0, 1] / np.var(xp, ddof=0)
    ahat = np.mean(yp) - np.mean(xp) * bhat
    yhat = ahat + bhat * x
    se = np.sqrt(np.sum((yp - yhat)**2) / (np.sum(~np.isnan(yhat)) - 2) / np.sum(w * x**2))
    pval = 2 * t.sf(abs(bhat / se), np.sum(~np.isnan(yhat)) - 2)
    
    return {
        "ahat": ahat,
        "bhat": bhat,
        "se": se,
        "pval": pval
    }

def mr_egger_regression_bootstrap(b_exp, b_out, se_exp, se_out, parameters):
    if np.sum(~np.isnan(b_exp) & ~np.isnan(b_out) & ~np.isnan(se_exp) & ~np.isnan(se_out)) < 3:
        return {
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nsnp": np.nan,
            "b_i": np.nan,
            "se_i": np.nan,
            "pval_i": np.nan,
            "mod": np.nan,
            "smod": np.nan,
            "dat": np.nan
        }
    
    nboot = parameters["nboot"]
    res = np.zeros((nboot + 1, 2))
    
    for i in range(nboot):
        xs = np.random.normal(b_exp, se_exp)
        ys = np.random.normal(b_out, se_out)
        
        ys *= np.sign(xs)
        xs = np.abs(xs)
        
        r = linreg(xs, ys, 1 / se_out**2)
        
        res[i, 0] = r["ahat"]
        res[i, 1] = r["bhat"]
    
    return {
        "b": np.mean(res[:, 1], axis=0),
        "se": np.std(res[:, 1], axis=0),
        "pval": np.sum(np.sign(np.mean(res[:, 1], axis=0)) * res[:, 1] < 0) / nboot,
        "nsnp": len(b_exp),
        "b_i": np.mean(res[:, 0], axis=0),
        "se_i": np.std(res[:, 0], axis=0),
        "pval_i": np.sum(np.sign(np.mean(res[:, 0], axis=0)) * res[:, 0] < 0) / nboot
    }
