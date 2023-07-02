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
