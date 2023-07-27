import pandas as pd
import numpy as np
from scipy.stats import t,norm,chi2,stats
from scipy.optimize import minimize
import scipy.stats
import statsmodels.api as sm
from harmonise_data import *
from extract_instruments import *
from extract_outcome_data import *
import statsmodels.api as sm
#from scipy.optimize import WeightedLeastSquares
def mr_egger_regression(b_exp, b_out, se_exp, se_out, parameters):
    '''
  Description: The harmonise_data function is used to perform harmonization of genetic effect sizes and alleles for Mendelian Randomization (MR) analysis.
  It ensures that the effect of a Single Nucleotide Polymorphism (SNP) on both the exposure and outcome variables is measured relative to the same allele.
  Parameters:
  -exposure_data (pandas.DataFrame): DataFrame containing the exposure data with columns such as "SNP," "effect_size," "standard_error," "effect_allele," and "other_allele."
  -outcome_data (pandas.DataFrame): DataFrame containing the outcome data with columns such as "SNP," "effect_size," "standard_error," "effect_allele," and "other_allele."
  Returns:
  harmonized_data (pandas.DataFrame): Data frame with harmonized genetic effect sizes and alleles, where the effect sizes are relative to the same
  allele for both exposure and outcome 
  data.
    '''
    if len(b_exp) != len(b_out) or len(se_exp) != len(se_out) or len(b_exp) != len(se_out):
        return {
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nsnp": np.nan,
            "b_i": np.nan,
            "se_i": np.nan,
            "pval_i": np.nan,
            "Q": np.nan,
            "Q_df": np.nan,
            "Q_pval": np.nan,
            "mod": np.nan,
            "smod": np.nan,
            "dat": np.nan
        }

    nulllist = {
        "b": np.nan,
        "se": np.nan,
        "pval": np.nan,
        "nsnp": np.nan,
        "b_i": np.nan,
        "se_i": np.nan,
        "pval_i": np.nan,
        "Q": np.nan,
        "Q_df": np.nan,
        "Q_pval": np.nan,
        "mod": np.nan,
        "smod": np.nan,
        "dat": np.nan
    }

    if np.sum(~np.isnan(b_exp) & ~np.isnan(b_out) & ~np.isnan(se_exp) & ~np.isnan(se_out)) < 3:
        return nulllist

    def sign0(x):
        x[x == 0] = 1
        return np.sign(x)

    to_flip = sign0(b_exp) == -1
    b_out = b_out * sign0(b_exp)
    b_exp = np.abs(b_exp)
    dat = {
        "b_out": b_out,
        "b_exp": b_exp,
        "se_exp": se_exp,
        "se_out": se_out,
        "flipped": to_flip
    }

    X = sm.add_constant(b_exp)
    global weights
    weights = 1 / se_out ** 2
    mod = sm.WLS(b_out, X, weights=weights)
    results = mod.fit()
    # weights = 1 / se_out**2
    # X = np.column_stack((np.ones_like(b_exp), b_exp))
    # Y = b_out
    # W = np.diag(weights)

    # model = WeightedLeastSquares(X, Y, weights=W)
    # results = model.fit()
    #beta = results.params
    # se = np.sqrt(results.cov[1, 1])


    b = results.params[1]
    se = results.bse[1]
    pval = 2 * stats.t.sf(abs(b / se), df=len(b_exp) - 2)

    b_i = results.params[0]
    se_i = results.bse[0]
    pval_i = 2 * stats.t.sf(abs(b_i / se_i),df=len(b_exp) - 2)

    Q = results.scale * (len(b_exp) - 2)
    Q_df = len(b_exp) - 2
    Q_pval = stats.chi2.sf(Q, Q_df)

    return {
        # "b": results.params[1],
        # "se": results.bse[1],
        # "pval": results.pvalues[1],
        "b": b,
        "se": se,
        "pval": pval,
        "nsnp": len(b_exp),
        "b_i": b_i,
        "se_i": se_i,
        "pval_i": pval_i,
        "Q": Q,
        "Q_df": Q_df,
        "Q_pval": Q_pval,
        "mod": results,
        "dat": dat
    }
