import pandas as pd
import numpy as np
from scipy.stats import t,norm,chi2,stats
from scipy.optimize import minimize
import statsmodels.api as sm
from harmonise_data import *
from extract_instruments import *
from extract_outcome_data import *
import warnings
def default_parameters():
    """
    Returns a dictionary containing default values for various parameters used in statistical analyses,
    particularly in the context of Mendelian Randomization (MR) methods.

    Returns:
    - A dictionary containing the default parameter values for statistical analyses.
    """
    return {
        "test_dist": "z",
        "nboot": 1000,
        "Cov": 0,
        "penk": 20,
        "phi": 1,
        "alpha": 0.05,
        "Qthresh": 0.05,
        "over.dispersion": True,
        "loss.function": "huber",
        "shrinkage": False
    }


def mr_method_list():
    """
    Generates a DataFrame containing information about various methods used in Mendelian Randomization (MR) analysis.

    Returns:
    - A pandas DataFrame containing the following information for each MR method:
        - 'obj': (string) The identifier or object name representing the MR method.
        - 'name': (string) The human-readable name of the MR method.
        - 'PubmedID': (string) The PubMed ID associated with the reference publication for the MR method, if available.
        - 'Description': (string) A brief description of the MR method and its characteristics.
        - 'use_by_default': (bool) A flag indicating whether the MR method is used by default in the analysis. If True,
                            the method is used by default; otherwise, it's not included by default.
        - 'heterogeneity_test': (bool) A flag indicating whether the MR method includes a heterogeneity test. If True,
                                a heterogeneity test is included; otherwise, it's not part of the method.
    """
    a = [
        {
            'obj': 'mr_wald_ratio',
            'name': 'Wald ratio',
            'PubmedID': '',
            'Description': '',
            'use_by_default': True,
            'heterogeneity_test': False
        },
        {
            'obj': 'mr_two_sample_ml',
            'name': 'Maximum likelihood',
            'PubmedID': '',
            'Description': '',
            'use_by_default': False,
            'heterogeneity_test': True
        },
        {
            'obj': 'mr_egger_regression',
            'name': 'MR Egger',
            'PubmedID': '26050253',
            'Description': '',
            'use_by_default': True,
            'heterogeneity_test': True
        },
        {
            'obj': 'mr_egger_regression_bootstrap',
            'name': 'MR Egger (bootstrap)',
            'PubmedID': '26050253',
            'Description': '',
            'use_by_default': False,
            'heterogeneity_test': False
        },
        {
            'obj': 'mr_simple_median',
            'name': 'Simple median',
            'PubmedID': '',
            'Description': '',
            'use_by_default': False,
            'heterogeneity_test': False
        },
        {
            'obj': 'mr_weighted_median',
            'name': 'Weighted median',
            'PubmedID': '',
            'Description': '',
            'use_by_default': True,
            'heterogeneity_test': False
        },
        {
            'obj': 'mr_penalised_weighted_median',
            'name': 'Penalised weighted median',
            'PubmedID': '',
            'Description': '',
            'use_by_default': False,
            'heterogeneity_test': False
        },
        {
            'obj': 'mr_ivw',
            'name': 'Inverse variance weighted',
            'PubmedID': '',
            'Description': '',
            'use_by_default': True,
            'heterogeneity_test': True
        },
        {
            'obj': 'mr_ivw_mre',
            'name': 'Inverse variance weighted (multiplicative random effects)',
            'PubmedID': '',
            'Description': '',
            'use_by_default': False,
            'heterogeneity_test': False
        },
        {
            'obj': 'mr_ivw_fe',
            'name': 'Inverse variance weighted (fixed effects)',
            'PubmedID': '',
            'Description': '',
            'use_by_default': False,
            'heterogeneity_test': False
        },
        {
            'obj': 'mr_raps',
            'name': 'Robust adjusted profile score (RAPS)',
            'PubmedID': '',
            'Description': '',
            'use_by_default': False,
            'heterogeneity_test': False
        },
        {
            'obj': 'mr_sign',
            'name': 'Sign concordance test',
            'PubmedID': '',
            'Description': 'Tests for concordance of signs between exposure and outcome',
            'use_by_default': False,
            'heterogeneity_test': False
        },
        {
            'obj': 'mr_uwr',
            'name': 'Unweighted regression',
            'PubmedID': '',
            'Description': "Doesn't use any weights",
            'use_by_default': False,
            'heterogeneity_test': True
        }
    ]
    
    df = pd.DataFrame(a)
    df['heterogeneity_test'] = df['heterogeneity_test'].astype(bool)
    df['use_by_default'] = df['use_by_default'].astype(bool)
    return df

def mr(dat, parameters=default_parameters(), method_list=mr_method_list()):
    """
    Perform Mendelian Randomization (MR) analysis on the provided dataset using specified methods.

    Args:
    - dat (pandas DataFrame): The input dataset containing harmonised data of exposure and outcome.
    - parameters (dict, optional): A dictionary containing the default parameters for MR analysis.
                                    Default parameters are used if not provided.
    - method_list (pandas DataFrame, optional): A DataFrame containing information about MR methods.
                                                Default method_list is used if not provided.

    Returns:
    - mr_tab (pandas DataFrame): A DataFrame containing the results of MR analysis for each exposure-outcome pair.
                                 The DataFrame includes columns such as 'outcome', 'exposure', 'method', 'nsnp', 'b', 'se', and 'pval'.
                                 'nsnp' represents the number of SNPs used for analysis, and 'b', 'se', 'pval' are the estimated effect size,
                                 standard error, and p-value, respectively, for each method.

    """
    global mr_keep
    mr_tab = pd.DataFrame()
    for (id_exposure, id_outcome), x in dat.groupby(["id.exposure", "id.outcome"]):
        #x = x1.loc['mr_keep']
        if x.shape[0] == 0:
            print(f"No SNPs available for MR analysis of '{id_exposure}' on '{id_outcome}'")
            continue
        else:
            print(f"Analysing '{id_exposure}' on '{id_outcome}'")
        
        method_list = mr_method_list()['obj'].tolist()
        print (method_list)
        res = {}
        for meth in method_list:
            res[meth] = globals()[meth](x["beta.exposure"], x["beta.outcome"], x["se.exposure"], x["se.outcome"], parameters)
            #res = list(map(lambda meth: globals()[meth](x["beta.exposure"], x["beta.outcome"], x["se.exposure"], x["se.outcome"], parameters), method_list))

        methl = mr_method_list()
        mr_tab = mr_tab.append(pd.DataFrame({
            "outcome": x["outcome"].iloc[0],
            "exposure": x["exposure"].iloc[0],
            "method": methl.loc[methl['obj'].isin(method_list), "name"].tolist(),
            "nsnp": [res[m]["nsnp"] for m in method_list],
            "b": [res[m]["b"] for m in method_list],
            "se": [res[m]["se"] for m in method_list],
            "pval": [res[m]["pval"] for m in method_list]
        }))

    mr_tab = mr_tab.dropna(subset=["b", "se", "pval"])
    return mr_tab

def mr_wald_ratio(b_exp, b_out, se_exp, se_out, parameters):
    """
    Perform Mendelian Randomization (MR) analysis using the Wald ratio method.

    Args:
    - b_exp (numpy array): Array of effect sizes (beta) for the exposure variable.
    - b_out (numpy array): Array of effect sizes (beta) for the outcome variable.
    - se_exp (numpy array): Array of standard errors for the exposure variable effect sizes.
    - se_out (numpy array): Array of standard errors for the outcome variable effect sizes.
    - parameters (dict): A dictionary containing additional parameters for the MR analysis.

    Returns:
    - Dictionary containing the MR analysis results with keys: 'b', 'se', 'pval', and 'nsnp'.
      'b': Estimated effect size (beta) from the Wald ratio method.
      'se': Standard error of the estimated effect size.
      'pval': P-value for the effect size estimate.
      'nsnp': Number of SNPs used in the analysis.
    """
    if len(b_exp) > 1:
        return {'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nsnp': np.nan}
    
    b = b_out / b_exp
    se = se_out / np.abs(b_exp)
    pval = norm.sf(np.abs(b) / se) * 2
    
    return {'b': b, 'se': se, 'pval': pval, 'nsnp': 1}

def mr_two_sample_ml(b_exp, b_out, se_exp, se_out, parameters):
    """
    Perform Mendelian Randomization (MR) analysis using the two-sample Maximum Likelihood (ML) method.

    Args:
    - b_exp (numpy array): Array of effect sizes (beta) for the exposure variable.
    - b_out (numpy array): Array of effect sizes (beta) for the outcome variable.
    - se_exp (numpy array): Array of standard errors for the exposure variable effect sizes.
    - se_out (numpy array): Array of standard errors for the outcome variable effect sizes.
    - parameters (dict): A dictionary containing additional parameters for the MR analysis.

    Returns:
    - Dictionary containing the MR analysis results with keys: 'b', 'se', 'pval', 'nsnp', 'Q', 'Q_df', and 'Q_pval'.
      'b': Estimated causal effect size from the two-sample ML method.
      'se': Standard error of the estimated causal effect size.
      'pval': P-value for the causal effect estimate.
      'nsnp': Number of SNPs used in the analysis.
      'Q': The test statistic for heterogeneity (Q statistic) in the MR analysis.
      'Q_df': Degrees of freedom associated with the Q statistic.
      'Q_pval': P-value for the Q statistic.
    """
    valid_indices = np.where(~np.isnan(b_exp) & ~np.isnan(b_out) & ~np.isnan(se_exp) & ~np.isnan(se_out))
    if len(valid_indices[0]) < 2:
        return {"b": np.nan, "se": np.nan, "pval": np.nan, "nsnp": np.nan, "Q": np.nan, "Q_df": np.nan, "Q_pval": np.nan}

    def loglikelihood(param):
        """
        Calculate the log-likelihood function for the two-sample Maximum Likelihood (ML) method in Mendelian Randomization (MR) analysis.

        Args:
        - param (numpy array): Array containing the parameter values to be used in the log-likelihood function.

        Returns:
        - float: The value of the log-likelihood function.
        """
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


def mr_egger_regression(b_exp, b_out, se_exp, se_out, parameters):
    '''
    Input: Beta values of exposure and outcomes\n
    Output: 
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
    weights = 1 / se_out ** 2

    model = sm.WLS(b_out, X, weights=weights)
    results = model.fit()

    b = results.params[1]
    se = results.bse[1]
    pval = 2 * stats.t.sf(abs(b / se), len(b_exp) - 2)

    b_i = results.params[0]
    se_i = results.bse[0]
    pval_i = 2 * stats.t.sf(abs(b_i / se_i), len(b_exp) - 2)

    Q = results.scale * (len(b_exp) - 2)
    Q_df = len(b_exp) - 2
    Q_pval = stats.chi2.sf(Q, Q_df)

    return {
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

def mr_simple_median(b_exp, b_out, se_exp, se_out, parameters):
    if np.sum(~np.isnan(b_exp) & ~np.isnan(b_out) & ~np.isnan(se_exp) & ~np.isnan(se_out)) < 3:
        return {
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nsnp": np.nan
        }
    
    b_iv = b_out / b_exp
    b = weighted_median(b_iv, np.repeat(1 / len(b_exp), len(b_exp)))
    se = weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, np.repeat(1 / len(b_exp), len(b_exp)), parameters["nboot"])
    pval = 2 * norm.sf(abs(b / se))
    
    return {
        "b": b,
        "se": se,
        "pval": pval,
        "nsnp": len(b_exp)
    }

def weighted_median(b_iv, weights):
    betaIV_order = b_iv[np.argsort(b_iv)]
    weights_order = weights[np.argsort(b_iv)]
    weights_sum = np.cumsum(weights_order) - 0.5 * weights_order
    weights_sum = weights_sum / np.sum(weights_order)
    below = np.max(np.where(weights_sum < 0.5))
    b = betaIV_order[below] + (betaIV_order[below+1] - betaIV_order[below]) * (0.5 - weights_sum[below]) / (weights_sum[below+1] - weights_sum[below])
    return b

def weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, weights, nboot):
    res = np.zeros(nboot)
    
    for i in range(nboot):
        idx = np.random.choice(len(b_exp), size=len(b_exp), replace=True)
        b_iv = b_out[idx] / b_exp[idx]
        res[i] = weighted_median(b_iv, weights[idx])
    
    return np.std(res, ddof=0)

def weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, weights, nboot):
    med = np.zeros(nboot)
    
    for i in range(nboot):
        b_exp_boot = np.random.normal(loc=b_exp, scale=se_exp, size=len(b_exp))
        b_out_boot = np.random.normal(loc=b_out, scale=se_out, size=len(b_out))
        betaIV_boot = b_out_boot / b_exp_boot
        med[i] = weighted_median(betaIV_boot, weights)
    
    return np.std(med)

def mr_penalised_weighted_median(b_exp, b_out, se_exp, se_out, parameters):
    if np.sum(~np.isnan(b_exp) & ~np.isnan(b_out) & ~np.isnan(se_exp) & ~np.isnan(se_out)) < 3:
        return {
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nsnp": np.nan
        }
    
    betaIV = b_out / b_exp
    betaIVW = np.sum(b_out * b_exp * se_out**-2) / np.sum(b_exp**2 * se_out**-2)
    VBj = (se_out**2) / (b_exp**2) + (b_out**2) * ((se_exp**2) / (b_exp**4))
    weights = 1 / VBj
    bwm = mr_weighted_median(b_exp, b_out, se_exp, se_out, parameters)
    penalty = norm.sf(weights * (betaIV - bwm["b"])**2)
    pen_weights = weights * np.minimum(1, penalty * parameters["penk"])
    b = weighted_median(betaIV, pen_weights)
    se = weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, pen_weights, parameters["nboot"])
    pval = 2 * norm.sf(abs(b / se))
    
    return {
        "b": b,
        "se": se,
        "pval": pval,
        "nsnp": len(b_exp)
    }

def mr_median(dat, parameters):
    if "mr_keep" in dat.columns:
        dat = dat[dat["mr_keep"]]
    
    if len(dat) < 3:
        print("Need at least 3 SNPs")
        return None
    
    b_exp = dat["beta.exposure"]
    b_out = dat["beta.outcome"]
    se_exp = dat["se.exposure"]
    se_out = dat["se.outcome"]
    
    sm = mr_simple_median(b_exp, b_out, se_exp, se_out, parameters)
    wm = mr_weighted_median(b_exp, b_out, se_exp, se_out, parameters)
    pm = mr_penalised_weighted_median(b_exp, b_out, se_exp, se_out, parameters)
    
    res = {
        "id.exposure": dat["id.exposure"].iloc[0],
        "id.outcome": dat["id.outcome"].iloc[0],
        "method": ["Simple median", "Weighted median", "Penalised median"],
        "nsnp": len(b_exp),
        "b": [sm["b"], wm["b"], pm["b"]],
        "se": [sm["se"], wm["se"], pm["se"]],
        "ci_low": [sm["b"] - norm.ppf(1-parameters["alpha"]/2) * sm["se"],
                   wm["b"] - norm.ppf(1-parameters["alpha"]/2) * wm["se"],
                   pm["b"] - norm.ppf(1-parameters["alpha"]/2) * pm["se"]],
        "ci_upp": [sm["b"] + norm.ppf(1-parameters["alpha"]/2) * sm["se"],
                   wm["b"] + norm.ppf(1-parameters["alpha"]/2) * wm["se"],
                   pm["b"] + norm.ppf(1-parameters["alpha"]/2) * pm["se"]],
        "pval": [sm["pval"], wm["pval"], pm["pval"]]
    }
    
    return pd.DataFrame(res)

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

def mr_uwr(b_exp, b_out, se_exp, se_out, parameters):
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
    
    ivw_res = sm.OLS(b_out, sm.add_constant(b_exp)).fit()
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

def mr_ivw_mre(b_exp, b_out, se_exp, se_out, parameters):
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
    se = ivw_res.bse[1]
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

def mr_ivw_fe(b_exp, b_out, se_exp, se_out, parameters):
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
    se = ivw_res.bse[1] / ivw_res.mse_resid**0.5
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
# e = extract_instruments("ieu-a-2")
# o = extract_outcome_data(snps=e["SNP"], outcomes=["ieu-a-7"])

# dat = harmonise_data(e , o)

# res = mr(dat)
# print (res)

