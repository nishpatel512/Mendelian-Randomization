import pandas as pd
import numpy as np
from scipy.stats import t,norm,chi2,stats
from scipy.optimize import minimize
import scipy.stats
import statsmodels.api as sm
from harmonise_data import *
from extract_instruments import *
from extract_outcome_data import *

def weighted_median(b_iv, weights):
'''
    Description: The weighted_median function calculates the weighted median of a given set of beta values (b_iv) using provided weights.
    It is used in Mendelian Randomization (MR) analysis to estimate causal effects by combining instrumental variable (IV) effect estimates.
    
    Parameters:
    -b_iv (array-like): Array of instrumental variable effect estimates (betaIV) for each SNP.
    -weights (array-like): Array of weights corresponding to each betaIV estimate.
    
    Returns:
    b (float): Weighted median estimate of the instrumental variable effect (b_iv) using the given weights.
'''
    betaIV_order = b_iv[np.argsort(b_iv)]
    weights_order = weights[np.argsort(b_iv)]
    weights_sum = np.cumsum(weights_order) - 0.5 * weights_order
    weights_sum = weights_sum / np.sum(weights_order)
    below = np.max(np.where(weights_sum < 0.5))
    b = betaIV_order[below] + (betaIV_order[below+1] - betaIV_order[below]) * (0.5 - weights_sum[below]) / (weights_sum[below+1] - weights_sum[below])
    return b
