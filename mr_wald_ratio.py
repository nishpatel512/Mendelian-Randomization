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
  '''
  Description: This function calculates Mendelian Randomization (MR) estimates using the Wald ratio method.
  
  Parameters:
  - b_exp (array-like or scalar): Beta values of the exposure variable.
  - b_out (array-like or scalar): Beta values of the outcome variable.
  - se_exp (array-like or scalar): Standard errors of the exposure variable.
  - se_out (array-like or scalar): Standard errors of the outcome variable.
  - parameters (dict): Dictionary containing additional parameters for MR analysis.

  Returns:
  result (dict): Dictionary containing the MR estimates, standard errors, p-values, and the number of SNPs used for the analysis.
  '''
  if not np.isscalar(b_exp) and np.size(b_exp) > 1:
      return {'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nsnp': np.nan}

  b = b_out / b_exp
  se = se_out / np.abs(b_exp)
  pval = norm.sf(np.abs(b) / se) * 2

  return {'b': b, 'se': se, 'pval': pval, 'nsnp': 1}


