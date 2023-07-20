import pandas as pd
import numpy as np
from format_data import *
import pandas as pd
import numpy as np
from format_data import *

def convert_outcome_to_exposure(outcome_dat):
    """
    Convert outcome data to exposure data format.

    Parameters:
        outcome_dat (pandas.DataFrame): The outcome data.

    Returns:
        pandas.DataFrame: The converted exposure data.
    """
    id = outcome_dat[~outcome_dat.duplicated('outcome')][['outcome', 'id.outcome']]
    exposure_dat = format_data(outcome_dat,
                               beta_col="beta.outcome",
                               se_col="se.outcome",
                               pval_col="pval.outcome",
                               phenotype_col="outcome",
                               effect_allele_col="effect_allele.outcome",
                               other_allele_col="other_allele.outcome",
                               eaf_col="eaf.outcome",
                               units_col="units.outcome")
    exposure_dat = exposure_dat.merge(id, left_on='exposure', right_on='outcome')
    exposure_dat = exposure_dat.drop(columns=['id.exposure'])
    exposure_dat.rename(columns={'id.outcome': 'id.exposure'}, inplace=True)
    return exposure_dat
