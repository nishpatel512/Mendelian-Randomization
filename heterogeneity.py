import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from harmonise_data import *
from mr import *

def mr_heterogeneity(dat, parameters=None, method_list=None):
    """
    Function to perform heterogeneity analysis on a given dataframe.

    Parameters:
    - dat: The input dataframe.
    - parameters: Optional parameter values. If not provided, default parameters will be used.
    - method_list: Optional list of methods. If not provided, a default method list will be used.

    Returns:
    - het_tab: A dataframe containing the results of the heterogeneity analysis.
    """
    if parameters is None:
        # Define default_parameters() function in Python to set default parameters
        parameters = default_parameters()
    if method_list is None:
        # Define mr_method_list() function in Python to get the method list
        method_list = mr_method_list()

    def analyze_group(x1):
        """
        Function to analyze a group within the input dataframe.

        Parameters:
        - x1: The group dataframe.

        Returns:
        - het_tab: A dataframe containing the results of the analysis for the group.
        """
        x = x1[x1["mr_keep"]]
        if len(x) < 2:
            print(f"Not enough SNPs available for Heterogeneity analysis of '{x1['id.exposure'].iloc[0]}' on '{x1['id.outcome'].iloc[0]}'")
            return None

        method_list = mr_method_list()['obj'].tolist()
        res = []
        for meth in method_list:
            result = globals()[meth](x["beta.exposure"], x["beta.outcome"], x["se.exposure"], x["se.outcome"], parameters)
            if isinstance(result, dict) and "Q" in result and "Q_df" in result and "Q_pval" in result:
                result.update({"method": meth})
                res.append(result)

        het_tab = pd.DataFrame(res)
        het_tab = het_tab.dropna(subset=["Q", "Q_df", "Q_pval"])
        return het_tab

    # Create an empty list to store individual results for each group
    result_list = []

    # Analyze each group and add the DataFrame to the result_list
    for _, group_df in dat.groupby(["id.exposure", "id.outcome"]):
        analyzed_group = analyze_group(group_df)
        if analyzed_group is not None:
            result_list.append(analyzed_group)

    # Combine results for all groups into a single DataFrame
    combined_res = pd.concat(result_list, ignore_index=True)
    return combined_res

def mr_pleiotropy_test(dat):
    """
    This function performs a pleiotropy analysis on groups within a DataFrame.
    It takes a DataFrame `dat` as input.

    Parameters:
    - dat: The input DataFrame containing the data for analysis.

    Returns:
    - ptab: The results of the pleiotropy analysis as a DataFrame.
    """
    def analyze_group(x1):
        """
        This function performs the pleiotropy analysis on a single group within the input DataFrame.
        It takes a DataFrame `x1` as input.

        Parameters:
        - x1: The input DataFrame containing the data for a single group.

        Returns:
        - out: The results of the pleiotropy analysis for the group as a DataFrame.
        """
        x = x1[x1["mr_keep"]]
        if len(x) < 2:
            print(f"Not enough SNPs available for pleiotropy analysis of '{x1['id.exposure'].iloc[0]}' on '{x1['id.outcome'].iloc[0]}'")
            return None

        res = mr_egger_regression(x["beta.exposure"], x["beta.outcome"], x["se.exposure"], x["se.outcome"], default_parameters())
        out = pd.DataFrame({
            "outcome": x["outcome"].iloc[0],
            "exposure": x["exposure"].iloc[0],
            "egger_intercept": res["b_i"],
            "se": res["se_i"],
            "pval": res["pval_i"]
        }, index=[0])

        return out

    ptab = dat.groupby(["id.exposure", "id.outcome"]).apply(analyze_group)
    return ptab.reset_index(drop=True)

# Example usage
e = extract_instruments("ieu-a-2")
o = extract_outcome_data(snps=e["SNP"], outcomes=["ieu-a-7"])

dat = harmonise_data(e, o)

res = mr_heterogeneity(dat)

print(res)