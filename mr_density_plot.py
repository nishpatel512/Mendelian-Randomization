import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from blank_plot import blank_plot
from mr_singlesnp import mr_singlesnp
from extract_instruments import extract_instruments
from extract_outcome_data import extract_outcome_data
from harmonise_data import harmonise_data
from mr import *
def get_method_obj(method):
    """
    Retrieve the method object from the list of MR methods.

    Parameters:
    - method (str): Name of the MR method.

    Returns:
    - method_obj (object): The method object associated with the provided method name.
    """
    method_index = mr_method_list().index(method)
    method_obj = mr_method_list()[method_index]
    return method_obj

def mr_density_plot(singlesnp_results, mr_results, exponentiate=False, bandwidth="nrd0"):
    """
    Create density plots for MR estimates using a histogram and per SNP MR estimates.

    Parameters:
    - singlesnp_results (pandas DataFrame): Dataframe containing the per-SNP MR results.
    - mr_results (pandas DataFrame): Dataframe containing the MR results for different methods.
    - exponentiate (bool, optional): Whether to exponentiate the MR estimates or not. Default is False.
    - bandwidth (str, optional): Bandwidth method for kernel density estimation. Default is "nrd0".

    Returns:
    - res (list): A list of matplotlib.pyplot objects containing the generated density plots.
    """
    res = []
    for (id_exposure, id_outcome), d in singlesnp_results.groupby(["id.exposure", "id.outcome"]):
        if np.sum(~d["SNP"].str.contains("All")) < 2:
            return blank_plot("Insufficient number of SNPs")
        d["SNP"] = d["SNP"].astype(str)

        d2 = d[~d["SNP"].str.contains("All - ")]
        # d1 = mr_results[(mr_results["id.exposure"] == d2["id.exposure"].iloc[0]) & (mr_results["id.outcome"] == d2["id.outcome"].iloc[0])]

        xint = 0
        if exponentiate:
            d["b"] = np.exp(d["b"])
            d["up"] = np.exp(d["up"])
            d["lo"] = np.exp(d["lo"])
            xint = 1

        plt.figure()
        plt.hist(d2["b"], weights=1/d2["se"], density=True, bins=20, label="Histogram")
        plt.axvline(x=xint, linestyle="dotted", label="Exponentiated Reference Line")
        plt.scatter(x=d2["b"], y=np.zeros(len(d2)), color="red", s=1/d2["se"], label="Per SNP MR")

        for method, b in zip(mr_results["method"], mr_results["b"]):
            if method == "MR Egger":
                plt.axvline(x=b, label="MR Egger Estimate")
                break  # Assuming MR Egger appears only once in the mr_results

        plt.xlabel("Per SNP MR estimate")
        plt.ylabel("Density")
        plt.legend()
        plt.show()
        res.append(plt)
    return res
# exp = extract_instruments("ieu-a-2")
# out = extract_outcome_data(snps=exp["SNP"], outcomes=["ieu-a-7"])

# data = harmonise_data(exp, out)
# res = mr_singlesnp(data)
# result = mr(data)
# mr_density_plot(res, result)

