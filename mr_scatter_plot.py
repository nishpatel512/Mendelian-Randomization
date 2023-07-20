import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mr_egger_regression import *
import mr_egger_regression_bootstrap
from mr import *
from extract_instruments import *
import blank_plot
from extract_outcome_data import *
import warnings



def mr_scatter_plot(mr_results, dat):
'''
  Description: This function creates scatter plots to compare SNP effects on the exposure variable with SNP effects on the
  outcome variable for Mendelian Randomization (MR) analysis results.
  
  Parameters:
  - mr_results (pandas.DataFrame): DataFrame containing MR analysis results.
  - dat (pandas.DataFrame): DataFrame containing harmonized data.
  
  Returns:
  res (list): A list of matplotlib.pyplot objects containing the generated scatter plots.
'''
    res = []

    for group, d in dat.groupby(["id.exposure", "id.outcome"]):
        if len(d) < 2 or sum(d["mr_keep"]) == 0:
            res.append(blank_plot("Insufficient number of SNPs"))
            continue

        d = d[d["mr_keep"]]

        index = d["beta.exposure"] < 0
        d.loc[index, "beta.exposure"] *= -1
        d.loc[index, "beta.outcome"] *= -1
        mrres = mr_results[(mr_results["exposure"] == d["id.exposure"].iloc[0]) & (mr_results["outcome"] == d["id.outcome"].iloc[0])]
        mrres["a"] = 0

        if "MR Egger" in mrres["method"].values:
            temp = mr_egger_regression(d["beta.exposure"], d["beta.outcome"], d["se.exposure"], d["se.outcome"], default_parameters())
            mrres.loc[mrres["method"] == "MR Egger", "a"] = temp.b_i

        if "MR Egger (bootstrap)" in mrres["method"].values:
            temp = mr_egger_regression_bootstrap(d["beta.exposure"], d["beta.outcome"], d["se.exposure"], d["se.outcome"], default_parameters())
            mrres.loc[mrres["method"] == "MR Egger (bootstrap)", "a"] = temp.b_i

        plt.figure()
        plt.scatter(d["beta.exposure"], d["beta.outcome"])
        plt.errorbar(x=d["beta.exposure"], y=d["beta.outcome"], xerr=d["se.exposure"], yerr=d["se.outcome"], fmt="none", color="grey")
        plt.errorbar(x=d["beta.exposure"], y=d["beta.outcome"], xerr=d["se.exposure"], yerr=d["se.outcome"], fmt="none", color="grey", alpha=0.5, capsize=0)
        for _, row in mrres.iterrows():
            plt.plot(d["beta.exposure"], row["a"] + row["b"] * d["beta.exposure"], label=row["method"])
        plt.xlabel(f"SNP effect on {d['exposure'].iloc[0]}")
        plt.ylabel(f"SNP effect on {d['outcome'].iloc[0]}")
        plt.legend(title="MR Test", loc="upper right")

        res.append(plt.gca())
    return res

e = extract_instruments("ieu-a-2")
o = extract_outcome_data(snps=e["SNP"], outcomes=["ieu-a-7"])

dat = harmonise_data(e , o)

res = mr(dat)
print (mr_scatter_plot(res,dat))
