import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from blank_plot import blank_plot
from mr_singlesnp import mr_singlesnp
from extract_instruments import extract_instruments
from extract_outcome_data import extract_outcome_data
from harmonise_data import harmonise_data
def mr_forest_plot(singlesnp_results, exponentiate=False):
    '''
  Description: This function creates forest plots for Mendelian Randomization (MR) analysis results based on the output from the mr_singlesnp function.
  
  Parameters:
  - singlesnp_results (pandas.DataFrame): DataFrame containing MR analysis results from the mr_singlesnp function.
  - exponentiate (bool, optional): Whether to exponentiate the effect sizes. Defaults to False.
  
  Returns:
  res (list): A list of matplotlib.pyplot objects containing the generated forest plots.
    '''
    res = []

    for group, d in singlesnp_results.groupby(["id.exposure", "id.outcome"]):
        if sum(~d["SNP"].str.contains("All")) < 2:
            res.append(blank_plot("Insufficient number of SNPs"))
            continue

        d.loc[d["SNP"].str.contains("All - Inverse variance weighted"), "SNP"] = "All - IVW"
        d.loc[d["SNP"].str.contains("All - MR Egger"), "SNP"] = "All - Egger"

        am = d.loc[d["SNP"].str.contains("All"), "SNP"].values
        d["up"] = pd.to_numeric(d["b"]) + 1.96 * pd.to_numeric(d["se"])
        d["lo"] = pd.to_numeric(d["b"]) - 1.96 * pd.to_numeric(d["se"])
        d["tot"] = 0.01
        d.loc[d["SNP"].isin(am), "tot"] = 1

        d["SNP"] = d["SNP"].astype(str)  # Convert SNP column to string
        nom = d.loc[~d["SNP"].isin(am), "SNP"].copy()
        if not nom.empty:
            nom = nom.iloc[np.argsort(pd.to_numeric(d.loc[~d["SNP"].isin(am), "b"]).values)].reset_index(drop=True)
        d = d.append(d.iloc[-1])
        d.iloc[-2, d.columns.get_loc("SNP")] = ""
        d.iloc[-2, d.columns.get_loc("b")] = np.nan
        d.iloc[-2, d.columns.get_loc("up")] = np.nan
        d.iloc[-2, d.columns.get_loc("lo")] = np.nan
        d["SNP"] = pd.Categorical(d["SNP"], categories=np.concatenate((am, [""] + nom)), ordered=True)

        xint = 0
        if exponentiate:
            d["b"] = np.exp(pd.to_numeric(d["b"]))
            d["up"] = np.exp(pd.to_numeric(d["up"]))
            d["lo"] = np.exp(pd.to_numeric(d["lo"]))
            xint = 1

        plt.figure()
        plt.axvline(x=xint, linestyle="dotted")
        #plt.errorbar(x=d["b"].astype(float), y=d["SNP"], xerr=[d["lo"].astype(float), d["up"].astype(float)], fmt="none", color="black", capsize=2)
        #plt.scatter(x=d["b"].astype(float), y=d["SNP"], color="black")
        plt.hlines(y=np.where(d["SNP"] == "", True, False), xmin=0, xmax=1, color="grey")
        plt.ylabel("")
        plt.xlabel(f"MR effect size for\n'{d['exposure'].iloc[0]}' on '{d['outcome'].iloc[0]}'")

        res.append(plt.gca())

    return res


exp = extract_instruments("ieu-a-2")
out = extract_outcome_data(snps=exp["SNP"], outcomes=["ieu-a-7"])

dat = harmonise_data(exp, out)

res = mr_singlesnp(dat)

plt = mr_forest_plot(res)
print(plt)
