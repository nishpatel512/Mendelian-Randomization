import pandas as pd
import numpy as np
from mr import mr_method_list
from mr_ivw import mr_ivw
from mr_wald_ratio import mr_wald_ratio
from extract_instruments import extract_instruments
from extract_outcome_data import extract_outcome_data
from harmonise_data import harmonise_data
from mr_egger_regression import *

def default_parameters():
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

def mr_singlesnp(dat, parameters=None, single_method="mr_wald_ratio", all_method=("mr_ivw", "mr_egger_regression")):
    '''
      Description: This function performs Mendelian Randomization (MR) analysis for each exposure-outcome pair.
      Parameters:
      - b_exp (array-like): Beta values of the exposure variable.
      - b_out (array-like): Beta values of the outcome variable.
      - se_exp (array-like): Standard errors of the exposure variable.
      - se_out (array-like): Standard errors of the outcome variable.
      - parameters (dict): Dictionary containing additional parameters for MR analysis.
      Returns:
      result (DataFrame): A pandas DataFrame containing the MR estimates, standard errors, and p-values for each SNP within each
      exposure-outcome pair.
    '''
    if parameters is None:
        parameters = default_parameters()

    if "samplesize.outcome" not in dat.columns:
        dat["samplesize.outcome"] = np.nan

    assert "outcome" in dat.columns, "Column 'outcome' is missing in the input data"
    assert "exposure" in dat.columns, "Column 'exposure' is missing in the input data"
    assert "beta.exposure" in dat.columns, "Column 'beta.exposure' is missing in the input data"
    assert "beta.outcome" in dat.columns, "Column 'beta.outcome' is missing in the input data"
    assert "se.exposure" in dat.columns, "Column 'se.exposure' is missing in the input data"
    assert "se.outcome" in dat.columns, "Column 'se.outcome' is missing in the input data"

    res = pd.DataFrame(columns=["exposure", "outcome", "id.exposure", "id.outcome", "samplesize", "SNP", "b", "se", "p"])
    
    for group, x in dat.groupby(["id.exposure", "id.outcome"]):
        x = x[x['mr_keep']]
        nsnp = len(x)

        if nsnp == 0:
            x = x.iloc[0]
            d = pd.DataFrame({
                "SNP": "No available data",
                "b": np.nan,
                "se": np.nan,
                "p": np.nan,
                "samplesize": np.nan,
                "outcome": x["outcome"],
                "exposure": x["exposure"]
            }, index=[0])
            res = res.append(d, ignore_index=True)
            continue
        #print (res)
        l = []

        for i in range(nsnp):
            l.append(mr_wald_ratio(x["beta.exposure"].iloc[i], x["beta.outcome"].iloc[i], x["se.exposure"].iloc[i], x["se.outcome"].iloc[i], parameters))

        nom = []
        for method in all_method:
            if method == "mr_ivw":
                l.append(mr_ivw(x["beta.exposure"], x["beta.outcome"], x["se.exposure"], x["se.outcome"], parameters))
            else:
                l.append(mr_egger_regression(x["beta.exposure"], x["beta.outcome"], x["se.exposure"], x["se.outcome"], parameters))
            nom.append("All - " + mr_method_list()[mr_method_list()["obj"] == method]["name"].values[0])

        d = pd.DataFrame({
            "SNP": np.concatenate((x["SNP"].astype(str).values, nom)),
            "b": [y["b"] for y in l],
            "se": [y["se"] for y in l],
            "p": [y["pval"] for y in l],
            "id.exposure": x["id.exposure"].iloc[0],
            "id.outcome": x["id.outcome"].iloc[0],
            "samplesize": x["samplesize.outcome"].iloc[0],
            "outcome": x["outcome"].iloc[0],
            "exposure": x["exposure"].iloc[0]
        })
        res = res.append(d, ignore_index=True)

    res = res[["exposure", "outcome", "id.exposure", "id.outcome", "samplesize", "SNP", "b", "se", "p"]]
    return res

# exp = extract_instruments("ieu-a-2")
# out = extract_outcome_data(snps=exp["SNP"], outcomes=["ieu-a-7"])

# dat = harmonise_data(exp, out)

# res = mr_singlesnp(dat)
# print(res)
