import pandas as pd
import numpy as np
from mv_extract_exposures import *
from extract_outcome_data import *
from harmonise_data import *

def mv_harmonise_data(exposure_dat, outcome_dat, harmonise_strictness=2):
    assert all(col in exposure_dat.columns for col in ["SNP", "id.exposure", "exposure", "effect_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure"]), "Missing required columns in exposure_dat"
    nexp = len(exposure_dat["id.exposure"].unique())
    assert nexp > 1, "Number of exposures (nexp) must be greater than 1"
    tab = exposure_dat["SNP"].value_counts()
    keepsnp = tab[tab == nexp].index.tolist()
    exposure_dat = exposure_dat[exposure_dat["SNP"].isin(keepsnp)].copy()

    exposure_beta = exposure_dat.pivot(index="SNP", columns="id.exposure", values="beta.exposure").reset_index()
    exposure_pval = exposure_dat.pivot(index="SNP", columns="id.exposure", values="pval.exposure").reset_index()
    exposure_se = exposure_dat.pivot(index="SNP", columns="id.exposure", values="se.exposure").reset_index()

    dat = harmonise_data(exposure_dat[exposure_dat["id.exposure"] == exposure_dat["id.exposure"].iloc[0]].copy(), outcome_dat, action=harmonise_strictness)
    dat = dat[dat["mr_keep"]]
    dat["SNP"] = dat["SNP"].astype(str)

    exposure_beta = exposure_beta[exposure_beta["SNP"].isin(dat["SNP"])]
    exposure_pval = exposure_pval[exposure_pval["SNP"].isin(dat["SNP"])]
    exposure_se = exposure_se[exposure_se["SNP"].isin(dat["SNP"])]

    exposure_beta = exposure_beta.drop("SNP", axis=1).to_numpy()
    exposure_pval = exposure_pval.drop("SNP", axis=1).to_numpy()
    exposure_se = exposure_se.drop("SNP", axis=1).to_numpy()

    exposure_beta = np.asarray(exposure_beta, dtype=float)
    exposure_pval = np.asarray(exposure_pval, dtype=float)
    exposure_se = np.asarray(exposure_se, dtype=float)

    exposure_beta = np.nan_to_num(exposure_beta)
    exposure_pval = np.nan_to_num(exposure_pval)
    exposure_se = np.nan_to_num(exposure_se)


    outcome_beta = dat["beta.outcome"].to_numpy()
    outcome_se = dat["se.outcome"].to_numpy()
    outcome_pval = dat["pval.outcome"].to_numpy()

    expname = exposure_dat[~exposure_dat.duplicated("id.exposure")][["id.exposure", "exposure"]].copy()
    outname = outcome_dat[~outcome_dat.duplicated("id.outcome")][["id.outcome", "outcome"]].copy()

    return {"exposure_beta": exposure_beta,
            "exposure_pval": exposure_pval,
            "exposure_se": exposure_se,
            "outcome_beta": outcome_beta,
            "outcome_pval": outcome_pval,
            "outcome_se": outcome_se,
            "expname": expname,
            "outname": outname}
# id_exposure = ["ieu-a-299","ieu-a-300"]
# exp = mv_extract_exposures(id_exposure)
# out = extract_outcome_data(snps=exp["SNP"], outcomes=["ieu-a-7"])

# dat = mv_harmonise_data(exp,out)
# print (dat)