import pandas as pd
from format_data import *
import random

def read_exposure_data(filename, clump=False, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta",
                       se_col="se", eaf_col="eaf", effect_allele_col="effect_allele",
                       other_allele_col="other_allele", pval_col="pval", units_col="units",
                       ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize",
                       gene_col="gene", id_col="id", min_pval=1e-200, log_pval=False,
                       chr_col="chr", pos_col="pos"):
    exposure_dat = pd.read_csv(filename, sep=sep)
    exposure_dat = format_data(exposure_dat, "exposure", None, phenotype_col, snp_col, beta_col, se_col, eaf_col,
                               effect_allele_col, other_allele_col, pval_col, units_col, ncase_col, ncontrol_col,
                               samplesize_col, gene_col, id_col, min_pval, log_pval, chr_col, pos_col)
    exposure_dat["data_source.exposure"] = "textfile"
    # if clump:
    #     exposure_dat = clump_data(exposure_dat)
    return exposure_dat


# def clump_data(dat, clump_kb=10000, clump_r2=0.001, clump_p1=1, clump_p2=1, pop="EUR"):
#     pval_column = "pval.exposure"

#     if not isinstance(dat, pd.DataFrame):
#         raise ValueError("Expecting data frame returned from format_data")

#     if "pval.exposure" in dat.columns and "pval.outcome" in dat.columns:
#         print("pval.exposure and pval.outcome columns present. Using pval.exposure for clumping.")
#     elif "pval.exposure" not in dat.columns and "pval.outcome" in dat.columns:
#         print("pval.exposure column not present, using pval.outcome column for clumping.")
#         pval_column = "pval.outcome"
#     elif "pval.exposure" not in dat.columns:
#         print("pval.exposure not present, setting clumping p-value to 0.99 for all variants")
#         dat["pval.exposure"] = 0.99
#     else:
#         pval_column = "pval.exposure"

#     if "id.exposure" not in dat.columns:
#         dat["id.exposure"] = random_string(1)

#     d = pd.DataFrame({"rsid": dat["SNP"], "pval": dat[pval_column], "id": dat["id.exposure"]})
#     out = ld_clump(d, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p1, pop=pop)
#     keep = (dat["SNP"] + dat["id.exposure"]).isin(out["rsid"] + out["id"])
#     return dat[keep]
