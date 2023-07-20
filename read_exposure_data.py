import pandas as pd
from format_data import *
import random

def read_exposure_data(filename, clump=False, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta",
                       se_col="se", eaf_col="eaf", effect_allele_col="effect_allele",
                       other_allele_col="other_allele", pval_col="pval", units_col="units",
                       ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize",
                       gene_col="gene", id_col="id", min_pval=1e-200, log_pval=False,
                       chr_col="chr", pos_col="pos"):
  '''
  Description: This function reads exposure data from a text file and processes it into a formatted DataFrame. The function allows for
  customization of column names and formatting options to ensure the data is appropriately structured. If required, it can also perform
  clumping of the data, which groups variants that are in linkage disequilibrium (LD) with each other based on specified parameters.

  Parameters:
    - filename (str): The path to the text file containing the exposure data.
    - clump (bool, optional): A flag to enable clumping of the data. Default is False.
    - sep (str, optional): The delimiter used in the text file to separate columns. Default is a space (" ").
    - phenotype_col (str, optional): The column name for the phenotype data in the text file. Default is "Phenotype".
    - snp_col (str, optional): The column name for SNP identifiers in the text file. Default is "SNP".
    - beta_col (str, optional): The column name for beta values (effect sizes) of the exposure variable in the text file. Default is "beta".
    - se_col (str, optional): The column name for standard errors of the exposure variable in the text file. Default is "se".
    - eaf_col (str, optional): The column name for the effect allele frequency (EAF) in the text file. Default is "eaf".
    - effect_allele_col (str, optional): The column name for the effect allele in the text file. Default is "effect_allele".
    - other_allele_col (str, optional): The column name for the other allele in the text file. Default is "other_allele".
    - pval_col (str, optional): The column name for p-values in the text file. Default is "pval".
    - units_col (str, optional): The column name for units of measurement in the text file. Default is "units".
    - ncase_col (str, optional): The column name for the number of cases in the text file. Default is "ncase".
    - ncontrol_col (str, optional): The column name for the number of controls in the text file. Default is "ncontrol".
    - samplesize_col (str, optional): The column name for the total sample size in the text file. Default is "samplesize".
    - gene_col (str, optional): The column name for gene information in the text file. Default is "gene".
    - id_col (str, optional): The column name for unique identifiers in the text file. Default is "id".
    - min_pval (float, optional): The minimum p-value threshold for the data. Default is 1e-200.
    - log_pval (bool, optional): A flag indicating if p-values in the text file are given in log-scale. Default is False.
    - chr_col (str, optional): The column name for chromosome information in the text file. Default is "chr".
    - pos_col (str, optional): The column name for base pair position information in the text file. Default is "pos".
  
  Returns:
    A DataFrame containing the formatted exposure data, with an additional column "data_source.exposure" indicating the data source.
  '''

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
