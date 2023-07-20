import pandas as pd
from format_data import *
def read_outcome_data(filename, snps=None, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", units_col="units", ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize", gene_col="gene", id_col="id", min_pval=1e-200, log_pval=False, chr_col="chr", pos_col="pos"):
    outcome_dat = pd.read_csv(filename, sep=sep)
'''
  Description: This function reads outcome data from a text file and formats it for Mendelian Randomization (MR) analysis.

  Parameters:
  - filename (str): The path to the text file containing the outcome data.
  - snps (list or None): A list of SNP IDs to include in the analysis. If None, all SNPs in the file will be used.
  - sep (str): The delimiter used in the text file.
  - phenotype_col (str): The column name for the outcome phenotype.
  - snp_col (str): The column name for SNP IDs.
  - beta_col (str): The column name for beta values.
  - se_col (str): The column name for standard errors of beta values.
  - eaf_col (str): The column name for effect allele frequencies.
  - effect_allele_col (str): The column name for the effect allele.
  - other_allele_col (str): The column name for the other allele.
  - pval_col (str): The column name for p-values.
  - units_col (str): The column name for units of measurement.
  - ncase_col (str): The column name for the number of cases (for binary outcomes).
  - ncontrol_col (str): The column name for the number of controls (for binary outcomes).
  - samplesize_col (str): The column name for the total sample size.
  - gene_col (str): The column name for gene names or identifiers (optional).
  - id_col (str): The column name for individual or study identifiers (optional).
  - min_pval (float): The minimum p-value to use for clumping (default is 1e-200).
  - log_pval (bool): If True, assumes p-values are in log scale (default is False).
  - chr_col (str): The column name for chromosome information (optional).
  - pos_col (str): The column name for base pair positions (optional).

  Returns:
  outcome_dat (pd.DataFrame): A formatted DataFrame containing the outcome data for MR analysis.
'''
    outcome_dat = format_data(
        outcome_dat,
        type="outcome",
        snps=snps,
        phenotype_col=phenotype_col,
        snp_col=snp_col,
        beta_col=beta_col,
        se_col=se_col,
        eaf_col=eaf_col,
        effect_allele_col=effect_allele_col,
        other_allele_col=other_allele_col,
        pval_col=pval_col,
        units_col=units_col,
        ncase_col=ncase_col,
        ncontrol_col=ncontrol_col,
        samplesize_col=samplesize_col,
        gene_col=gene_col,
        id_col=id_col,
        min_pval=min_pval,
        log_pval=log_pval,
        chr_col=chr_col,
        pos_col=pos_col
    )
    outcome_dat['data_source.outcome'] = "textfile"
    return outcome_dat
