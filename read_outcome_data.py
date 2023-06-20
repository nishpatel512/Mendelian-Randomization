import pandas as pd
from format_data import *
def read_outcome_data(filename, snps=None, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", units_col="units", ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize", gene_col="gene", id_col="id", min_pval=1e-200, log_pval=False, chr_col="chr", pos_col="pos"):
    outcome_dat = pd.read_csv(filename, sep=sep)
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