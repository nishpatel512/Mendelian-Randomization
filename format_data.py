import pandas as pd

def format_data(dat, type="exposure", snps=None, header=True,
                phenotype_col="Phenotype", snp_col="SNP", beta_col="beta",
                se_col="se", eaf_col="eaf", effect_allele_col="effect_allele",
                other_allele_col="other_allele", pval_col="pval", units_col="units",
                ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize",
                gene_col="gene", id_col="id", min_pval=1e-200, z_col="z",
                info_col="info", chr_col="chr", pos_col="pos", log_pval=False):

    all_cols = [phenotype_col, snp_col, beta_col, se_col, eaf_col,
                effect_allele_col, other_allele_col, pval_col, units_col,
                ncase_col, ncontrol_col, samplesize_col, gene_col, id_col,
                z_col, info_col, chr_col, pos_col]
    present_cols = []
    for col in all_cols:
        if col in dat.columns:
            present_cols.append(col)
    if len(present_cols) == 0:
        raise ValueError("None of the specified columns are present")
    #print(dat.columns)
    dat = dat[present_cols]
    #print(dat.columns)

    if snp_col not in dat.columns:
        raise ValueError("SNP column not found")

    dat.loc[:, "SNP"] = dat[snp_col]
    #dat.drop(columns=[snp_col], inplace=True)
    
    dat.loc[:, "SNP"] = dat["SNP"].str.lower().str.replace(r"\s+", "", regex=True)
    dat = dat[~dat["SNP"].isna()]
    #print(dat.columns)
    if snps is not None:
        dat = dat[dat["SNP"].isin(snps)]

    if phenotype_col not in dat.columns:
        print(f"No phenotype name specified, defaulting to '{type}'.")
        dat[type] = type
    else:
        dat[type] = dat[phenotype_col]
        if phenotype_col != type:
            dat = dat.drop(columns=[phenotype_col])

    if log_pval:
        dat[pval_col] = 10**-dat[pval_col]

    def remove_duplicates(x):
        x = x.copy()
        dup = x["SNP"].duplicated()
        if dup.any():
            print(f"Duplicated SNPs present in exposure data for phenotype '{x[type].iloc[0]}'. Just keeping the first instance:")
            print("\n".join(x.loc[dup, "SNP"]))
            x = x[~dup]
        return x

    dat = dat.groupby(type, group_keys=True).apply(remove_duplicates).reset_index(drop=True)
    # snp_col = "SNP" 
    mr_cols_required = [snp_col, beta_col, se_col, effect_allele_col]
    mr_cols_desired = [other_allele_col, eaf_col]
    #print(dat.columns)
    if not all(col in dat.columns for col in mr_cols_required):
        print("The following columns are not present and are required for MR analysis:")
        print("\n".join(col for col in mr_cols_required if col not in dat.columns))
        dat["mr_keep.outcome"] = False
    else:
        dat["mr_keep.outcome"] = True

    if not all(col in dat.columns for col in mr_cols_desired):
        print("The following columns are not present but are helpful for harmonization:")
        print("\n".join(col for col in mr_cols_desired if col not in dat.columns))

    def coerce_numeric_column(x, col_name):
        x[col_name] = pd.to_numeric(x[col_name], errors="coerce")
        return x
    

    dat = dat.groupby(type, group_keys=True).apply(coerce_numeric_column, col_name=beta_col).reset_index(drop=True)
    dat = dat.groupby(type, group_keys=True).apply(coerce_numeric_column, col_name=se_col).reset_index(drop=True)

    if eaf_col in dat.columns:
        dat = dat.groupby(type, group_keys=True).apply(coerce_numeric_column, col_name=eaf_col).reset_index(drop=True)

    dat = dat[dat["mr_keep.outcome"]]

    dat[beta_col] = pd.to_numeric(dat[beta_col], errors="raise")
    dat[se_col] = pd.to_numeric(dat[se_col], errors="raise")

    if eaf_col in dat.columns:
        dat[eaf_col] = pd.to_numeric(dat[eaf_col], errors="raise")

    if pval_col in dat.columns:
        dat[pval_col] = pd.to_numeric(dat[pval_col], errors="raise")
        dat = dat[dat[pval_col] >= min_pval]

    if units_col in dat.columns:
        dat[units_col] = dat[units_col].str.lower()

    if ncase_col in dat.columns and ncontrol_col in dat.columns:
        dat[ncase_col] = pd.to_numeric(dat[ncase_col], errors="coerce")
        dat[ncontrol_col] = pd.to_numeric(dat[ncontrol_col], errors="coerce")
        dat[samplesize_col] = dat[ncase_col] + dat[ncontrol_col]

    if info_col in dat.columns:
        dat[info_col] = pd.to_numeric(dat[info_col], errors="coerce")

    if chr_col in dat.columns and pos_col in dat.columns:
        dat[chr_col] = pd.to_numeric(dat[chr_col], errors="coerce")
        dat[pos_col] = pd.to_numeric(dat[pos_col], errors="coerce")

    if gene_col in dat.columns:
        dat[gene_col] = dat[gene_col].str.upper().str.strip()

    return dat
