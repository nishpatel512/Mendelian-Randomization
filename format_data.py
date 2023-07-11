import pandas as pd
import re
import numpy as np
from scipy.stats import norm
import random
import string

# Function to format the data for analysis
def format_data(dat, type="exposure", snps=None, header=True,
                phenotype_col="Phenotype", snp_col="SNP",
                beta_col="beta", se_col="se", eaf_col="eaf",
                effect_allele_col="effect_allele",
                other_allele_col="other_allele", pval_col="pval",
                units_col="units", ncase_col="ncase",
                ncontrol_col="ncontrol", samplesize_col="samplesize",
                gene_col="gene", id_col="id", min_pval=1e-200,
                z_col="z", info_col="info", chr_col="chr",
                pos_col="pos", log_pval=False):
    # Define the list of all column names
    all_cols = [phenotype_col, snp_col, beta_col, se_col, eaf_col,
                effect_allele_col, other_allele_col, pval_col,
                units_col, ncase_col, ncontrol_col, samplesize_col,
                gene_col, id_col, z_col, info_col, chr_col, pos_col]

    # Check if the specified columns are present in the data
    present_cols = []
    for col in all_cols:
        if col in dat.columns:
            present_cols.append(col)
    if len(present_cols) == 0:
        raise ValueError("None of the specified columns are present")
    
    # Check if the SNP column is present in the data
    if snp_col not in dat.columns:
        raise ValueError("SNP column not found")

    # Clean and preprocess the SNP column
    dat.rename(columns={snp_col: "SNP"}, inplace=True)
    dat["SNP"] = dat["SNP"].str.lower()
    dat.loc[:, "SNP"] = dat["SNP"].str.replace(r"\s+", "", regex=True)

    # Remove rows with missing SNP values
    dat = dat.dropna(subset=["SNP"])

    # Filter data based on specified SNPs
    if snps is not None:
        dat = dat[dat["SNP"].isin(snps)]

    # Handle phenotype column
    if phenotype_col not in dat.columns:
        print(f"No phenotype name specified, defaulting to '{type}'.")
        dat[type] = type
    else:
        dat[type] = dat[phenotype_col]
        if phenotype_col != type:
            dat = dat.drop(columns=[phenotype_col])

    # Convert p-values from log scale if specified
    if log_pval:
        dat["pval"] = 10 ** -dat[pval_col]

    # Function to remove duplicate SNPs within each phenotype
    def remove_duplicates(x):
        dup = x.duplicated(subset="SNP")
        if dup.any():
            print(f"Duplicated SNPs present in exposure data for phenotype '{x[type].iloc[0]}'. Just keeping the first instance:")
            print("\n".join(x[dup]["SNP"]))
            x = x[~dup]
        return x

    # Remove duplicate SNPs within each phenotype
    dat = dat.groupby(type , group_keys = True).apply(remove_duplicates)
    snp_col="SNP"

    # Define required and desired columns for MR analysis
    mr_cols_required = [snp_col, beta_col, se_col, effect_allele_col]
    mr_cols_desired = [other_allele_col, eaf_col]

    # Check if the required columns are present.
    missing_required_cols = set(mr_cols_required) - set(dat.columns)
    if missing_required_cols:
        print("The following columns are not present and are required for MR analysis:")
        print("\n".join(missing_required_cols))
        dat["mr_keep.outcome"] = False
    else:
        dat["mr_keep.outcome"] = True

    missing_desired_cols = set(mr_cols_desired) - set(dat.columns)
    if missing_desired_cols:
        print("The following columns are not present but are helpful for harmonization:")
        print("\n".join(missing_desired_cols))
    
    #Check beta
    i = np.where(dat.columns == beta_col)[0]
    if len(i) > 0:
        dat.columns.values[i[0]] = "beta.outcome"
        if not pd.api.types.is_numeric_dtype(dat["beta.outcome"]):
            print("beta column is not numeric. Coercing...")
            dat["beta.outcome"] = pd.to_numeric(dat["beta.outcome"], errors='coerce')
        index = ~np.isfinite(dat["beta.outcome"])
        index[np.isnan(index)] = True
        dat.loc[index, "beta.outcome"] = np.nan

    # Check se
    se_col_index = dat.columns.get_loc(se_col)
    if not pd.isnull(se_col_index):
        dat.rename(columns={se_col: "se.outcome"}, inplace=True)
        if not pd.api.types.is_numeric_dtype(dat["se.outcome"]):
            print("se column is not numeric. Coercing...")
            dat["se.outcome"] = pd.to_numeric(dat["se.outcome"], errors="coerce")
        index = (~np.isfinite(dat["se.outcome"])) | (dat["se.outcome"] <= 0)
        dat.loc[index.fillna(False), "se.outcome"] = pd.NA

    # Check eaf
    eaf_col_index = dat.columns.get_loc(eaf_col)
    if not pd.isnull(eaf_col_index):
        dat.rename(columns={eaf_col: "eaf.outcome"}, inplace=True)
        if not pd.api.types.is_numeric_dtype(dat["eaf.outcome"]):
            print("eaf column is not numeric. Coercing...")
            dat["eaf.outcome"] = pd.to_numeric(dat["eaf.outcome"], errors="coerce")
        index = (~np.isfinite(dat["eaf.outcome"])) | (dat["eaf.outcome"] <= 0) | (dat["eaf.outcome"] >= 1)
        dat.loc[index.fillna(False), "eaf.outcome"] = pd.NA

     # Check effect_allele
    if effect_allele_col in dat.columns:
        dat = check_column_character(dat, effect_allele_col, "effect_allele.outcome")
        #dat = dat.rename(columns={effect_allele_col: "effect_allele.outcome"})


    # Check other_allele
    if other_allele_col in dat.columns:
        dat = check_column_character(dat, other_allele_col, "other_allele.outcome")

        
    # Check pval
        pval_col_index = dat.columns.get_loc(pval_col)
        if not pd.isnull(pval_col_index):
            dat.rename(columns={pval_col: "pval.outcome"}, inplace=True)
            if not pd.api.types.is_numeric_dtype(dat["pval.outcome"]):
                print("pval column is not numeric. Coercing...")
                dat["pval.outcome"] = pd.to_numeric(dat["pval.outcome"], errors="coerce")
            index = (~np.isfinite(dat["pval.outcome"])) | (dat["pval.outcome"] < 0) | (dat["pval.outcome"] > 1)
            dat.loc[index.fillna(False), "pval.outcome"] = pd.NA
            index = dat["pval.outcome"] < min_pval
            dat.loc[index.fillna(False), "pval.outcome"] = min_pval

            dat["pval_origin.outcome"] = "reported"
            if any(dat["pval.outcome"].isna()):
                if ("beta.outcome" in dat.columns) and ("se.outcome" in dat.columns):
                    index = dat["pval.outcome"].isna()
                    dat.loc[index.fillna(False), "pval.outcome"] = np.abs(dat.loc[index.fillna(False), "beta.outcome"]) / dat.loc[index.fillna(False), "se.outcome"]
                    dat.loc[index.fillna(False), "pval_origin.outcome"] = "inferred"
        
        # If no pval column then create it from beta and se if available
        if ("beta.outcome" in dat.columns) and ("se.outcome" in dat.columns) and ("pval.outcome" not in dat.columns):
            print("Inferring p-values")
            dat["pval.outcome"] = norm.sf(abs(dat["beta.outcome"]) / dat["se.outcome"]) * 2
            dat["pval_origin.outcome"] = "inferred"

        if ncase_col in dat.columns:
            dat.rename(columns={ncase_col: "ncase.outcome"}, inplace=True)
            if not pd.api.types.is_numeric_dtype(dat["ncase.outcome"]):
                print(ncase_col, " column is not numeric")
                dat["ncase.outcome"] = pd.to_numeric(dat["ncase.outcome"], errors="coerce")

        if ncontrol_col in dat.columns:
            dat.rename(columns={ncontrol_col: "ncontrol.outcome"}, inplace=True)
            if not pd.api.types.is_numeric_dtype(dat["ncontrol.outcome"]):
                print(ncontrol_col, " column is not numeric")
                dat["ncontrol.outcome"] = pd.to_numeric(dat["ncontrol.outcome"], errors="coerce")

        if samplesize_col in dat.columns:
            dat.rename(columns={samplesize_col: "samplesize.outcome"}, inplace=True)
            if not pd.api.types.is_numeric_dtype(dat["samplesize.outcome"]):
                print(samplesize_col, " column is not numeric")
                dat["samplesize.outcome"] = pd.to_numeric(dat["samplesize.outcome"], errors="coerce")

            if ("ncontrol.outcome" in dat.columns) and ("ncase.outcome" in dat.columns):
                index = dat["samplesize.outcome"].isna() & ~dat["ncase.outcome"].isna() & ~dat["ncontrol.outcome"].isna()
                if any(index):
                    print("Generating sample size from ncase and ncontrol")
                    dat.loc[index.fillna(False), "samplesize.outcome"] = dat.loc[index.fillna(False), "ncase.outcome"] + dat.loc[index.fillna(False), "ncontrol.outcome"]
        elif ("ncontrol.outcome" in dat.columns) and ("ncase.outcome" in dat.columns):
            print("Generating sample size from ncase and ncontrol")
            dat["samplesize.outcome"] = dat["ncase.outcome"] + dat["ncontrol.outcome"]

        if gene_col in dat.columns:
            dat.rename(columns={gene_col: "gene.outcome"}, inplace=True)

        if info_col in dat.columns:
            dat.rename(columns={info_col: "info.outcome"}, inplace=True)

        if z_col in dat.columns:
            dat.rename(columns={z_col: "z.outcome"}, inplace=True)

        if chr_col in dat.columns:
            dat.rename(columns={chr_col: "chr.outcome"}, inplace=True)

        if pos_col in dat.columns:
            dat.rename(columns={pos_col: "pos.outcome"}, inplace=True)

        if units_col in dat.columns:
            dat.rename(columns={units_col: "units.outcome"}, inplace=True)
            dat["units.outcome_dat"] = dat["units.outcome"].astype(str)
            temp = check_units(dat, type, "units.outcome")
            if any(temp["ph"]):
                dat[type] = dat[type] + " (" + dat["units.outcome"] + ")"

        # Create id column
        if id_col in dat.columns:
            dat.rename(columns={id_col: "id.outcome"}, inplace=True)
            dat["id.outcome"] = dat["id.outcome"].astype(str)
        else:
            dat["id.outcome"] = create_ids(dat[type])

        if any(dat["mr_keep.outcome"]):
            mrcols = ["SNP", "beta.outcome", "se.outcome", "effect_allele.outcome"]
            mrcols_present = [col for col in mrcols if col in dat.columns]
            dat["mr_keep.outcome"] = dat["mr_keep.outcome"] & dat[mrcols_present].apply(lambda x: ~x.isna()).all(axis=1)
            if any(~dat["mr_keep.outcome"]):
                print("The following SNP(s) are missing required information for the MR tests and will be excluded:")
                print("\n".join(dat.loc[~dat["mr_keep.outcome"], "SNP"].tolist()))

        if not all(dat["mr_keep.outcome"]):
            print("None of the provided SNPs can be used for MR analysis; they are missing required information.")

        # Add in missing MR cols
        missing_mr_cols = ["SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome"]
        for col in missing_mr_cols:
            if col not in dat.columns:
                dat[col] = pd.NA

        dat.columns = [col.replace("outcome", type) for col in dat.columns]
        dat.reset_index(drop=True, inplace=True)
        return dat
    
def check_units(x, id, col):
    """
    Check if more than one type of unit is specified for each group.
    """
    temp = x.groupby(id).apply(lambda x1: pd.DataFrame({'ph': [False]} if len(x1[col].unique()) <= 1 else {'ph': [True]}))
    if any(temp['ph']):
        print("More than one type of unit specified for", x[id].iloc[0])
        x1 = x.copy()
        x1['ph'] = temp['ph']
        return x1
    return temp


def random_string(n=1, length=6):
    """
    Generate random strings of alphanumeric characters.
    """
    random_strings = []
    for _ in range(n):
        random_string = ''.join(random.choices(string.digits + string.ascii_letters, k=length))
        random_strings.append(random_string)
    return random_strings

def create_ids(x):
    """
    Create categorical IDs and replace them with random strings.
    """
    a = x.astype('category')
    a = a.cat.set_categories(random_string(len(a.cat.categories)))
    a = a.astype('str')
    return a.tolist()

def check_column_character(dat, col_name, new_col_name):
    """
    Check if a column contains character data and coerce it if necessary.
    Convert the column values to uppercase and filter out invalid values.
    """
    if not pd.api.types.is_string_dtype(dat[col_name]):
        print(f"{col_name} column is not character data. Coercing...")
        dat[new_col_name] = dat[col_name].astype(str)
    else:
        dat[new_col_name] = dat[col_name]

    dat[new_col_name] = dat[new_col_name].str.upper()
    index = ~dat[new_col_name].str.match(r"^[ACTG]+$") & ~dat[new_col_name].isin(["D", "I"])
    if any(index):
        print("effect_allele column has some values that are not A/C/T/G or an indel comprising only these characters or D/I. These SNPs will be excluded.")
        dat.loc[index, new_col_name] = np.nan
        dat.loc[index, "mr_keep.outcome"] = False
    dat = dat.drop(columns=[col_name])
    
    return dat
