import ieugwaspy as igd
import pandas as pd
import ieugwaspy.backwards as backwards
import ieugwaspy.query as query
import numpy as np
import scipy.stats as stats
from read_outcome_data import *
from extract_instruments import *

def extract_outcome_data(snps, outcomes, proxies=True, rsq=0.8, align_alleles=1, palindromes=1, maf_threshold=0.3, access_token=None, splitsize=10000, proxy_splitsize=500):
    '''
    Extracts outcome data for the given SNPs and outcomes using the ieugwaspy library.

    Args:
    - snps (list): List of SNPs for which outcome data is to be extracted.
    - outcomes (list): List of outcomes for which data is to be extracted.
    - proxies (bool, optional): Flag indicating whether to find proxies for missing SNPs. Default is True.
    - rsq (float, optional): The r-squared threshold for LD proxy search. Default is 0.8.
    - align_alleles (int, optional): Flag indicating whether to align alleles for proxies. Default is 1.
    - palindromes (int, optional): Flag indicating whether to consider palindromic SNPs for proxies. Default is 1.
    - maf_threshold (float, optional): The minor allele frequency threshold for proxies. Default is 0.3.
    - access_token (str, optional): Access token for authentication. Default is None.
    - splitsize (int, optional): Chunk size for splitting the data extraction process. Default is 10000.
    - proxy_splitsize (int, optional): Chunk size for splitting the proxy search process. Default is 500.

    Returns:
    firstpass (pandas.DataFrame): DataFrame containing the extracted outcome data.
    '''
    outcomes = list(set(outcomes))
    snps = list(set(snps))
    firstpass = extract_outcome_data_internal(snps, outcomes, proxies=False, access_token=access_token, splitsize=splitsize)
    proxies = False
    if proxies:
        for i in range(len(outcomes)):
            if firstpass is None:
                missedsnps = snps
            else:
                missedsnps = [snp for snp in snps if snp not in firstpass[firstpass['id.outcome'] == outcomes[i]]['SNP']]
            #print (missedsnps)
            if len(missedsnps) > 0:
                print(f"Finding proxies for {len(missedsnps)} SNPs in outcome {outcomes[i]}")
                temp = extract_outcome_data_internal(missedsnps, [outcomes[i]], proxies=True, rsq=rsq, align_alleles=align_alleles, palindromes=palindromes, maf_threshold=maf_threshold, access_token=access_token, splitsize=proxy_splitsize)
                if temp is not None:
                    firstpass = pd.concat([firstpass, temp], ignore_index=True)

    return firstpass


def extract_outcome_data_internal(snps, outcomes, proxies=True, rsq=0.8, align_alleles=1, palindromes=1, maf_threshold=0.3, access_token=None, splitsize=10000):
    '''
    Internal function to extract outcome data for the given SNPs and outcomes using the ieugwaspy library.

    Args:
    - snps (list): List of SNPs for which outcome data is to be extracted.
    - outcomes (list): List of outcomes for which data is to be extracted.
    - proxies (bool, optional): Flag indicating whether to find proxies for missing SNPs. Default is True.
    - rsq (float, optional): The r-squared threshold for LD proxy search. Default is 0.8.
    - align_alleles (int, optional): Flag indicating whether to align alleles for proxies. Default is 1.
    - palindromes (int, optional): Flag indicating whether to consider palindromic SNPs for proxies. Default is 1.
    - maf_threshold (float, optional): The minor allele frequency threshold for proxies. Default is 0.3.
    - access_token (str, optional): Access token for authentication. Default is None.
    - splitsize (int, optional): Chunk size for splitting the data extraction process. Default is 10000.

    Returns:
    d (pandas.DataFrame): DataFrame containing the extracted outcome data.
    '''
    snps = list(set(snps))
    print(f"Extracting data for {len(snps)} SNP(s) from {len(set(outcomes))} GWAS(s)")
    outcomes = list(set(outcomes))  

    if (proxies == False):
        proxies = 0
    elif (proxies==True):
        proxies = 1
    else:
        raise ValueError("'proxies' argument should be True or False")

    if (len(snps) < splitsize and len(outcomes) < splitsize) or (len(outcomes) < splitsize and len(snps) < splitsize):
        d = query.associations(variantlist=snps, id=outcomes, proxies=proxies, r2=rsq, align_alleles=align_alleles, palindromes=palindromes, maf_threshold=maf_threshold, access_token=access_token)
        d = pd.DataFrame(d)
        #print (d)
        
    elif len(snps) > len(outcomes):
        n = len(snps)
        splits = pd.DataFrame({"snps": snps, "chunk_id": np.repeat(range(1, int(np.ceil(n/splitsize)) + 1), repeats=splitsize)[:n]})
        d = []
        for i in range(len(outcomes)):
            print(f"{i+1} of {len(outcomes)} outcomes")
            temp = splits.groupby("chunk_id").apply(lambda x: x.assign(chunk_id_first=x["chunk_id"].iloc[0]))
            out = query.associations(variants=temp["snps"].tolist(), id=outcomes[i], proxies=proxies, r2=rsq, align_alleles=align_alleles, palindromes=palindromes, maf_threshold=maf_threshold, access_token=access_token)
            if not isinstance(out, pd.DataFrame):
                out = pd.DataFrame()
            d.append(out)

        d = pd.concat(d, ignore_index=True)

    else:
        n = len(outcomes)
        splits = pd.DataFrame({"outcomes": outcomes, "chunk_id": np.repeat(range(1, int(np.ceil(n/splitsize)) + 1), repeats=splitsize)[:n]})
        d = []
        for i in range(len(snps)):
            print(f"{i+1} of {len(snps)} snps")
            temp = splits.groupby("chunk_id").apply(lambda x: x.assign(chunk_id_first=x["chunk_id"].iloc[0]))
            out = query.associations(variants=snps[i], id=temp["outcomes"].tolist(), proxies=proxies, r2=rsq, align_alleles=align_alleles, palindromes=palindromes, maf_threshold=maf_threshold, access_token=access_token)
            if not isinstance(out, pd.DataFrame):
                out = pd.DataFrame()
            d.append(out)

        d = pd.concat(d, ignore_index=True)

    if d.shape[0] == 0 or pd.isnull(d).all().all():
        return None

    d = format_d(d)

    if d.shape[0] > 0:
        d["data_source.outcome"] = "igd"
        return d
    else:
        return None


def cleanup_outcome_data(d):
    '''
    Cleans up the outcome data by replacing invalid values with NaN.

    Args:
    d (pandas.DataFrame): DataFrame containing the outcome data.

    Returns:
    d (pandas.DataFrame): Cleaned DataFrame.
    '''
    d.loc[d["se.outcome"] <= 0, "se.outcome"] = np.nan
    d.loc[(d["eaf.outcome"] <= 0) | (d["eaf.outcome"] >= 1), "eaf.outcome"] = np.nan
    d.loc[d["beta.outcome"] == -9, "beta.outcome"] = np.nan
    return d


def get_se(eff, pval):
    '''
    Calculates the standard error (SE) based on the effect size and p-value.

    Args:
    eff (float): Effect size.
    pval (float): P-value.

    Returns:
    se (float): Standard error.
    '''
    return abs(eff) / abs(stats.norm.ppf(pval / 2))


def format_d(d):
    '''
    Formats the outcome data DataFrame to match the required format.

    Args:
    d (pandas.DataFrame): DataFrame containing the outcome data.

    Returns:
    d (pandas.DataFrame): Formatted DataFrame.
    '''
    d1 = pd.DataFrame({
        "SNP": d["rsid"].astype(str),
        "chr": d["chr"].astype(str),
        "pos": d["position"].astype(str),
        "beta.outcome": pd.to_numeric(d["beta"]),
        "se.outcome": pd.to_numeric(d["se"]),
        "samplesize.outcome": pd.to_numeric(d["n"]),
        "pval.outcome": pd.to_numeric(d["p"]),
        "eaf.outcome": pd.to_numeric(d["eaf"]),
        "effect_allele.outcome": d["ea"].astype(str),
        "other_allele.outcome": d["nea"].astype(str),
        "outcome": d["trait"].astype(str),
        "id.outcome": d["id"].astype(str)
    })

    if "proxy" in d.columns:
        p = pd.DataFrame({
            "proxy.outcome": d["proxy"],
            "target_snp.outcome": d["target_snp"],
            "proxy_snp.outcome": d["proxy_snp"],
            "target_a1.outcome": d["target_a1"],
            "target_a2.outcome": d["target_a2"],
            "proxy_a1.outcome": d["proxy_a1"],
            "proxy_a2.outcome": d["proxy_a2"]
        }, dtype=str)
        d = pd.concat([d1, p], axis=1)

        d = d[~d.duplicated(subset=["proxy_snp.outcome"])]

    else:
        d = d1

    if d.shape[0] == 0:
        print("No matches")
        return d

    d["originalname.outcome"] = d["outcome"]
    if "consortium.outcome" in d.columns and "year.outcome" in d.columns:
        d["outcome.deprecated"] = d["outcome"] + " || " + d["consortium.outcome"] + " || " + d["year.outcome"]
    else:
        d["outcome.deprecated"] = d["outcome"]
    d["outcome"] = d["outcome"] + " || id:" + d["id.outcome"]

    rem = pd.isna(d["beta.outcome"]) & pd.isna(d["pval.outcome"])
    d = d[~rem]

    index = (pd.isna(d["se.outcome"]) | (d["se.outcome"] == 0)) & (~pd.isna(d["beta.outcome"]) & ~pd.isna(d["pval.outcome"]))
    if index.any():
        d.loc[index, "se.outcome"] = get_se(d.loc[index, "beta.outcome"], d.loc[index, "pval.outcome"])

    d = cleanup_outcome_data(d)

    mrcols = ["beta.outcome", "se.outcome", "effect_allele.outcome"]
    d["mr_keep.outcome"] = ~d[mrcols].isna().any(axis=1)
    if (~d["mr_keep.outcome"]).any():
        print("The following SNP(s) are missing required information for the MR tests and will be excluded\n", "\n".join(d[~d["mr_keep.outcome"]]["SNP"]))
    if (~d["mr_keep.outcome"]).all():
        print("None of the provided SNPs can be used for MR analysis, they are missing required information.")

    return d



# d = extract_instruments("ieu-a-2")
# #print(d)
# # #snps = d['rsid']
# # #outcomes = "ieu-a-7"
# outcome_dat = extract_outcome_data(snps=d["SNP"], outcomes=["ieu-a-7"])
# print (outcome_dat)
