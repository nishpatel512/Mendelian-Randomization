import ieugwaspy as igd
import pandas as pd
import ieugwaspy.backwards as backwards
from format_data import *

def extract_instruments(list_ids):
    '''
    Input : Takes a specified outcome as input based on id
    Output : Returns a dataframe with only independent significant associations.
    Description : This function searches for GWAS significant SNPs for a specified set of outcomes.
    '''
    #print (list_ids)
    actual_ids = backwards.legacy_ids([list_ids]) # Converting the list_ids to actual_ids using the legacy_ids function from the backwards module
    #print (actual_ids)
    tophit = igd.tophits(actual_ids)
    #print(tophit)
    d = pd.DataFrame(tophit)
    #print (d)
    d['phenotype'] = d['trait'] + " || id:" + d["id"] # Creating a new column 'phenotype' in the DataFrame by concatenating 'trait' and 'id'
    d = format_data(
        d,
        type="exposure",
        snps=None, header=True,
        phenotype_col="phenotype",
        snp_col="rsid",
        chr_col="chr",
        pos_col="position",
        beta_col="beta",
        se_col="se",
        eaf_col="eaf",
        effect_allele_col="ea",
        other_allele_col="nea",
        pval_col="p",
        samplesize_col="n",
        min_pval=1e-200,
        id_col="id"
    )
    d['data_source.exposure'] = 'igd' # Adding a new column 'data_source.exposure' with the value 'igd' to know the source of the data
    return d

#a = (extract_instruments("ieu-a-2"))
# print (a)
