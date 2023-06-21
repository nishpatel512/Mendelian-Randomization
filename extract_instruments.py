import ieugwaspy as igd
import pandas as pd
import ieugwaspy.backwards as backwards
from format_data import *

def extract_instruments(list_ids):
    actual_ids = backwards.legacy_ids([list_ids])
    #print (actual_ids)
    tophit = igd.tophits(actual_ids)
    #print(tophit)
    d = pd.DataFrame(tophit)
    #print (d)
    d['phenotype'] = d['trait'] + " || id:" + d["id"]
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
    d['data_source.exposure'] = 'igd'
    return d

#print(extract_instruments("ieu-a-2"))
