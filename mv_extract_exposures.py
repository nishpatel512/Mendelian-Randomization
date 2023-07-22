import pandas as pd
from extract_instruments import *
from extract_outcome_data import *
from harmonise_data import *
from mv_format_data import *
def convert_outcome_to_exposure(outcome_dat):
    id = outcome_dat[~outcome_dat.duplicated('outcome')][['outcome', 'id.outcome']]
    exposure_dat = mv_format_data(
        outcome_dat,
        beta_col="beta.outcome",
        se_col="se.outcome",
        pval_col="pval.outcome",
        phenotype_col="outcome",
        effect_allele_col="effect_allele.outcome",
        other_allele_col="other_allele.outcome",
        eaf_col="eaf.outcome",
        units_col="units.outcome"
    )
    #print (exposure_dat)
    exposure_dat = exposure_dat.merge(id, left_on='exposure', right_on='outcome')
    exposure_dat = exposure_dat.drop(columns=['id.exposure','samplesize.exposure','originalname.exposure','exposure.deprecated','data_source.exposure','outcome'])
    exposure_dat.rename(columns={'id.outcome': 'id.exposure'}, inplace=True)
    return exposure_dat

def mv_extract_exposures(id_exposure, clump_r2=0.001, clump_kb=10000, harmonise_strictness=2, access_token=None, find_proxies=True, force_server=False, pval_threshold=5e-8, pop="EUR"):
    assert len(id_exposure) > 1, "Length of id_exposure should be greater than 1"
    L=[]
    for id in id_exposure:
        # print(id)
        exp_dat = extract_instruments(id)
        L.append(exp_dat)
    exposure_dat = pd.concat(L)

    # Get best instruments for each exposure
    temp = exposure_dat.copy()
    temp['id.exposure'] = 1
    temp = temp.sort_values('pval.exposure', ascending=True)
    temp = temp.drop_duplicates('SNP')
    #temp = clump_data(temp, clump_p1=pval_threshold, clump_r2=clump_r2, clump_kb=clump_kb, pop=pop)
    exposure_dat = exposure_dat[exposure_dat['SNP'].isin(temp['SNP'])]

    # Get effects of each instrument from each exposure
    d1 = extract_outcome_data(snps = exposure_dat['SNP'], outcomes = id_exposure)
    # print (d1)
    assert len(d1['id.outcome'].unique()) == len(id_exposure), "Number of unique IDs in outcome data does not match number of unique id_exposure values"
    d1 = d1[d1['mr_keep.outcome']]
    d2 = d1[d1['id.outcome'] != id_exposure[0]]
    # Filter 'd1' based on the condition 'id.outcome == id_exposure[1]'
    subset_d1 = d1[d1['id.outcome'] == id_exposure[0]]
    # Create a new DataFrame 'd1' by assigning 'subset_d1' to it
    d1 = subset_d1.copy()
    #print (d1)
    d1 = convert_outcome_to_exposure(d1)
    #print (d1)
    #print (d2)
    # Harmonise against the first id
    d = harmonise_data(d1, d2, action=harmonise_strictness)
    #print (d)
    # Only keep SNPs that are present in all
    tab = d['SNP'].value_counts()
    keepsnps = tab[tab == len(id_exposure) - 1].index
    d = d[d['SNP'].isin(keepsnps)]

    # Reshape exposures
    dh1 = d[d['id.outcome'] == id_exposure[1]][['SNP', 'exposure', 'id.exposure', 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure']]
    #print (dh1)
    dh2 = d[['SNP', 'outcome', 'id.outcome', 'effect_allele.outcome', 'other_allele.outcome', 'eaf.outcome', 'beta.outcome', 'se.outcome', 'pval.outcome']]
    dh2.columns = [col.replace('outcome', 'exposure') for col in dh2.columns]
    dh = pd.concat([dh1, dh2])

    return dh

# id_exposure = ["ieu-a-299","ieu-a-300"]
# res = mv_extract_exposures(id_exposure)
# print (res)

