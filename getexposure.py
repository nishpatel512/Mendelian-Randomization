import ieugwaspy as igd
import pandas as pd

data = igd.gwasinfo()
df = pd.DataFrame.from_dict(data, orient='index')
df.head()


#getting exposure instruments from GWAS
def extract_instruments(outcomes, p1=5e-8, clump=True, p2=5e-8, r2=0.001, kb=10000, access_token=igd.check_access_token(), force_server=False):

  # .Deprecated("igd::tophits()")
  outcomes = igd.legacy_ids(unique(outcomes))

  d = igd.tophits(outcomes, pval=p1, clump=clump, r2=r2, kb=kb, force_server=force_server, access_token=access_token)

  # d$phenotype.deprecated <- paste0(d$trait, " || ", d$consortium, " || ", d$year, " || ", d$unit)
  if nrow(d) == 0:
    return None

  d["phenotype"] = d["trait"] + " || id:" + d["id"]
  d = format_data(
    d,
    type="exposure",
    snps=None,
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
  d["data_source.exposure"] = "igd"
  return d