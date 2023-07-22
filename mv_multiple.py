import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.stats import t,norm,chi2,stats
from mv_harmonise_data import *

def mv_multiple(mvdat, intercept=False, instrument_specific=False, pval_threshold=5e-8, plots=False):
    beta_outcome = mvdat['outcome_beta']
    beta_exposure = mvdat['exposure_beta']
    pval_exposure = mvdat['exposure_pval']
    w = 1 / np.power(mvdat['outcome_se'], 2)

    nexp = beta_exposure.shape[1]
    effs = np.zeros(nexp)
    se = np.zeros(nexp)
    pval = np.zeros(nexp)
    nsnp = np.zeros(nexp)
    p = []

    nom = [str(i) for i in range(beta_exposure.shape[1])]
    nom2 = [mvdat['expname'].loc[mvdat['expname']['id.exposure'] == name, 'exposure'].values[0] if len(mvdat['expname'].loc[mvdat['expname']['id.exposure'] == name]) > 0 else 'Unknown' for name in nom]

    for i in range(nexp):
        index = pval_exposure[:, i] < pval_threshold

        if not intercept:
            if instrument_specific:
                mod = sm.WLS(beta_outcome[index], beta_exposure[index], weights=w[index])
            else:
                mod = sm.WLS(beta_outcome, beta_exposure, weights=w)
        else:
            if instrument_specific:
                mod = sm.WLS(beta_outcome[index], sm.add_constant(beta_exposure[index]), weights=w[index])
            else:
                mod = sm.WLS(beta_outcome, sm.add_constant(beta_exposure), weights=w)

        if instrument_specific and np.sum(index) <= (nexp + int(intercept)):
            effs[i] = np.nan
            se[i] = np.nan
        else:
            results = mod.fit()
            effs[i] = results.params[int(intercept) + i]
            se[i] = results.bse[int(intercept) + i]

        pval[i] = 2 * (1 - stats.norm.cdf(abs(effs[i]) / se[i]))
        nsnp[i] = np.sum(index)

        # Make scatter plot
        d = pd.DataFrame({'outcome': beta_outcome, 'exposure': beta_exposure[:, i]})
        flip = np.sign(d['exposure']) == -1
        d.loc[flip, 'outcome'] *= -1
        d['exposure'] = np.abs(d['exposure'])

        if plots:
            fig, ax = plt.subplots()
            ax.scatter(d.loc[index, 'exposure'], d.loc[index, 'outcome'])
            ax.plot([0, 1], [0, effs[i]], 'r-', lw=2)
            ax.set_xlabel(f'SNP effect on {nom2[i]}')
            ax.set_ylabel('Marginal SNP effect on outcome')
            p.append(fig)

    result = pd.DataFrame({'id.exposure': nom, 'id.outcome': mvdat['outname']['id.outcome'], 'outcome': mvdat['outname']['outcome'], 'nsnp': nsnp, 'b': effs, 'se': se, 'pval': pval})
    result = pd.merge(mvdat['expname'], result)
    out = {'result': result}

    if plots:
        out['plots'] = p

    return out

# id_exposure = ["ieu-a-299","ieu-a-300"]
# exp = mv_extract_exposures(id_exposure)
# out = extract_outcome_data(snps=exp["SNP"], outcomes=["ieu-a-7"])

# mvdat = mv_harmonise_data(exp,out)

# res = mv_multiple(mvdat)
# print (res)


# def mv_multiple(mvdat, intercept=False, instrument_specific=False, pval_threshold=5e-8, plots=False):
#     beta_outcome = mvdat["outcome_beta"]
#     beta_exposure = pd.DataFrame(mvdat["exposure_beta"], columns=mvdat["expname"]["id.exposure"])
#     pval_exposure = pd.DataFrame(mvdat["exposure_pval"], columns=mvdat["expname"]["id.exposure"])
#     w = 1 / mvdat["outcome_se"]**2

#     nexp = beta_exposure.shape[1]
#     effs = np.zeros(nexp)
#     se = np.zeros(nexp)
#     pval = np.zeros(nexp)
#     nsnp = np.zeros(nexp, dtype=int)
#     p = {} if plots else None

#     nom = list(beta_exposure.columns)
#     nom2 = mvdat["expname"].set_index("id.exposure").loc[nom, "exposure"].values

#     for i in range(nexp):
#         index = pval_exposure.iloc[:, i] < pval_threshold

#         if not intercept:
#             if instrument_specific:
#                 mod = np.linalg.lstsq(beta_exposure.loc[index, nom].T, beta_outcome.loc[index], rcond=None)
#             else:
#                 mod = np.linalg.lstsq(beta_exposure[nom].T, beta_outcome, rcond=None)
#         else:
#             if instrument_specific:
#                 mod = np.linalg.lstsq(np.column_stack((np.ones(sum(index)), beta_exposure.loc[index, nom].T)), beta_outcome.loc[index], rcond=None)
#             else:
#                 mod = np.linalg.lstsq(np.column_stack((np.ones(beta_exposure.shape[0]), beta_exposure[nom].T)), beta_outcome, rcond=None)

#         if instrument_specific and sum(index) <= (nexp + int(intercept)):
#             effs[i] = np.nan
#             se[i] = np.nan
#         else:
#             effs[i] = mod[0][int(intercept) + i]
#             se[i] = np.sqrt(np.diag(mod[1]))[int(intercept) + i]

#         pval[i] = 2 * (1 - stats.norm.cdf(abs(effs[i]) / se[i]))

#         nsnp[i] = sum(index)

#         if plots:
#             d = pd.DataFrame({"outcome": beta_outcome.values.flatten(), "exposure": beta_exposure.iloc[:, i].values})
#             flip = np.sign(d["exposure"]) == -1
#             d.loc[flip, "outcome"] *= -1
#             d["exposure"] = abs(d["exposure"])

#             p[i] = plt.scatter(d.loc[index, "exposure"], d.loc[index, "outcome"])
#             plt.plot(d.loc[index, "exposure"], effs[i] * d.loc[index, "exposure"], color="red")
#             plt.xlabel(f"SNP effect on {nom2[i]}")
#             plt.ylabel("Marginal SNP effect on outcome")
#             plt.show()

#     result = pd.DataFrame({
#         "id.exposure": nom,
#         "id.outcome": mvdat["outname"]["id.outcome"],
#         "outcome": mvdat["outname"]["outcome"],
#         "nsnp": nsnp,
#         "b": effs,
#         "se": se,
#         "pval": pval
#     })

#     result = pd.merge(mvdat["expname"], result, on="id.exposure")

#     out = {"result": result}
#     if plots:
#         out["plots"] = p

#     return out

id_exposure = ["ieu-a-299","ieu-a-300"]
exp = mv_extract_exposures(id_exposure)
out = extract_outcome_data(snps=exp["SNP"], outcomes=["ieu-a-7"])

mvdat = mv_harmonise_data(exp,out)

res = mv_multiple(mvdat)
print (res)