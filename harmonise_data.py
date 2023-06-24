import pandas as pd
import numpy as np
from extract_instruments import *
from extract_outcome_data import *
def harmonise_data(exposure_dat, outcome_dat, action=2):
    if action not in range(1, 4):
        raise ValueError("Action argument must be between 1 and 3 (inclusive).")

    check_required_columns(exposure_dat, "exposure")
    check_required_columns(outcome_dat, "outcome")

    res_tab = pd.merge(outcome_dat, exposure_dat, on="SNP")
    ncombinations = res_tab["id.outcome"].nunique()

    if isinstance(action, int):
        action = [action] * ncombinations
    elif len(action) != ncombinations:
        raise ValueError("Action argument must be of length 1 or number of unique id.outcome values.")

    res_tab = harmonise_cleanup_variables(res_tab)

    if res_tab.shape[0] == 0:
        return res_tab

    d = pd.DataFrame({"id.outcome": res_tab["id.outcome"].unique(), "action": action})
    res_tab = pd.merge(res_tab, d, on="id.outcome")

    combs = res_tab[~res_tab[["id.exposure", "id.outcome"]].duplicated()]

    fix_tab = []
    mr_cols = ["beta.exposure", "beta.outcome", "se.exposure", "se.outcome"]
    
    for i in range(combs.shape[0]):
        x = res_tab[(res_tab["id.exposure"] == combs.iloc[i]["id.exposure"]) & 
                    (res_tab["id.outcome"] == combs.iloc[i]["id.outcome"])]
        
        print(f"Harmonising {x['exposure'].iloc[0]} ({x['id.exposure'].iloc[0]}) and {x['outcome'].iloc[0]} ({x['id.outcome'].iloc[0]})")
        
        x = harmonise(x, 0.08, x["action"].iloc[0])
        x.loc[x[mr_cols].isna().any(axis=1), "mr_keep"] = False
        
        x_attrs = {
            "candidate_variants": exposure_dat["id.exposure"].eq(x["id.exposure"].iloc[0]).sum(),
            "variants_absent_from_reference": exposure_dat["id.exposure"].eq(x["id.exposure"].iloc[0]).sum() - x.shape[0],
            "total_variants": x.shape[0],
            "total_variants_for_mr": x["mr_keep"].sum(),
            "proxy_variants": x["proxy.outcome"].sum() if "proxy.outcome" in x.columns else 0
        }
        
        x_attrs = {key: value for key, value in x_attrs.items() if not pd.isnull(value)}
        x = x.assign(log=x_attrs)
        fix_tab.append(x)
    
    jlog = pd.concat([x.attrs["log"] for x in fix_tab]).reset_index(drop=True)
    fix_tab = pd.concat(fix_tab).reset_index(drop=True)
    fix_tab.attrs["log"] = jlog
    
    if "samplesize.outcome" not in fix_tab.columns:
        fix_tab["samplesize.outcome"] = np.nan
    fix_tab = fix_tab.drop(['log'],axis=1)
    return fix_tab


def harmonise_cleanup_variables(res_tab):
    res_tab["beta.exposure"] = pd.to_numeric(res_tab["beta.exposure"])
    res_tab["beta.outcome"] = pd.to_numeric(res_tab["beta.outcome"])
    res_tab["eaf.exposure"] = pd.to_numeric(res_tab["eaf.exposure"])
    res_tab.loc[res_tab["eaf.outcome"].isin(["NR", "NR "]), "eaf.outcome"] = np.nan
    res_tab["eaf.outcome"] = pd.to_numeric(res_tab["eaf.outcome"])

    res_tab["effect_allele.exposure"] = res_tab["effect_allele.exposure"].str.upper()
    res_tab["other_allele.exposure"] = res_tab["other_allele.exposure"].str.upper()
    res_tab["effect_allele.outcome"] = res_tab["effect_allele.outcome"].str.upper()
    res_tab["other_allele.outcome"] = res_tab["other_allele.outcome"].str.upper()
    res_tab.loc[res_tab["other_allele.outcome"] == "", "other_allele.outcome"] = np.nan

    return res_tab


def harmonise_make_snp_effects_positive(res_tab):
    pos_change = res_tab["beta.exposure"] < 0
    res_tab.loc[pos_change, "eaf.exposure"] = 1 - res_tab.loc[pos_change, "eaf.exposure"]
    res_tab.loc[pos_change, "beta.exposure"] *= -1

    eff_allele_change = res_tab.loc[pos_change, "effect_allele.exposure"]
    oth_allele_change = res_tab.loc[pos_change, "other_allele.exposure"]
    res_tab.loc[pos_change, "effect_allele.exposure"] = oth_allele_change
    res_tab.loc[pos_change, "other_allele.exposure"] = eff_allele_change

    return res_tab
def check_palindromic(A1, A2):
    return ((A1 == "T") & (A2 == "A")) | ((A1 == "A") & (A2 == "T")) | ((A1 == "G") & (A2 == "C")) | ((A1 == "C") & (A2 == "G"))


def flip_alleles(x):
    if isinstance(x, pd.Series):
        return x.apply(flip_alleles)
    else:
        x = str(x).upper()
        flipped_alleles = {'C': 'G', 'G': 'C', 'A': 'T', 'T': 'A'}
        return ''.join(flipped_alleles.get(base, base) for base in x)

def recode_indels_22(A1, A2, B1, B2):
    ncA1 = [len(a) for a in A1]
    ncA2 = [len(a) for a in A2]
    ncB1 = [len(b) for b in B1]
    ncB2 = [len(b) for b in B2]

    i1 = [i for i, (a1, a2, b1, b2) in enumerate(zip(ncA1, ncA2, B1, B2)) if a1 > a2 and b1 == "I" and b2 == "D"]
    for i in i1:
        B1[i] = A1[i]
        B2[i] = A2[i]

    i1 = [i for i, (a1, a2, b1, b2) in enumerate(zip(ncA1, ncA2, B1, B2)) if a1 < a2 and b1 == "I" and b2 == "D"]
    for i in i1:
        B1[i] = A2[i]
        B2[i] = A1[i]

    i1 = [i for i, (a1, a2, b1, b2) in enumerate(zip(ncA1, ncA2, B1, B2)) if a1 > a2 and b1 == "D" and b2 == "I"]
    for i in i1:
        B1[i] = A2[i]
        B2[i] = A1[i]

    i1 = [i for i, (a1, a2, b1, b2) in enumerate(zip(ncA1, ncA2, B1, B2)) if a1 < a2 and b1 == "D" and b2 == "I"]
    for i in i1:
        B1[i] = A1[i]
        B2[i] = A2[i]

    keep = [True] * len(A1)
    for i, (a1, a2, b1, b2) in enumerate(zip(ncA1, ncA2, B1, B2)):
        if a1 > 1 and a1 == a2 and (b1 == "D" or b1 == "I"):
            keep[i] = False
    for i, (a1, a2) in enumerate(zip(A1, A2)):
        if a1 == a2:
            keep[i] = False
    for i, (b1, b2) in enumerate(zip(B1, B2)):
        if b1 == b2:
            keep[i] = False

    return pd.DataFrame({"A1": A1, "A2": A2, "B1": B1, "B2": B2, "keep": keep})


def recode_indels_21(A1, A2, B1):
    ncA1 = [len(a) for a in A1]
    ncA2 = [len(a) for a in A2]
    ncB1 = [len(b) for b in B1]

    B2 = [None] * len(B1)

    i1 = [i for i, (a1, a2, b1) in enumerate(zip(ncA1, ncA2, B1)) if a1 > a2 and b1 == "I"]
    for i in i1:
        B1[i] = A1[i]
        B2[i] = A2[i]

    i1 = [i for i, (a1, a2, b1) in enumerate(zip(ncA1, ncA2, B1)) if a1 < a2 and b1 == "I"]
    for i in i1:
        B1[i] = A2[i]
        B2[i] = A1[i]

    i1 = [i for i, (a1, a2, b1) in enumerate(zip(ncA1, ncA2, B1)) if a1 > a2 and b1 == "D"]
    for i in i1:
        B1[i] = A2[i]
        B2[i] = A1[i]

    i1 = [i for i, (a1, a2, b1) in enumerate(zip(ncA1, ncA2, B1)) if a1 < a2 and b1 == "D"]
    for i in i1:
        B1[i] = A1[i]
        B2[i] = A2[i]

    keep = [True] * len(A1)
    for i, (a1, a2) in enumerate(zip(A1, A2)):
        if a1 == "I" and a2 == "D":
            keep[i] = False
        if a1 == "D" and a2 == "I":
            keep[i] = False
        if a1 == a2 and (ncA1[i] > 1 and ncA1[i] == ncA2[i]) and (B1[i] == "D" or B1[i] == "I"):
            keep[i] = False
        if a1 == a2:
            keep[i] = False

    return pd.DataFrame({"A1": A1, "A2": A2, "B1": B1, "B2": B2, "keep": keep})


def recode_indels_12(A1, B1, B2):
    ncA1 = [len(a) for a in A1]
    ncB1 = [len(b) for b in B1]
    ncB2 = [len(b) for b in B2]

    A2 = [None] * len(A1)

    i1 = [i for i, (b1, b2) in enumerate(zip(B1, B2)) if b1 == "I" and b2 is not None]
    for i in i1:
        A1[i] = B1[i]
        A2[i] = B2[i]

    i1 = [i for i, (b1, b2) in enumerate(zip(B1, B2)) if b1 == "I" and b2 is not None]
    for i in i1:
        A2[i] = B1[i]
        A1[i] = B2[i]

    i1 = [i for i, (b1, b2) in enumerate(zip(B1, B2)) if b1 == "D" and b2 is not None]
    for i in i1:
        A2[i] = B1[i]
        A1[i] = B2[i]

    i1 = [i for i, (b1, b2) in enumerate(zip(B1, B2)) if b1 == "D" and b2 is not None]
    for i in i1:
        A1[i] = B1[i]
        A2[i] = B2[i]

    keep = [True] * len(A1)
    for i, (b1, b2) in enumerate(zip(B1, B2)):
        if b1 == "I" and b2 == "D":
            keep[i] = False
        if b1 == "D" and b2 == "I":
            keep[i] = False
        if b1 == b2 and (ncB1[i] > 1 and ncB1[i] == ncB2[i]) and (A1[i] == "D" or A1[i] == "I"):
            keep[i] = False
        if b1 == b2:
            keep[i] = False

    return pd.DataFrame({"A1": A1, "A2": A2, "B1": B1, "B2": B2, "keep": keep})


def harmonise_22(SNP, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance, action):
    if len(SNP) == 0:
        return pd.DataFrame()

    jlog = {}
    jlog['alleles'] = '2-2'

    indel_index = (A1.str.len() > 1) | (A2.str.len() > 1) | (A1 == 'D') | (A1 == 'I')
    temp = recode_indels_22(A1[indel_index], A2[indel_index], B1[indel_index], B2[indel_index])

    A1.loc[indel_index] = temp['A1']
    A2.loc[indel_index] = temp['A2']
    B1.loc[indel_index] = temp['B1']
    B2.loc[indel_index] = temp['B2']

    status1 = (A1 == B1) & (A2 == B2)
    to_swap = (A1 == B2) & (A2 == B1)
    jlog['switched_alleles'] = to_swap.sum()

    if action == 2:
        Btemp = B1.loc[to_swap].copy()
        B1.loc[to_swap] = B2.loc[to_swap]
        B2.loc[to_swap] = Btemp
        betaB.loc[to_swap] = betaB.loc[to_swap] * -1
        fB.loc[to_swap] = 1 - fB.loc[to_swap]

    status1 = (A1 == B1) & (A2 == B2)
    palindromic = check_palindromic(A1, A2)

    i = (~palindromic) & (~status1)
    B1.loc[i] = flip_alleles(B1.loc[i])
    B2.loc[i] = flip_alleles(B2.loc[i])
    status1 = (A1 == B1) & (A2 == B2)
    jlog['flipped_alleles_basic'] = i.sum()

    i = (~palindromic) & (~status1)
    to_swap = (A1 == B2) & (A2 == B1)
    Btemp = B1.loc[to_swap].copy()
    B1.loc[to_swap] = B2.loc[to_swap]
    B2.loc[to_swap] = Btemp
    betaB.loc[to_swap] = betaB.loc[to_swap] * -1
    fB.loc[to_swap] = 1 - fB.loc[to_swap]

    status1 = (A1 == B1) & (A2 == B2)
    remove = (~status1)
    remove[indel_index][~temp['keep']] = True

    minf = 0.5 - tolerance
    maxf = 0.5 + tolerance
    tempfA = fA.copy()
    tempfB = fB.copy()
    tempfA[tempfA.isna()] = 0.5
    tempfB[tempfB.isna()] = 0.5
    ambiguousA = (tempfA > minf) & (tempfA < maxf)
    ambiguousB = (tempfB > minf) & (tempfB < maxf)

    if action == 2:
        status2 = ((tempfA < 0.5) & (tempfB > 0.5)) | ((tempfA > 0.5) & (tempfB < 0.5)) & palindromic
        to_swap = status2 & ~remove
        betaB.loc[to_swap] = betaB.loc[to_swap] * -1
        fB.loc[to_swap] = 1 - fB.loc[to_swap]
        jlog['flipped_alleles_palindrome'] = to_swap.sum()
    else:
        jlog['flipped_alleles_palindrome'] = 0

    d = pd.DataFrame({
        'SNP': SNP,
        'effect_allele.exposure': A1,
        'other_allele.exposure': A2,
        'effect_allele.outcome': B1,
        'other_allele.outcome': B2,
        'beta.exposure': betaA,
        'beta.outcome': betaB,
        'eaf.exposure': fA,
        'eaf.outcome': fB,
        'remove': remove,
        'palindromic': palindromic,
        'ambiguous': (ambiguousA | ambiguousB) & palindromic
    })
    d.attrs['log'] = jlog
    return d


def harmonise_21(SNP, A1, A2, B1, betaA, betaB, fA, fB, tolerance, action):
    if len(SNP) == 0:
        return pd.DataFrame()

    jlog = {}
    jlog['alleles'] = '2-1'

    n = len(A1)
    B2 = pd.Series([None] * n)
    ambiguous = pd.Series([False] * n)
    palindromic = check_palindromic(A1, A2)
    remove = palindromic

    indel_index = (A1.str.len() > 1) | (A2.str.len() > 1) | (A1 == 'D') | (A1 == 'I')
    temp = recode_indels_21(A1[indel_index], A2[indel_index], B1[indel_index])

    A1.loc[indel_index] = temp['A1']
    A2.loc[indel_index] = temp['A2']
    B1.loc[indel_index] = temp['B1']
    B2.loc[indel_index] = temp['B2']
    remove[indel_index][~temp['keep']] = True

    status1 = A1 == B1
    minf = 0.5 - tolerance
    maxf = 0.5 + tolerance

    tempfA = fA.copy()
    tempfB = fB.copy()
    tempfA[tempfA.isna()] = 0.5
    tempfB[tempfB.isna()] = 0.5

    freq_similar1 = (tempfA < minf) & (tempfB < minf) | (tempfA > maxf) & (tempfB > maxf)
    ambiguous[status1 & ~freq_similar1] = True

    B2[status1] = A2[status1]

    to_swap = A2 == B1
    jlog['switched_alleles'] = to_swap.sum()
    freq_similar2 = (tempfA < minf) & (tempfB > maxf) | (tempfA > maxf) & (tempfB < minf)

    ambiguous[to_swap & ~freq_similar2] = True
    B2[to_swap] = B1[to_swap]
    B1[to_swap] = A1[to_swap]
    betaB[to_swap] = betaB[to_swap] * -1
    fB[to_swap] = 1 - fB[to_swap]

    to_flip = (A1 != B1) & (A2 != B1)
    jlog['flipped_alleles_no_oa'] = to_flip.sum()

    ambiguous[to_flip] = True

    B1[to_flip] = flip_alleles(B1[to_flip])
    status1 = A1 == B1
    B2[status1] = A2[status1]

    to_swap = A2 == B1
    B2[to_swap] = B1[to_swap]
    B1[to_swap] = A1[to_swap]
    betaB[to_swap] = betaB[to_swap] * -1
    fB[to_swap] = 1 - fB[to_swap]

    d = pd.DataFrame({
        'SNP': SNP,
        'effect_allele.exposure': A1,
        'other_allele.exposure': A2,
        'effect_allele.outcome': B1,
        'other_allele.outcome': B2,
        'beta.exposure': betaA,
        'beta.outcome': betaB,
        'eaf.exposure': fA,
        'eaf.outcome': fB,
        'remove': remove,
        'palindromic': palindromic,
        'ambiguous': ambiguous | palindromic
    })
    d.attrs['log'] = jlog
    return d

def harmonise_12(SNP, A1, B1, B2, betaA, betaB, fA, fB, tolerance, action):
    if len(SNP) == 0:
        return pd.DataFrame()
    
    jlog = {}
    jlog['alleles'] = "1-2"
    
    n = len(A1)
    A2 = pd.Series([None] * n)
    ambiguous = pd.Series([False] * n)
    palindromic = check_palindromic(B1, B2)
    remove = palindromic
    
    indel_index = (B1.str.len() > 1) | (B2.str.len() > 1) | (B1 == 'D') | (B1 == 'I')
    temp = recode_indels_21(A1[indel_index], B1[indel_index], B2[indel_index])
    
    A1[indel_index] = temp['A1']
    A2[indel_index] = temp['A2']
    B1[indel_index] = temp['B1']
    B2[indel_index] = temp['B2']
    remove[indel_index][~temp['keep']] = True
    
    status1 = A1 == B1
    minf = 0.5 - tolerance
    maxf = 0.5 + tolerance
    
    tempfA = fA.copy()
    tempfB = fB.copy()
    tempfA[tempfA.isna()] = 0.5
    tempfB[tempfB.isna()] = 0.5
    
    freq_similar1 = (tempfA < minf) & (tempfB < minf) | (tempfA > maxf) & (tempfB > maxf)
    ambiguous[status1 & ~freq_similar1] = True
    
    A2[status1] = B2[status1]
    
    to_swap = A1 == B2
    jlog['switched_alleles'] = to_swap.sum()
    
    freq_similar2 = (tempfA < minf) & (tempfB > maxf) | (tempfA > maxf) & (tempfB < minf)
    
    ambiguous[to_swap & ~freq_similar2] = True
    A2[to_swap] = A1[to_swap]
    A1[to_swap] = B1[to_swap]
    betaA[to_swap] = betaA[to_swap] * -1
    fA[to_swap] = 1 - fA[to_swap]
    
    to_flip = (A1 != B1) & (A1 != B2)
    jlog['flipped_alleles_no_oa'] = to_flip.sum()
    
    ambiguous[to_flip] = True
    
    A1[to_flip] = flip_alleles(A1[to_flip])
    status1 = A1 == B1
    A2[status1] = B2[status1]
    
    to_swap = B2 == A1
    B2[to_swap] = B1[to_swap]
    B1[to_swap] = A1[to_swap]
    betaB[to_swap] = betaB[to_swap] * -1
    fB[to_swap] = 1 - fB[to_swap]
    
    d = pd.DataFrame({
        'SNP': SNP,
        'effect_allele.exposure': A1,
        'other_allele.exposure': A2,
        'effect_allele.outcome': B1,
        'other_allele.outcome': B2,
        'beta.exposure': betaA,
        'beta.outcome': betaB,
        'eaf.exposure': fA,
        'eaf.outcome': fB,
        'remove': remove,
        'palindromic': palindromic,
        'ambiguous': ambiguous | palindromic
    })
    d.attrs['log'] = jlog
    return d


def harmonise_11(SNP, A1, B1, betaA, betaB, fA, fB, tolerance, action):
    if len(SNP) == 0:
        return pd.DataFrame()
    
    jlog = {}
    jlog['alleles'] = "1-1"
    
    n = len(A1)
    A2 = pd.Series([None] * n)
    B2 = pd.Series([None] * n)
    ambiguous = pd.Series([False] * n)
    palindromic = False
    
    status1 = A1 == B1
    remove = ~status1
    
    minf = 0.5 - tolerance
    maxf = 0.5 + tolerance
    
    tempfA = fA.copy()
    tempfB = fB.copy()
    tempfA[tempfA.isna()] = 0.5
    tempfB[tempfB.isna()] = 0.5
    
    freq_similar1 = (tempfA < minf) & (tempfB < minf) | (tempfA > maxf) & (tempfB > maxf)
    ambiguous[status1 & ~freq_similar1] = True
    
    d = pd.DataFrame({
        'SNP': SNP,
        'effect_allele.exposure': A1,
        'other_allele.exposure': A2,
        'effect_allele.outcome': B1,
        'other_allele.outcome': B2,
        'beta.exposure': betaA,
        'beta.outcome': betaB,
        'eaf.exposure': fA,
        'eaf.outcome': fB,
        'remove': remove,
        'palindromic': palindromic,
        'ambiguous': ambiguous | palindromic
    })
    d.attrs['log'] = jlog
    return d


def harmonise(dat, tolerance, action):
    dat['orig_SNP'] = dat['SNP']
    SNP_index = dat.groupby('SNP').cumcount() + 1
    dat['SNP'] = dat['SNP'] + "_" + SNP_index.astype(str)
    
    SNP = dat['SNP']
    A1 = dat['effect_allele.exposure']
    A2 = dat['other_allele.exposure']
    B1 = dat['effect_allele.outcome']
    B2 = dat['other_allele.outcome']
    betaA = dat['beta.exposure']
    betaB = dat['beta.outcome']
    fA = dat['eaf.exposure']
    fB = dat['eaf.outcome']
    
    dat = dat.drop(['effect_allele.exposure', 'other_allele.exposure', 'effect_allele.outcome', 'other_allele.outcome', 'beta.exposure', 'beta.outcome', 'eaf.exposure', 'eaf.outcome'], axis=1)

    i22 = (~A1.isna()) & (~A2.isna()) & (~B1.isna()) & (~B2.isna())
    i21 = (~A1.isna()) & (~A2.isna()) & (~B1.isna()) & (B2.isna())
    i12 = (~A1.isna()) & (A2.isna()) & (~B1.isna()) & (~B2.isna())
    i11 = (~A1.isna()) & (A2.isna()) & (~B1.isna()) & (B2.isna())

    d22 = harmonise_22(SNP[i22], A1[i22], A2[i22], B1[i22], B2[i22], betaA[i22], betaB[i22], fA[i22], fB[i22], tolerance, action)
    d21 = harmonise_21(SNP[i21], A1[i21], A2[i21], B1[i21], betaA[i21], betaB[i21], fA[i21], fB[i21], tolerance, action)
    d12 = harmonise_12(SNP[i12], A1[i12], B1[i12], B2[i12], betaA[i12], betaB[i12], fA[i12], fB[i12], tolerance, action)
    d11 = harmonise_11(SNP[i11], A1[i11], B1[i11], betaA[i11], betaB[i11], fA[i11], fB[i11], tolerance, action)
    #print (d22)
    log_data = []
    if 'log' in d22.attrs:
        log_data.append(pd.DataFrame(d22.attrs['log'],index=[0]))
    if 'log' in d21.attrs:
        log_data.append(pd.DataFrame(d21.attrs['log']))
    if 'log' in d12.attrs:
        log_data.append(pd.DataFrame(d12.attrs['log']))
    if 'log' in d11.attrs:
        log_data.append(pd.DataFrame(d11.attrs['log']))

    jlog = pd.concat(log_data, ignore_index=True)
    jlog = pd.concat([pd.DataFrame({'id.exposure': [dat['id.exposure'][0]], 'id.outcome': [dat['id.outcome'][0]]}), jlog], axis=1)

    d = pd.concat([d21, d22, d12, d11])
    d = pd.merge(d, dat, on='SNP', how='left')
    d['SNP'] = d['orig_SNP']
    d = d.drop('orig_SNP', axis=1)
    d = d.sort_values('id.outcome')
    d['mr_keep'] = True

    if action == 3:
        d.loc[d['palindromic'] | d['remove'] | d['ambiguous'], 'mr_keep'] = False
        if d['palindromic'].any():
            print("Removing the following SNPs for being palindromic:\n", ", ".join(d.loc[d['palindromic'], 'SNP']))
        if d['remove'].any():
            print("Removing the following SNPs for incompatible alleles:\n", ", ".join(d.loc[d['remove'], 'SNP']))
        jlog['incompatible_alleles'] = d['remove'].sum()
        if (d['ambiguous'] & ~d['palindromic']).any():
            print("Removing the following SNPs for having incompatible allele frequencies:\n", ", ".join(d.loc[d['ambiguous'] & ~d['palindromic'], 'SNP']))
        jlog['ambiguous_alleles'] = d['ambiguous'].sum()
    
    if action == 2:
        d.loc[d['remove'] | d['ambiguous'], 'mr_keep'] = False
        if d['remove'].any():
            print("Removing the following SNPs for incompatible alleles:\n", ", ".join(d.loc[d['remove'], 'SNP']))
        jlog['incompatible_alleles'] = d['remove'].sum()
        if d['ambiguous'].any():
            print("Removing the following SNPs for being palindromic with intermediate allele frequencies:\n", ", ".join(d.loc[d['ambiguous'], 'SNP']))
        jlog['ambiguous_alleles'] = d['ambiguous'].sum()
    
    if action == 1:
        d.loc[d['remove'], 'mr_keep'] = False
        if d['remove'].any():
            print("Removing the following SNPs for incompatible alleles:\n", ", ".join(d.loc[d['remove'], 'SNP']))
        jlog['incompatible_alleles'] = d['remove'].sum()

    d.attrs['log'] = jlog
    return d


def check_required_columns(dat, column_type="exposure"):
    required_columns = [
        "SNP",
        f"id.{column_type}",
        f"beta.{column_type}",
        f"se.{column_type}",
        f"effect_allele.{column_type}",
        f"other_allele.{column_type}"
    ]
    missing_columns = [col for col in required_columns if col not in dat.columns]
    if missing_columns:
        raise ValueError(f"The following required columns are missing from {column_type}: {', '.join(missing_columns)}")
    return None

# e = extract_instruments("ieu-a-2")
# o = extract_outcome_data(snps=e["SNP"], outcomes=["ieu-a-7"])

# dat = harmonise_data(e , o)

# print (dat)
