# Mandelian-Randomization
# Introduction
Mendelian randomization is a method of using measured variation in genes of known function to examine the causal effect of a modifiable exposure on disease in observational studies.
We are making use of python package called MRPackage which has all the functionalities to perfom MR analysis using data from local computer or data from GWAS.

# Methods:

  ```python
  available_outcomes()
  '''
    Description: Retrieves GWAS (Genome-Wide Association Study) information using the ieugwaspy library and returns
    it as a Pandas DataFrame.
    Returns: df (pandas.DataFrame): DataFrame containing GWAS information, with outcomes as the index.
  '''
  ```

  ```python
  blank_plot()
  '''
    Description: Create a blank plot with a centered text message.
    Parameters: message (str): The message to be displayed in the plot.
    Returns: matplotlib.axes._subplots.AxesSubplot: The blank plot with the message.
  '''
  ```

  ```python
  check_access_token()
  '''
    Description: Retrieves the access token for Google API authentication.
    Returns:str: The access token for API authentication.
  '''
  ```

  ```python
  extract_instruments(list_ids)
  '''
    Description : This function searches for GWAS significant SNPs for a specified set of outcomes.
    Parameters : Takes a specified outcome as input based on id
    Returns : Returns a dataframe with only independent significant associations.
  '''
  ```
  ```python
  extract_outcome_data()
  '''
    Description: Extracts outcome data for the given SNPs and outcomes using the ieugwaspy library.
    Paramters:
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
  ```

  ```python
  format_data(dat, type="exposure", snps=None, header=True, phenotype_col="Phenotype", snp_col="SNP",beta_col="beta",
  se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval",
  units_col="units", ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize", gene_col="gene", id_col="id",
  min_pval=1e-200, z_col="z", info_col="info", chr_col="chr",pos_col="pos", log_pval=False)
  '''
    Description: This function is used to format and preprocess genetic association data for Mendelian randomization (MR) analysis. It takes a
    pandas DataFrame containing genetic association data, which typically includes information about the exposure (e.g., risk factor) and
    outcome (e.g., disease) variables. The function performs several data cleaning and manipulation steps to ensure the data is in a suitable
    format for conducting MR analysis.

    Parameters:
    - dat: pandas DataFrame: Input data containing genetic association information, such as beta coefficients, standard errors, p-values,
      allele frequencies, sample sizes, etc.
    - type: str (optional, default: "exposure"): Specifies whether the data represents the "exposure" or "outcome" variable.
    - snps: list (optional, default: None): A list of specific SNPs to be included in the analysis. If provided, the function filters
      the data to include only the specified SNPs.
    - header: bool (optional, default: True): Specifies whether the input DataFrame has a header row.
    - phenotype_col, snp_col, beta_col, se_col, eaf_col, effect_allele_col, other_allele_col, pval_col, units_col, ncase_col, ncontrol_col,
      samplesize_col, gene_col, id_col, min_pval, z_col, info_col, chr_col, pos_col: str (optional): Column names corresponding to specific
      genetic association data in the input DataFrame.
    - log_pval: bool (optional, default: False): Specifies whether p-values are provided in logarithmic scale (log10).

    Returns:
    - dat: pandas DataFrame
    The formatted and preprocessed data, ready for Mendelian randomization analysis.
  '''
  ```

  ```python
  harmonise_data(exposure_dat, outcome_dat, action=2)
  '''
    Description: In order to perform MR the effect of a SNP on an outcome and exposure must be harmonised to be relative to the same allele.
    Parameters: Takes formated input data of Exposure and Outcome data
    Returns: Data frame with harmonised effects and alleles
  '''
  ```

  ```python
  mr(dat, parameters=default_parameters(), method_list=mr_method_list())
  '''
  Description: The harmonise_data function is used to perform harmonization of genetic effect sizes and alleles for Mendelian Randomization (MR) analysis.
  It ensures that the effect of a Single Nucleotide Polymorphism (SNP) on both the exposure and outcome variables is measured relative to the same allele.
  Parameters:
  -exposure_data (pandas.DataFrame): DataFrame containing the exposure data with columns such as "SNP," "effect_size," "standard_error," "effect_allele," and "other_allele."
  -outcome_data (pandas.DataFrame): DataFrame containing the outcome data with columns such as "SNP," "effect_size," "standard_error," "effect_allele," and "other_allele."
  Returns:
  harmonized_data (pandas.DataFrame): Data frame with harmonized genetic effect sizes and alleles, where the effect sizes are relative to the
  same allele for both exposure and outcome data.
  '''
  ```

  ```python
  mr_egger_regression(b_exp, b_out, se_exp, se_out, parameters)
  '''
    Description: The harmonise_data function is used to perform harmonization of genetic effect sizes and alleles for Mendelian Randomization (MR) analysis.
    It ensures that the effect of a Single Nucleotide Polymorphism (SNP) on both the exposure and outcome variables is measured relative to the same allele.
    Parameters:
    -exposure_data (pandas.DataFrame): DataFrame containing the exposure data with columns such as "SNP," "effect_size," "standard_error," "effect_allele," and "other_allele."
    -outcome_data (pandas.DataFrame): DataFrame containing the outcome data with columns such as "SNP," "effect_size," "standard_error," "effect_allele," and "other_allele."
    Returns:
    harmonized_data (pandas.DataFrame): Data frame with harmonized genetic effect sizes and alleles, where the effect sizes are relative to the same
    allele for both exposure and outcome 
    data.
  '''
  ```

  ```python
  mr_egger_regression_bootstrap(b_exp, b_out, se_exp, se_out, parameters)
  '''
    Description: This function performs Mendelian Randomization (MR) using Egger regression with bootstrap for statistical inference.
    It estimates causal effects of an exposure variable on an outcome variable by harmonizing genetic effect sizes and alleles.
    Parameters:
    -b_exp (array-like): Beta values of the exposure variable.
    -b_out (array-like): Beta values of the outcome variable.
    -se_exp (array-like): Standard errors of the exposure variable.
    -se_out (array-like): Standard errors of the outcome variable.
    -parameters (dict): Dictionary containing additional parameters for the bootstrap process. Should contain the "nboot" parameter
     indicating the number of bootstrap iterations.
    Returns:
    result (dict): Dictionary containing the MR estimates, standard errors, and p-values for the causal effect, as well as for the instrument strength.
  '''
  ```

  ```python
  mr_forest_plot(singlesnp_results, exponentiate=False)
  '''
    Description: This function creates forest plots for Mendelian Randomization (MR) analysis results based on the output from the mr_singlesnp function.
    Parameters:
    -singlesnp_results (pandas.DataFrame): DataFrame containing MR analysis results from the mr_singlesnp function.
    -exponentiate (bool, optional): Whether to exponentiate the effect sizes. Defaults to False.
    Returns:
    res (list): A list of matplotlib.pyplot objects containing the generated forest plots.
  '''
  ```

  ```python
  mr_ivw(b_exp, b_out, se_exp, se_out, parameters)
  '''
    Description: This function performs Mendelian Randomization (MR) using the Inverse Variance Weighted (IVW) method.
    Parameters:
    -b_exp (array-like): Beta values of the exposure variable.
    -b_out (array-like): Beta values of the outcome variable.
    -se_exp (array-like): Standard errors of the exposure variable.
    -se_out (array-like): Standard errors of the outcome variable.
    -parameters (dict): Additional parameters (Not used in this method).
    Returns:
    result (dict): Dictionary containing the IVW MR estimates, standard errors, and p-values for the causal effect.
  '''
  ```

  ```python
  mr_ivw_fe(b_exp, b_out, se_exp, se_out, parameters)
  '''
    Description: This function performs Mendelian Randomization (MR) using the Inverse Variance Weighted (IVW) method with
    fixed-effects standard error estimation.
    Parameters:
    -b_exp (array-like): Beta values of the exposure variable.
    -b_out (array-like): Beta values of the outcome variable.
    -se_exp (array-like): Standard errors of the exposure variable.
    -se_out (array-like): Standard errors of the outcome variable.
    -parameters (dict): Additional parameters (Not used in this method).
    Returns:
    result (dict): Dictionary containing the IVW MR estimates, fixed-effects standard errors, and p-values for the causal effect.
  '''
  ```

  ```python
  mr_ivw_mre(b_exp, b_out, se_exp, se_out, parameters)
  '''
    Description: This function performs Mendelian Randomization (MR) using the Inverse Variance Weighted (IVW) method
    with standard errors that account for multi-variant effects.
    Parameters:
    -b_exp (array-like): Beta values of the exposure variable.
    -b_out (array-like): Beta values of the outcome variable.
    -se_exp (array-like): Standard errors of the exposure variable.
    -se_out (array-like): Standard errors of the outcome variable.
    -parameters (dict): Additional parameters (Not used in this method).
    Returns:
    result (dict): Dictionary containing the IVW MR estimates, standard errors with multi-variant effects, and p-values for the causal effect.
  '''
  ```

  ```python
  mr_penalised_weighted_median(b_exp, b_out, se_exp, se_out, parameters)
    '''
    Description: This function performs Mendelian Randomization (MR) using the penalized weighted median method.
    Parameters:
    -b_exp (array-like): Beta values of the exposure variable.
    -b_out (array-like): Beta values of the outcome variable.
    -se_exp (array-like): Standard errors of the exposure variable.
    -se_out (array-like): Standard errors of the outcome variable.
    -parameters (dict): Dictionary containing additional parameters for the penalized weighted median method, including "penk" and "nboot".
    Returns:
    result (dict): Dictionary containing the MR estimates, standard errors, and p-values for the causal effect.
    '''
  ```

  ```python
  mr_raps(b_exp, b_out, se_exp, se_out, parameters)
    '''
    Description: This function performs Mendelian Randomization (MR) using the robust adjusted profile score (MR-RAPS) method.
    Parameters:
    -b_exp (array-like): Beta values of the exposure variable.
    -b_out (array-like): Beta values of the outcome variable.
    -se_exp (array-like): Standard errors of the exposure variable.
    -se_out (array-like): Standard errors of the outcome variable.
    -parameters (dict): Dictionary containing additional parameters for the MR-RAPS method, including "over_dispersion", "loss_function",
     and "shrinkage".
    Returns:
    result (dict): Dictionary containing the MR estimates, robust standard errors, and p-values for the causal effect.
    '''
  ```

  ```python
  mr_scatter_plot(mr_results, dat)
    '''
    Description: This function creates scatter plots to compare SNP effects on the exposure variable with SNP effects on the
    outcome variable for Mendelian Randomization (MR) analysis results.
    Parameters:
    -mr_results (pandas.DataFrame): DataFrame containing MR analysis results.
    -dat (pandas.DataFrame): DataFrame containing harmonized data.
    Returns:
    res (list): A list of matplotlib.pyplot objects containing the generated scatter plots.
    '''
  ```

  ```python
  mr_sign(b_exp, b_out, se_exp, se_out, parameters)
    '''
    Description: This function calculates Mendelian Randomization (MR) estimates using the sign instrument approach.
    Parameters:
    -b_exp (array-like): Beta values of the exposure variable.
    -b_out (array-like): Beta values of the outcome variable.
    -se_exp (array-like): Standard errors of the exposure variable.
    -se_out (array-like): Standard errors of the outcome variable.
    -parameters (dict): Additional parameters (Not used in this method).
    Returns:
    result (dict): Dictionary containing the MR estimates, standard errors, and p-values for the causal effect.
    '''
  ```

  ```python
 mr_simple_median(b_exp, b_out, se_exp, se_out, parameters)
    '''
    Description: This function calculates Mendelian Randomization (MR) estimates using the simple median method.
    Parameters:
    -b_exp (array-like): Beta values of the exposure variable.
    -b_out (array-like): Beta values of the outcome variable.
    -se_exp (array-like): Standard errors of the exposure variable.
    -se_out (array-like): Standard errors of the outcome variable.
    -parameters (dict): Dictionary containing additional parameters, including "nboot" for the number of bootstrap iterations.
    Returns:
    result (dict): Dictionary containing the MR estimates, standard errors, and p-values for the causal effect.
    '''
  ```

  ```python
  mr_singlesnp(b_exp, b_out, se_exp, se_out, parameters)
    '''
    Description: This function performs Mendelian Randomization (MR) analysis for each exposure-outcome pair.
    Parameters:
    - b_exp (array-like): Beta values of the exposure variable.
    - b_out (array-like): Beta values of the outcome variable.
    - se_exp (array-like): Standard errors of the exposure variable.
    - se_out (array-like): Standard errors of the outcome variable.
    - parameters (dict): Dictionary containing additional parameters for MR analysis.
    Returns:
    result (DataFrame): A pandas DataFrame containing the MR estimates, standard errors, and p-values for each SNP within each
    exposure-outcome pair.
    '''
  ```

  ```python
  mr_two_sample_ml(b_exp, b_out, se_exp, se_out, parameters)
    '''
    Description: This function calculates Mendelian Randomization (MR) estimates using the two-sample instrumental
    variable method with maximum likelihood estimation.   
    Parameters:
    - b_exp (array-like): Beta values of the exposure variable.
    - b_out (array-like): Beta values of the outcome variable.
    - se_exp (array-like): Standard errors of the exposure variable.
    - se_out (array-like): Standard errors of the outcome variable.
    - parameters (dict): Dictionary containing additional parameters for MR analysis.
    Returns:
    result (dict): Dictionary containing the MR estimates, standard errors, p-values, and other statistics.
    '''
  ```

  ```python
  mr_uwr(b_exp, b_out, se_exp, se_out, parameters)
    '''
    Description: This function calculates Mendelian Randomization (MR) estimates using the unweighted regression method.
    Parameters:
    - b_exp (array-like): Beta values of the exposure variable.
    - b_out (array-like): Beta values of the outcome variable.
    - se_exp (array-like): Standard errors of the exposure variable.
    - se_out (array-like): Standard errors of the outcome variable.
    - parameters (dict): Dictionary containing additional parameters for MR analysis.
    Returns:
    result (dict): Dictionary containing the MR estimates, standard errors, p-values, and other statistics.
    '''
  ```

  ```python
  mr_wald_ratio(b_exp, b_out, se_exp, se_out, parameters)
    '''
    Description: This function calculates Mendelian Randomization (MR) estimates using the Wald ratio method.
    
    Parameters:
    - b_exp (array-like or scalar): Beta values of the exposure variable.
    - b_out (array-like or scalar): Beta values of the outcome variable.
    - se_exp (array-like or scalar): Standard errors of the exposure variable.
    - se_out (array-like or scalar): Standard errors of the outcome variable.
    - parameters (dict): Dictionary containing additional parameters for MR analysis.

    Returns:
    result (dict): Dictionary containing the MR estimates, standard errors, p-values, and the number of SNPs used for the analysis.
    '''
  ```

  ```python
  mr_weighted_median(b_exp, b_out, se_exp, se_out, parameters)
    '''
    Description: This function calculates Mendelian Randomization (MR) estimates using the weighted median method.
    
    Parameters:
    - b_exp (array-like): Beta values of the exposure variable.
    - b_out (array-like): Beta values of the outcome variable.
    - se_exp (array-like): Standard errors of the exposure variable.
    - se_out (array-like): Standard errors of the outcome variable.
    - parameters (dict): Dictionary containing additional parameters for MR analysis.

    Returns:
    result (dict): Dictionary containing the MR estimates, standard errors, p-values, Q-statistic, degrees of freedom (Q_df),
                   and p-value of the Q-statistic (Q_pval), as well as the number of SNPs used for the analysis (nsnp).
    '''
  ```

  ```python
  read_exposure_data(filename, clump=False, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf",
  effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", units_col="units", ncase_col="ncase", ncontrol_col="ncontrol",
  samplesize_col="samplesize", gene_col="gene", id_col="id", min_pval=1e-200, log_pval=False, chr_col="chr", pos_col="pos")
    '''
    Description: This function reads exposure data from a text file and processes it into a formatted DataFrame. The function allows for
    customization of column names and formatting options to ensure the data is appropriately structured. If required, it can also perform
    clumping of the data, which groups variants that are in linkage disequilibrium (LD) with each other based on specified parameters.

    Parameters:
      -filename (str): The path to the text file containing the exposure data.
      -clump (bool, optional): A flag to enable clumping of the data. Default is False.
      -sep (str, optional): The delimiter used in the text file to separate columns. Default is a space (" ").
      -phenotype_col (str, optional): The column name for the phenotype data in the text file. Default is "Phenotype".
      -snp_col (str, optional): The column name for SNP identifiers in the text file. Default is "SNP".
      -beta_col (str, optional): The column name for beta values (effect sizes) of the exposure variable in the text file. Default is "beta".
      -se_col (str, optional): The column name for standard errors of the exposure variable in the text file. Default is "se".
      -eaf_col (str, optional): The column name for the effect allele frequency (EAF) in the text file. Default is "eaf".
      -effect_allele_col (str, optional): The column name for the effect allele in the text file. Default is "effect_allele".
      -other_allele_col (str, optional): The column name for the other allele in the text file. Default is "other_allele".
      -pval_col (str, optional): The column name for p-values in the text file. Default is "pval".
      -units_col (str, optional): The column name for units of measurement in the text file. Default is "units".
      -ncase_col (str, optional): The column name for the number of cases in the text file. Default is "ncase".
      -ncontrol_col (str, optional): The column name for the number of controls in the text file. Default is "ncontrol".
      -samplesize_col (str, optional): The column name for the total sample size in the text file. Default is "samplesize".
      -gene_col (str, optional): The column name for gene information in the text file. Default is "gene".
      -id_col (str, optional): The column name for unique identifiers in the text file. Default is "id".
      -min_pval (float, optional): The minimum p-value threshold for the data. Default is 1e-200.
      -log_pval (bool, optional): A flag indicating if p-values in the text file are given in log-scale. Default is False.
      -chr_col (str, optional): The column name for chromosome information in the text file. Default is "chr".
      -pos_col (str, optional): The column name for base pair position information in the text file. Default is "pos".
    Returns:
      A DataFrame containing the formatted exposure data, with an additional column "data_source.exposure" indicating the data source.
    '''
  ```

  ```python
  read_outcome_data(filename, snps=None, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf",
  effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", units_col="units", ncase_col="ncase", ncontrol_col="ncontrol",
  samplesize_col="samplesize", gene_col="gene", id_col="id", min_pval=1e-  200, log_pval=False, chr_col="chr", pos_col="pos")
    '''
    Description: This function reads outcome data from a text file and formats it for Mendelian Randomization (MR) analysis.

    Parameters:
    - filename (str): The path to the text file containing the outcome data.
    - snps (list or None): A list of SNP IDs to include in the analysis. If None, all SNPs in the file will be used.
    - sep (str): The delimiter used in the text file.
    - phenotype_col (str): The column name for the outcome phenotype.
    - snp_col (str): The column name for SNP IDs.
    - beta_col (str): The column name for beta values.
    - se_col (str): The column name for standard errors of beta values.
    - eaf_col (str): The column name for effect allele frequencies.
    - effect_allele_col (str): The column name for the effect allele.
    - other_allele_col (str): The column name for the other allele.
    - pval_col (str): The column name for p-values.
    - units_col (str): The column name for units of measurement.
    - ncase_col (str): The column name for the number of cases (for binary outcomes).
    - ncontrol_col (str): The column name for the number of controls (for binary outcomes).
    - samplesize_col (str): The column name for the total sample size.
    - gene_col (str): The column name for gene names or identifiers (optional).
    - id_col (str): The column name for individual or study identifiers (optional).
    - min_pval (float): The minimum p-value to use for clumping (default is 1e-200).
    - log_pval (bool): If True, assumes p-values are in log scale (default is False).
    - chr_col (str): The column name for chromosome information (optional).
    - pos_col (str): The column name for base pair positions (optional).

    Returns:
    outcome_dat (pd.DataFrame): A formatted DataFrame containing the outcome data for MR analysis.
    '''
  ```

  ```python
  weighted_median(b_iv, weights)
  '''
  Description: The weighted_median function calculates the weighted median of a given set of beta values (b_iv) using provided weights.
  It is used in Mendelian Randomization (MR) analysis to estimate causal effects by combining instrumental variable (IV) effect estimates.

  Parameters:
  -b_iv (array-like): Array of instrumental variable effect estimates (betaIV) for each SNP.
  -weights (array-like): Array of weights corresponding to each betaIV estimate.

  Returns:
  b (float): Weighted median estimate of the instrumental variable effect (b_iv) using the given weights.
  '''
  ```
  ```python
  weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, weights, nboot)
  '''
  Description: The weighted_median_bootstrap function calculates the standard error of the weighted median estimate of instrumental
  variable (IV) effect in Mendelian Randomization (MR) analysis using bootstrap resampling.

  Parameters:  
  -b_exp (array-like): Beta values of the exposure variable.
  -b_out (array-like): Beta values of the outcome variable.
  -se_exp (array-like): Standard errors of the exposure variable.
  -se_out (array-like): Standard errors of the outcome variable.
  -weights (array-like): Array of weights corresponding to each IV effect estimate.
  -nboot (int): Number of bootstrap iterations.

  Returns:
  se (float): Standard error of the weighted median estimate calculated using bootstrap resampling.
    '''
  ```
