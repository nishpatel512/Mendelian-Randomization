# Mandelian-Randomization
# Introduction
Mendelian randomization is a method of using measured variation in genes of known function to examine the causal effect of a modifiable exposure on disease in observational studies.
We are making use of python package called MRPackage which has all the functionalities to perfom MR analysis using data from local computer or data from GWAS.

# Methods:

  ```python
  available_outcomes()
  '''
    Description: Retrieves GWAS (Genome-Wide Association Study) information using the ieugwaspy library and returns it as a Pandas DataFrame.
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
  format_data()
  '''
    Description: This function is used to format and preprocess genetic association data for Mendelian randomization (MR) analysis. It takes a pandas DataFrame containing genetic     
    association data, which typically includes information about the exposure (e.g., risk factor) and outcome (e.g., disease) variables. The function performs several data cleaning and 
    manipulation steps to ensure the data is in a suitable format for conducting MR analysis.

    Parameters:
    - dat: pandas DataFrame
    Input data containing genetic association information, such as beta coefficients, standard errors, p-values, allele frequencies, sample sizes, etc.

    - type: str (optional, default: "exposure")
    Specifies whether the data represents the "exposure" or "outcome" variable.

    - snps: list (optional, default: None)
    A list of specific SNPs to be included in the analysis. If provided, the function filters the data to include only the specified SNPs.

    - header: bool (optional, default: True)
    Specifies whether the input DataFrame has a header row.

    - phenotype_col, snp_col, beta_col, se_col, eaf_col, effect_allele_col, other_allele_col, pval_col, units_col, ncase_col, ncontrol_col, samplesize_col, gene_col, id_col, min_pval, 
    z_col, info_col, chr_col, pos_col: str (optional)
    Column names corresponding to specific genetic association data in the input DataFrame.

    - log_pval: bool (optional, default: False)
    Specifies whether p-values are provided in logarithmic scale (log10).

    Returns:
    - dat: pandas DataFrame
    The formatted and preprocessed data, ready for Mendelian randomization analysis.
  '''
  ```

  ```python
  harmonise_data()
  '''
    Description: In order to perform MR the effect of a SNP on an outcome and exposure must be harmonised to be relative to the same allele.
    Parameters: Takes formated input data of Exposure and Outcome data
    Returns: Data frame with harmonised effects and alleles
  '''
  ```

  ```python
  mr()
  '''
    Description: The harmonise_data function is used to perform harmonization of genetic effect sizes and alleles for Mendelian Randomization (MR) analysis. It ensures that the effect of      a Single Nucleotide Polymorphism (SNP) on both the exposure and outcome variables is measured relative to the same allele.
    Parameters:
    exposure_data (pandas.DataFrame): DataFrame containing the exposure data with columns such as "SNP," "effect_size," "standard_error," "effect_allele," and "other_allele."
    outcome_data (pandas.DataFrame): DataFrame containing the outcome data with columns such as "SNP," "effect_size," "standard_error," "effect_allele," and "other_allele."
    Returns:
    harmonized_data (pandas.DataFrame): Data frame with harmonized genetic effect sizes and alleles, where the effect sizes are relative to the same allele for both exposure and outcome       data.
  '''
  ```

  ```python
  mr_egger_regression()
  '''
    Description: The harmonise_data function is used to perform harmonization of genetic effect sizes and alleles for Mendelian Randomization (MR) analysis. It ensures that the effect of      a Single Nucleotide Polymorphism (SNP) on both the exposure and outcome variables is measured relative to the same allele.
    Parameters:
    exposure_data (pandas.DataFrame): DataFrame containing the exposure data with columns such as "SNP," "effect_size," "standard_error," "effect_allele," and "other_allele."
    outcome_data (pandas.DataFrame): DataFrame containing the outcome data with columns such as "SNP," "effect_size," "standard_error," "effect_allele," and "other_allele."
    Returns:
    harmonized_data (pandas.DataFrame): Data frame with harmonized genetic effect sizes and alleles, where the effect sizes are relative to the same allele for both exposure and outcome data.
  '''
  ```

  ```python
  mr_egger_regression_bootstrap()
  '''
    Create a blank plot with a centered text message.
    Parameters:
        message (str): The message to be displayed in the plot.
    Returns:
        matplotlib.axes._subplots.AxesSubplot: The blank plot with the message.
  '''
  ```
