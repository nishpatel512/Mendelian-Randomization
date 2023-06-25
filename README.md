# Mandelian-Randomization
# Introduction
Mendelian randomization is a method of using measured variation in genes of known function to examine the causal effect of a modifiable exposure on disease in observational studies.
We are making use of python package called MRPackage which has all the functionalities to perfom MR analysis using data from local computer or data from GWAS.

# Methods:

  ```python
  available_outcomes()
  '''
    OUTPUT: a list of available instruments from the IEU-GWAS public database.
  '''
  ```
  ```python
  extract_instruments(list_ids)
  '''
    INPUT: List of instrument ids from the IEU_GWAS public database. Subset of available_outcomes() method output. 
    OUTPUT: Formatted dataframe of extracted instruments with ids from list_ids, to find exposure data.
  '''
  ```
  ```python
  extract_outcome_data()
  '''
    INPUT: 2 formatted dataframes from extract_instruments(from IEU_GWAS public database) and read_outcome_data(from local directory).
    OUTPUT: Formatted dataframe which contains outcome data
  '''
  ```
   ```python
  format_data()
  '''
    INPUT: 2 formatted dataframes from extract_instruments(from IEU_GWAS public database) and read_outcome_data(from local directory).
    OUTPUT: Formatted dataframe which contains outcome data
  '''
  ```
