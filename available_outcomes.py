import ieugwaspy as igd
import pandas as pd
def available_outcomes():
    '''
    Retrieves GWAS (Genome-Wide Association Study) information using the ieugwaspy library
    and returns it as a Pandas DataFrame.
    Returns:
    df (pandas.DataFrame): DataFrame containing GWAS information, with outcomes as the index.
    '''
    data = igd.gwasinfo()
    df = pd.DataFrame.from_dict(data, orient='index')
    return df


print (available_outcomes())
