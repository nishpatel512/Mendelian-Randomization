import ieugwaspy as igd
import pandas as pd
def available_outcomes():
    data = igd.gwasinfo()
    df = pd.DataFrame.from_dict(data, orient='index')
    return df


print (available_outcomes())