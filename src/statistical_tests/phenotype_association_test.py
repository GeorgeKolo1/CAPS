import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from src.statistical_tests.wallace_coefficient import CT
from typing import Optional

def AssociationTest(arr1, phenotype, outfile: Optional[str] = None):
    '''
    Function to test for subtype (within a subtype method) associations with a phenotype
    
    The method deconstructs each subtype and each phenotype (categorical - nominal) into a 2x2 contingency table where the two rows are each subtype and "other subtypes" and each column is each phenotype and "other phenotypes"
    This is done iteratively so that an effect size (odds ratio) can be compted for each subtype with each phenotype

    Args:
        arr1 (:pandas Series or numpy array: 'str'): A pandas Series or 1D numpy array that contains the different subtypes (as strings)
        phenotype (:pandas Series or numpy array: 'str'): A pandas Series of 1D numpy array that contains the different phenotypes (as strings)
        outfile (:str): A string that contains the directory and filename (prefix) in which to save the results

    Returns:
        
    '''
    if isinstance(arr1, np.ndarray) == False:
        arr1 = arr1.to_numpy()

    if isinstance(phenotype, np.ndarray) == False:
        phenotype = phenotype.to_numpy()

    for i in arr1:
        tmp_arr = np.where(arr1 == i, arr1, 'other')

        for y in phenotype:
            tmp_phenotype = np.where(phenotype == y, phenotype, 'other_phenotype')


        ct = pd.crosstab(tmp_arr, tmp_phenotype)
        res= fisher_exact(ct)
        OR = odds_ratio(ct)

        if outfile is not None:
            with open(outfile, 'w') as f:
                f.write(f'subtype {i}, phenotype {y}')
                f.write(f'statistic {res.statistic}, pvalue {res.pvalue}, odds ratio {OR}, oddsratio 95% CI {OR.confidence_interval}')

        return res.statistic, res.pvalue, OR, OR.confidence_interval


        