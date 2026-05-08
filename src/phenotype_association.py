import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from statistical_tests import CT
from firthmodels import FirthLogisticRegression
from firthmodels import detect_separation
from typing import Optional
import os


def AssociationTest(arr1: pd.Series | np.ndarray, phenotype: pd.Series | np.ndarray, seperation: bool) -> pd.DataFrame:
    '''
    Function to test for subtype (within a subtype method) associations with a phenotype

    The method deconstructs each subtype and each phenotype (categorical - nominal) into a 2x2 contingency table where the two rows are each subtype and "other subtypes" and each column is each phenotype and "other phenotypes"
    This is done iteratively so that an effect size (odds ratio) can be computed for each subtype with each phenotype

    Args:
        arr1 (:pandas Series or numpy array: 'str'): A pandas Series or 1D numpy array that contains the different subtypes (as strings)
        phenotype (:pandas Series or numpy array: 'str'): A pandas Series of 1D numpy array that contains the different phenotypes (as strings)
        outfile (:str): A string that contains the directory and filename (prefix) in which to save the results; Defaults to current working directory when not specified or left as None (default)

    Returns:

    Raises:


    '''
    if isinstance(arr1, np.ndarray) == False:
        arr1 = arr1.to_numpy()

    if isinstance(phenotype, np.ndarray) == False:
        phenotype = phenotype.to_numpy()

    results = []

    for i in np.unique(arr1):
        tmp_arr = np.where(arr1 == i, 1, 0)

        for y in np.unique(phenotype):
            tmp_phenotype = np.where(phenotype == y, 1, 0)

        if seperation == True:
            model = FirthLogisticRegression().fit(tmp_arr.reshape(-1, 1), tmp_phenotype.reshape(-1, 1))
            CI = model.conf_int()
            results.append({
                "subtype": i,
                "phenotype": y,
                "test": "FirthLR",
                "coefficient": model.coef_[0],
                "odds_ratio": np.exp(model.coef_[0]),
                "pvalue": model.pvalues_[0],
                "CI_low": np.exp(CI[1][0]),
                "CI_high": np.exp(CI[1][1])
            })

        else:
            ct = pd.crosstab(tmp_arr, tmp_phenotype)
            res_or = odds_ratio(ct)
            res = fisher_exact(ct)
            OR = res_or.statistic
            CI = res_or.confidence_interval()
            results.append({
                "subtype": i,
                "phenotype": y,
                "test": "FisherExact",
                "statistic": res.statistic,
                "odds_ratio": OR,
                "pvalue": res.pvalue,
                "CI_low": CI.low,
                "CI_high": CI.high
            })

    results_df = pd.DataFrame(results)

    return results_df

def SeperationInvestigator(arr1: pd.Series | np.ndarray, phenotype: pd.Series | np.ndarray) -> bool:
    '''
    Function to investigate whether there is seperation in the data that would require the use of Firth logistic regression to compute odds ratios for the association analysis

    Args:
        arr1 (:pandas Series or numpy array: 'str'): A pandas Series or 1D numpy array that contains the different subtypes (as strings)
        phenotype (:pandas Series or numpy array: 'str'): A pandas Series of 1D numpy array that contains the different phenotypes (as strings)

    Returns:
        bool: True if seperation is detected, False if not
    Raises:

    '''
    if isinstance(arr1, np.ndarray) == False:
        arr1 = arr1.to_numpy()

    if isinstance(phenotype, np.ndarray) == False:
        phenotype = phenotype.to_numpy()

    seperation_detected = False

    #list to store boolean values for whether seperation is detected for each subtype-phenotype combination; if any are True then seperation is detected and Firth logistic regression should be used to compute odds ratios
    bool_list = []

    for i in np.unique(arr1):
        tmp_arr = np.where(arr1 == i, 1, 0)

        for y in np.unique(phenotype):
            tmp_phenotype = np.where(phenotype == y, 1, 0)

            seperation = detect_separation(tmp_arr.reshape(-1, 1), tmp_phenotype.reshape(-1, 1))
            bool_list.append(seperation.separation)
        
    if any(bool_list) == True:
        seperation_detected = True
    else:
        seperation_detected = False

    return seperation_detected

