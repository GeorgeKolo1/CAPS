import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from src.statistical_tests.wallace_coefficient import CT
from firthmodels import FirthLogisticRegression
from firthmodels import detect_separation
from typing import Optional
import os


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
        
    Raises:

        
    '''
    if isinstance(arr1, np.ndarray) == False:
        arr1 = arr1.to_numpy()

    if isinstance(phenotype, np.ndarray) == False:
        phenotype = phenotype.to_numpy()


    for i in arr1:
        tmp_arr = np.where(arr1 == i, 1, 0)
        tmp_arr = tmp_arr.reshape(-1, 1)

        for y in phenotype:
            tmp_phenotype = np.where(phenotype == y, 1, 0)
            tmp_phenotype = tmp_phenotype.reshape(-1, 1)

        seperation = detect_separation(tmp_arr, tmp_phenotype)

        if seperation == True:
            print('Sepeartion Detected! Will use Firth Penalized Logistic Regression to compute effect size (odds ratio)')
            model = FirthLogisticRegression().fit(tmp_arr, tmp_phenotype)
            CI_low, CI_high = model.conf_int()

            if outfile is None:
                print("OUTFILE IS NONE! Will save to working directory...")
                tmp_result_dict = {"subtype": i, "coefficient" : model.coef_[0], "odds_ratio" : np.exp(model.coef_[0]), "pvalues" : model.pvalues_[0], "CI_low" : np.exp(CI_low), "CI_high" : np.exp(CI_high)}
                results_df = pd.DataFrame(data=tmp_result_dict)
                results_df.to_csv(outfile + f'_{y}_FirthLR.csv')

            else:
                tmp_result_dict = {"subtype": i, "coefficient" : model.coef_[0], "odds_ratio" : np.exp(model.coef_[0]), "pvalues" : model.pvalues_[0], "CI_low" : np.exp(CI_low), "CI_high" : np.exp(CI_high)}
                results_df = pd.DataFrame(data=tmp_result_dict)
                results_df.to_csv(outfile + f'_{y}_FirthLR.csv')
                

        else:
            ct = pd.crosstab(tmp_arr, tmp_phenotype)

            res_or = odds_ratio(ct)
            res = fisher_exact(ct)
            OR = res_or.statistic
            CI = res_or.confidence_interval()

            if outfile is None:
                print("OUTFILE IS NONE! Will save to working directory...")
                tmp_result_dict = {'subtype': i, "statistic" : res.statistic, 'pvalue' : res.pvalue, 'odds_ratio' : OR, "CI_low": CI.low, 'CI_high': CI.high}
                results_df = pd.DataFrame(data=tmp_result_dict)
                outfile = os.getcwd()
                results_df.to_csv(outfile + f'_{y}_OR.csv')

            else:
                tmp_result_dict = {'subtype': i, "statistic" : res.statistic, 'pvalue' : res.pvalue, 'odds_ratio' : OR, "CI_low": CI.low, 'CI_high': CI.high}
                results_df = pd.DataFrame(data=tmp_result_dict)
                outfile = os.getcwd()
                results_df.to_csv(outfile + f'_{y}_OR.csv')
        
    
    


        