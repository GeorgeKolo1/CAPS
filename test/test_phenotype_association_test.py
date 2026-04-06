import numpy as np
import pandas as pd
import os
from src.statistical_tests.phenotype_association import AssociationTest

def test_association_test(tmp_path):
    '''
    
    
    
    
    
    
    '''
    #Set up test data
    data = pd.read_csv('test/test_data_folder/test_data.csv')
    data = data.to_numpy()
    arr1 = data[:, 1]

    #arr2 to test perfect seperation causes Firth LR to be used
    arr2 = np.array(['A', 'B', 'A', 'B', 'A', 'B']) 
    np.random.seed(42)
    phenotype = np.random.choice(['A', 'B', 'C', 'D', 'E'], size=len(arr1))

    #Phenotype for perfect association testing
    phenotype_sep = np.array([1, 2, 1, 2, 1, 2])
    assert  len(arr1) == len(phenotype)

    outfile = os.path.join(tmp_path, 'test_output')

    df = AssociationTest(arr1, phenotype, outfile)
    df_sep = AssociationTest(arr2, phenotype_sep, outfile)

    assert os.path.exists(outfile)
    assert (df['odds_ratio'] >= df['CI_low']).all()
    assert (df['odds_ratio'] <= df ['CI_high']).all()
    assert (df_sep['test'] == "FirthLR")
    