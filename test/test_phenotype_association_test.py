import numpy as np
import pandas as pd
import os
from src.statistical_tests.phenotype_association_test import AssociationTest

def test_association_test(tmp_path):
    '''
    
    
    
    
    
    
    '''
    data = pd.read_csv('test/test_data_folder/test_data.csv')
    data = data.to_numpy()
    arr1 = data[:, 1]
    
    np.random.seed(42)
    phenotype = np.random.choice(['A', 'B', 'C', 'D', 'E'], size=len(arr1))
    assert  len(arr1) == len(phenotype)

    outfile = os.path.join(tmp_path, 'test_output.txt')
    statistic, pvalue, OR, CI = AssociationTest(arr1, phenotype, outfile)

    assert os.path.exists(outfile)
    assert statistic != 0
    assert pvalue != 0

    