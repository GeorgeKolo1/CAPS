import numpy as np
import pandas as pd
from src.statistical_tests import AdjustedWallace

def test_adjusted_wallace():
    '''Function that tests the adjusted wallace function
    
    This function tests the adjusted wallace function by using test data provided in the test_data_folder directory. It runs the adjusted wallace coefficient on two subtyping method classifications
    and asserts whether the adjusted wallace coefficient is between 0 and 1, and that it is equal to the known correct value of the output of the adjusted wallace coefficient when run using other means

    Args:

    Returns:
        A pytest assert statement to ensure that the test data is running correctly and the adjusted wallace coefficient is behaving as expected    
    '''
    
    data = pd.read_csv('test/test_data_folder/test_data.csv')
    data = data.to_numpy()

    arr1 = data[:, 0]
    arr2 = data[:, 1]

    assert len(arr1) == 325
    assert len(arr2) == 325

    AW_ab, AW_ba = AdjustedWallace(arr1, arr2)

    assert AW_ab != AW_ba
    assert AW_ab > 0 and AW_ab < 1
    assert AW_ba > 0 and AW_ba < 1
    assert round(AW_ba, 3) == 0.807
    assert round(AW_ab, 3) == 0.608

