import pandas as pd
import numpy as np
from src import wallace_coefficient as wc

def test_wallace():
    """Function to test the wallace coefficient function
    
    This function takes test data and runs the wallace coefficient function on it, then checks if 
        out is as expected 
        
    Args:
        
    Returns:
        A pytest assert statement that checks if the output of the wallace coefficient function is as expected
"""
    data = pd.read_csv('test/test_data_folder/test_data.csv')
    data = data.to_numpy()

    assert len(data[:, 0]) == 325
    assert len(data[:, 1]) == 325


    wc_ab, wc_ba = wc.Wallace(data[:, 0], data[:, 1])

    assert wc_ab != wc_ba
    assert wc_ab > 0 and wc_ab <= 1
    assert wc_ba > 0 and wc_ba <= 1
    assert round(wc_ab, 3) == 0.696
    assert round(wc_ba, 3) == 0.861