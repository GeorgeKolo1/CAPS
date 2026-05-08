import numpy as np
import pandas as pd
from src.input_output_processor import InputProcessor as ip

def test_input_processor():
    '''
    Function to test the input processor function
    
    This function tests the input processor function by using test data provided in the test_data_folder directory. It runs the input processor on the test data and asserts whether the output is a pandas dataframe, and that the dataframe has the expected columns and number of rows.

    Args:

    Returns:
        A pytest assert statement to ensure that the test data is running correctly and the input processor is behaving as expected    
    '''
    df = ip('test/test_data_folder/test_data.csv')

    assert isinstance(df, pd.DataFrame)
    assert len(df) == 325
    assert 'T type' in df.columns