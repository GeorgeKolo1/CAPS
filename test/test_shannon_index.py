import numpy as np
import pandas as pd
from src.statistical_tests.shannon_index import ShannonIndex


def test_shannon_index():
    '''
    Function to test the Shannon Index function

    This function uses test data that contains differnet subtypes of a subtyping method and tests them
    using the Shannon Index function.

    Args:
        None
    Returns:
        None
    Raises:
        Assert
    '''
    df = pd.read_csv('test/test_data_folder/test_data.csv')
    df = df.to_numpy()
    sm = df[:, 1]

    unique_subtypes, H, H_max, H_min = ShannonIndex(sm)

    assert unique_subtypes == 12
    assert H != 0
    assert H_max >= H
    assert H_min <= H
    assert H == 2.603
    assert H_max == 3.585
    assert H_min == 0.330

