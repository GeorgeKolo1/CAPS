from src.statistical_tests import DiscriminationIndex
import pandas as pd

def test_DiscriminationIndex():
    data = pd.read_csv('test/test_data_folder/test_data.csv')
    sm = data['T type']
    D, CI_low, CI_high = DiscriminationIndex(sm)

    assert D >= CI_low
    assert D <= CI_high
    assert D == 0.721
    assert CI_low == 0.676
    assert CI_high == 0.767