from src.statistical_tests import DiscriminationIndex
from src.statistical_tests import AdjustedWallace
from src.statistical_tests import Wallace

def test_discrimination_index_95_ci(test_data):
    D, CI_low, CI_high = DiscriminationIndex(test_data['T type'])

    assert D >= CI_low
    assert D <= CI_high

def test_discrimination_index_value(test_data):
    D, CI_low, CI_high = DiscriminationIndex(test_data['T type'])

    assert D == 0.721
    assert CI_low == 0.676
    assert CI_high == 0.767

def test_adjusted_wallace_function(subtype_arrays):

    arr1, arr2 = subtype_arrays

    AW_ab, AW_ba = AdjustedWallace(arr1, arr2)

    assert AW_ab != AW_ba
    assert AW_ab > 0 and AW_ab < 1
    assert AW_ba > 0 and AW_ba < 1


def test_adjusted_wallace_value(subtype_arrays):

    arr1, arr2 = subtype_arrays
    
    AW_ab, AW_ba = AdjustedWallace(arr1, arr2)

    assert round(AW_ba, 3) == 0.807
    assert round(AW_ab, 3) == 0.608


def test_wallace_function(subtype_arrays):

    arr1, arr2 = subtype_arrays

    wc_ab, wc_ba = Wallace(arr1, arr2)

    assert wc_ab != wc_ba
    assert wc_ab > 0 and wc_ab <= 1
    assert wc_ba > 0 and wc_ba <= 1


def test_wallace_value(subtype_arrays):

    arr1, arr2 = subtype_arrays

    wc_ab, wc_ba = Wallace(arr1, arr2)

    assert round(wc_ab, 3) == 0.696
    assert round(wc_ba, 3) == 0.861