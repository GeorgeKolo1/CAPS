import numpy as np
import pandas as pd

def ShannonIndex(sm):
    '''
    Function to compute the Shannon Index for a given subtype method

    Shannon Index is a measure of diversity that takes into account 
    the number of subtypes and the evenness of their distribution. 
    A higher Shannon Index indicates greater diversity, while a lower index indicates less diversity.
    It is calculated using the formula: H = -∑(p_i * log(p_i)), 
    where p_i is the proportion of each subtype in the dataset.

    Args:
        sm (:pandas Series or numpy array: 'str'): A pandas Series or 1D numpy array that contains the different subtypes (as strings)

    Returns:
        H (float): The Shannon Index for the given subtyping method
        H_max (float): The maximum possible Shannon Index for the given subtyping method
        H_min (float): The minimum possible Shannon Index for the given subtyping method

    Raises:
        ValueError: If the input is not a panadas Series or a numpy array
    '''
    if isinstance(sm, pd.Series) == False and isinstance(sm, np.ndarray) == False:
        raise ValueError('Input must be a pandas Series or numpy array')

    sm = pd.Series(sm).astype(str)

    # Compute the proportion of each subtype
    unique, counts = np.unique(sm, return_counts=True)

    if sum(counts) != len(sm):
        raise ValueError('The input contains values that are not strings. The sum of counts does not equal the number of individuals')

    proportions = counts / len(sm)

    #Take the negative sum of the proportions multiplied by the natural log of the proportions
    H = -np.sum(proportions * np.log(proportions))

    S = len(unique)
    Q = len(sm)

    #Compute the maximum possible Shannon Index using natural log of Q (unique subtypes)
    H_max = np.log(S)

    #Compute the minimum possible Shannon Index as (ln(Q) - (Q - S + 1)) * (ln(Q - S + 1) / Q) where S is the number of unique subtypes and Q is the number of individuals
    H_min = ((np.log(Q) - (Q - S + 1)) * (np.log((Q - S + 1)) / Q))

    return len(unique), round(H, 3), round(H_max, 3), round(H_min, 3)