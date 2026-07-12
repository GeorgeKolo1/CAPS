import numpy as np
import math
import pandas as pd

def DiscriminationIndex(sm):
    if isinstance(sm, pd.DataFrame) == True:
        return("This is a dataframe, series or np array required")
    elif isinstance(sm, pd.Series) == True:
        sm.to_numpy
    
    #discrimination index calculation (Hunter and Gaston)

    N = len(sm)
    s  = np.unique(sm, return_counts=True)
    summation_list =[]
    for i in s[1]:
        y = (i * (i-1))
        summation_list.append(y)
    D = 1 - (sum(summation_list)) / (N * (N-1))
    #95% CI calculation (Grundmann, Hori & Tanner, 2001)

    summation_pi_j3 = []
    summation_pi_j2 = []
    for i in s[1]:
        pi_j = i/N
        pi_j3 = pi_j ** 3
        pi_j2 = pi_j ** 2
        summation_pi_j3.append(pi_j3)
        summation_pi_j2.append(pi_j2)

    pi_j3 = sum(summation_pi_j3)
    pi_j2 = (sum(summation_pi_j2) ** 2)
    S2 = (4/N) * (pi_j3 - pi_j2)
    CI_low = (D - (2 * math.sqrt(S2)))
    CI_high = (D + (2 * math.sqrt(S2)))

    return round(D, 3), round(CI_low, 3), round(CI_high, 3)

#------------------------------------------------------------------------------------------------------------------------------------------

# Compute the contingency table for two arrays (subtyping methods)
def CT(arr1, arr2):
    
    ct_ab = pd.crosstab(arr1, arr2)
    ct_ab = ct_ab.to_numpy()
    ct_ba = pd.crosstab(arr2, arr1)
    ct_ba = ct_ba.to_numpy()

    return ct_ab, ct_ba

def Wallace(arr1, arr2):
    ''' Compute the wallace coefficient for two arrays (subtyping methods) in both directions'''

    ct_ab, ct_ba = CT(arr1, arr2)

    sum_ct_ab = np.sum((ct_ab * (ct_ab - 1)))
    row_sum_ab = np.sum(ct_ab, axis=1)
    denominator_ab = np.sum(row_sum_ab * (row_sum_ab - 1))

    sum_ct_ba = np.sum((ct_ba * (ct_ba - 1)))
    row_sum_ba = np.sum(ct_ba, axis=1)
    denominator_ba = np.sum(row_sum_ba * (row_sum_ba - 1))

    wallace_ab = sum_ct_ab / denominator_ab
    wallace_ba = sum_ct_ba / denominator_ba


    return wallace_ab, wallace_ba
                
def AdjustedWallace(arr1, arr2):
    '''Calculates the adjusted wallace coefficient 
    
    This function computes the adjusted Wallace coefficient which is = W(a->b) - Wi(a->b) / (1 - Wi(a->b))
    Where W(a->b) is the wallace coefficient in direction a->b and Wi(a->b) is the wallace coefficient under independence,
    in other words that wallace coefficient if the two subtyping methods were completely independent.
        
    Args:
        arr1 - pandas series or numpy array containing the subtypes of subtyping method a
        arr2 - pandas or numpy array containing the subtypes of subtyping method b
        wallace_ab - wallace coefficient in diretion a->b
        wallace_ba - wallace coefficient in diretion b->a 
        
    Returns:
        AW_ab (float) - adjusted wallace coefficient in direction a->b
        AW_ba (float) - adjusted wallace coefficient in direction b->a

    Raises:
        ValueError: 1.If the input arrays do not have the same number of samples
                    2. Are not one-dimensional
        
    '''
    if arr1.shape[0] != arr2.shape[0]:
        raise ValueError("Input arrays must have the same number of samples")
    if len(arr1.shape) != 1 or len(arr2.shape) != 1:
        raise ValueError("Input arrays must be one-dimensional")
    if arr1.dtype != arr2.dtype:
        arr1 = arr1.astype(str)
        arr2 = arr2.astype(str)
    
    Da, Da_low, Da_high = DiscriminationIndex(arr1)
    Db, Db_low, Db_high = DiscriminationIndex(arr2)

    wallace_ab, wallace_ba = Wallace(arr1, arr2)

    Wi_ab = 1 - Db
    Wi_ba = 1 - Da

    AW_ab = (wallace_ab - Wi_ab) / (1 - Wi_ab)
    AW_ba = (wallace_ba - Wi_ba) / (1- Wi_ba)

    return AW_ab, AW_ba