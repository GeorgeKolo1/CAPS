import numpy as np
import pytest
import os
from phenotype_association import AssociationTest, SeperationInvestigator

@pytest.fixture
def perfect_separation():
    subtypes = np.array(['A', 'B', 'A', 'B', 'A', 'B'])
    phenotype = np.array([1, 2, 1, 2, 1, 2])
    return subtypes, phenotype

@pytest.fixture
def no_separation():
    rng = np.random.default_rng(42)
    subtypes = rng.choice(['A', 'B', 'C'], size=100)
    phenotype = rng.choice(['X', 'Y'], size=100)
    seperator = False
    return subtypes, phenotype, seperator

def test_separation_investigator_detects_perfect_separation(perfect_separation, no_separation):
    
    arr, phenotype = perfect_separation
    arr2, phenotype2 = no_separation

    assert SeperationInvestigator(arr, phenotype) is True
    assert SeperationInvestigator(arr2, phenotype2) is False


def test_firth_computation(perfect_separation):
    
    arr, phenotype, seperator = perfect_separation
    df = AssociationTest(arr, phenotype, seperator)


    assert (df['odds_ratio'] >= df['CI_low']).all()
    assert (df['odds_ratio'] <= df ['CI_high']).all()
    assert (df['test'] == "FirthLR")


def test_OR_computation(no_separation):
    
    arr, phenotype, seperator = no_separation
    df = AssociationTest(arr, phenotype, seperator)

    assert (df['odds_ratio'] >= df['CI_low']).all()
    assert (df['odds_ratio'] <= df ['CI_high']).all()
    assert (df['test'] == "FisherExact")
    