import numpy as np
import pytest
import os
from phenotype_association import AssociationTest, SeperationInvestigator

@pytest.fixture
def perfect_separation():
    subtypes = np.array(['A', 'B', 'A', 'B', 'A', 'B'])
    phenotype = np.array([1, 2, 1, 2, 1, 2])
    seperator= True
    return subtypes, phenotype, seperator

@pytest.fixture
def quasi_complete_separation():
    """
    45 samples with quasi-complete separation in one cell.

    Subtype A is perfectly associated with phenotype X (triggers detection),
    while B and C have realistic mixed phenotypes so the other coefficients
    aren't degenerate. This is the regime Firth was designed for.
    """
    rng = np.random.default_rng(42)
    subtypes = np.array(
        ['A'] * 15
        + ['B'] * 15
        + ['C'] * 15
    )
    phenotype = np.concatenate([
        np.array(['X'] * 15),                              # A → always X (separated)
        rng.choice(['X', 'Y'], size=15, p=[0.3, 0.7]),     # B → mostly Y
        rng.choice(['X', 'Y'], size=15, p=[0.5, 0.5]),     # C → balanced
    ])
    return subtypes, phenotype, True

@pytest.fixture
def no_separation():
    rng = np.random.default_rng(42)
    subtypes = rng.choice(['A', 'B', 'C'], size=100)
    phenotype = rng.choice(['X', 'Y'], size=100)
    seperator = False
    return subtypes, phenotype, seperator

def test_separation_investigator_detects_perfect_separation(perfect_separation, no_separation):
    
    arr, phenotype, _ = perfect_separation
    arr2, phenotype2, _ = no_separation

    assert SeperationInvestigator(arr, phenotype) is True
    assert SeperationInvestigator(arr2, phenotype2) is False

def test_separation_investigator_detects_quasi_separation(quasi_complete_separation):
    arr, phenotype, _ = quasi_complete_separation
    assert SeperationInvestigator(arr, phenotype) is True

def test_firth_used_with_quasi_separation(quasi_complete_separation):
    arr, phenotype, sep = quasi_complete_separation
    df = AssociationTest(arr, phenotype, sep)
    assert (df['test'] == "FirthLR").all()

def test_firth_odds_ratio_within_confidence_interval(quasi_complete_separation):
    arr, phenotype, sep = quasi_complete_separation
    df = AssociationTest(arr, phenotype, sep)
    assert (df['CI_low'] <= df['odds_ratio']).all()
    assert (df['odds_ratio'] <= df['CI_high']).all()

def test_firth_computation(perfect_separation):
    
    arr, phenotype, seperator = perfect_separation
    df = AssociationTest(arr, phenotype, seperator)

    assert (df['test'] == "FirthLR").all()

def test_OR_computation(no_separation):
    
    arr, phenotype, seperator = no_separation
    df = AssociationTest(arr, phenotype, seperator)

    assert (df['odds_ratio'] >= df['CI_low']).all()
    assert (df['odds_ratio'] <= df ['CI_high']).all()
    assert (df['test'] == "FisherExact").all()
    