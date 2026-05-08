import numpy as np
import pytest
import os
from src.phenotype_association import AssociationTest, SeperationInvestigator

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
    return subtypes, phenotype


def test_association_test(tmp_path):
 
    

def test_separation_investigator_detects_perfect_separation(perfect_separation):
    arr, phenotype = perfect_separation
    assert SeperationInvestigator(arr, phenotype) is True
# ...etc

    outfile = os.path.join(tmp_path, 'test_output')

    df = AssociationTest(arr1, phenotype, outfile)
    df_sep = AssociationTest(arr2, phenotype_sep, outfile)

    assert os.path.exists(outfile)
    assert (df['odds_ratio'] >= df['CI_low']).all()
    assert (df['odds_ratio'] <= df ['CI_high']).all()
    assert (df_sep['test'] == "FirthLR")
    