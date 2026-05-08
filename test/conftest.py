import pandas as pd
import pytest
from pathlib import Path

TEST_DATA_DIR = Path(__file__).parent / "test_data_folder"

@pytest.fixture
def test_data():
    return pd.read_csv(TEST_DATA_DIR / "test_data.csv")

@pytest.fixture
def subtype_arrays(test_data):
    arr = test_data.to_numpy()

    return arr[:, 0], arr[:, 1]