import numpy as np
import pandas as pd
import pytest
from input_output_processor import DataProcessor, OutputProcessor

@pytest.fixture
def sample_dataframe():
    """Small in-memory DataFrame with a phenotype column and two subtype columns."""
    return pd.DataFrame({
        "S type": ["A", "B", "A", "C", "B"],
        "T type": ["X", "Y", "X", "Y", "X"],
        "phenotype": ["R", "S", "R", "S", "R"],
    })


@pytest.fixture
def csv_file(tmp_path, sample_dataframe):
    """Path to a CSV file written from sample_dataframe."""
    path = tmp_path / "input.csv"
    sample_dataframe.to_csv(path, index=False)
    return str(path)

@pytest.fixture
def tsv_file(tmp_path, sample_dataframe):
    """Path to a TSV file written from sample_dataframe."""
    path = tmp_path / "input.tsv"
    sample_dataframe.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def unsupported_extension_file(tmp_path):
    """Path to a file with an unsupported extension (for the TypeError branch)."""
    path = tmp_path / "input.xlsx"
    path.write_text("irrelevant")
    return str(path)


@pytest.fixture
def results_dataframe():
    """Minimal results DataFrame for OutputProcessor tests."""
    return pd.DataFrame({
        "subtype": ["A", "B"],
        "phenotype": ["sick", "healthy"],
        "odds_ratio": [1.5, 0.8],
        "pvalue": [0.04, 0.21],
    })



def test_data_processor_reads_csv(csv_file):
    '''Function to test the behaviour of the data processor function.'''
    
    df, phenotype = DataProcessor(csv_file, 'phenotype')

    assert isinstance(df, pd.DataFrame)
    assert len(df) == 5
    assert 'T type' in df.columns
    assert 'S type' in df.columns
    assert 'phenotype' not in df.columns
    assert len(phenotype) == 5

def test_data_processor_reads_tsv(tsv_file):
    '''Function to test the behaviour of the data processor function.'''
    
    df, phenotype = DataProcessor(tsv_file, 'phenotype')

    assert isinstance(df, pd.DataFrame)
    assert len(df) == 5
    assert 'T type' in df.columns
    assert 'S type' in df.columns
    assert 'phenotype' not in df.columns
    assert len(phenotype) == 5

def test_data_processor_unsupported_extension(unsupported_extension_file):
    '''Function to test the behaviour of the data processor function when an unsupported file extension is provided.'''
    
    with pytest.raises(ValueError, match="Input file must be in CSV or TSV format"):
        DataProcessor(unsupported_extension_file, 'phenotype')


def test_data_processor_missing_phenotype(csv_file):
    '''Function to test the behaviour of the data processor function when the specified phenotype column is missing.'''
    
    with pytest.raises(ValueError, match="Phenotype column 'missing_column' not found in input file"):
        DataProcessor(csv_file, 'missing_column')

def test_output_processor_saves_csv(tmp_path, results_dataframe):
    '''Function to test the behaviour of the output processor function.'''
    
    output_file = tmp_path
    OutputProcessor(results_dataframe, output_file)

    saved_file = tmp_path / "results.csv"
    assert saved_file.exists()

    loaded_df = pd.read_csv(saved_file)
    pd.testing.assert_frame_equal(loaded_df, results_dataframe)


def test_output_processor_saves_csv_with_prefix(tmp_path, results_dataframe):
    '''Function to test the behaviour of the output processor function.'''

    output_file = tmp_path / "run"
    OutputProcessor(results_dataframe, output_file)

    saved_file = tmp_path / "run_results.csv"
    assert saved_file.exists()

    loaded_df = pd.read_csv(saved_file)
    pd.testing.assert_frame_equal(loaded_df, results_dataframe)