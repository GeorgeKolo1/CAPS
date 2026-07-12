import pandas as pd
from pathlib import Path

def DataProcessor(input, phenotype):
    '''
    Function to process the input data file and divorce the phenotype column from the DataFrame containing subtyping method data.
    
    Parameters:
    input (str): Path to the input file (CSV or TSV) containing subtyping data and phenotype data for association analysis.
    
    Returns:
    pd.DataFrame: A cleaned DataFrame ready for analysis.
    pd.Series: A cleaned pandas series ready for analysis.
    '''
    if input.endswith('.tsv'):
        df = pd.read_csv(input, sep='\t')
    elif input.endswith('.csv'):
        df = pd.read_csv(input)
    else:
        raise ValueError("Input file must be in CSV or TSV format")
    
    
    if phenotype not in df.columns:
        raise ValueError(f"Phenotype column '{phenotype}' not found in input file")
    else:
        phenotype_series = df[phenotype]
        df = df.drop(columns=[phenotype])
    
    return df, phenotype_series

def OutputProcessor(results_df, outfile):
    '''
    Function to process the results DataFrame and save it to a specified output file.
    
    Parameters:
    results_df (pd.DataFrame): A DataFrame containing the results of the association analysis.
    outfile (str): A string that contains the directory and filename (prefix) in which to save the results; 
    Defaults to current working directory when not specified or left as None (default)'''

    outfile = Path(outfile)
    if outfile.is_dir():
        outpath = outfile / "results.csv"
    else:
        outpath = outfile.with_name(outfile.name + "_results.csv")
   
    results_df.to_csv(outpath, index=False)