import argparse
from src.statistical_tests import discrimination_index as di, adjusted_wallace as aw
import pandas as pd


def get_args():
    '''
    
    '''
    parser = argparse.ArgumentParser(description='CAPS is a subtyping comparison command-line tool')

    subtype = parser.add_argument_group('Subtyping method comparison')
    subtype.add_argument('--input', '-i', type=str, required=True, help='Path to the input file (CSV or TSV) containing subtyping data')
    subtype.add_argument('--output', '-o', type=str, required=True, help='Path to the output file where results will be saved')
    subtype.add_argument('--comparator', '-c', default=None, help='name of column (subtyping method) to use as comparator. ONLY NEEDED IF NUMBER OF SUBTYPING METHODS > 2' )

    phenotype = parser.add_argument_group('Phenotype association analysis')
    phenotype.add_argument('--phenotype', '-p', type=str, help='Path to the file (CSV or TSV) containing phenotype data for association analysis')
    phenotype.add_argument_group.add_argument('--subtypes', '-s', help='path to file (CSV or TSV) containing subtypes to perform phenotype association analysis for, can be a single column (one subtyping method) or multiple columns (multiple subtyping methods)')

    return parser.parse_args()

def main():
    args = get_args()

    df = pd.read_csv(args.input, sep='\t' if args.input.endswith('.tsv') else ',')

    DI, DI_low, DI_high = di.DiscriminationIndex(df)

    if args.comparator is not None:
        if args.comparator not in df.columns:
            raise ValueError(f"Comparator column '{args.comparator}' not found in input file")
        comparator_col = df[args.comparator]

        AW_ab, AW_ba = aw.AdjustedWallace(df.iloc[:, 0], comparator_col)
    else:
        AW_ab, AW_ba = aw.AdjustedWallace(df.iloc[: 0], df.iloc[:, 1])


    results = pd.DataFrame({
        'Discrimination Index': [DI],
        '95% CI Lower': [DI_low],
        '95% CI Upper': [DI_high]
    })

    if args.output.endswith('.tsv'):
        results.to_csv(args.output, sep='\t', index=False)
    else:
        results.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()

