import argparse
from src import input_output_processor as iop
from src import statistical_tests as di
from src import adjusted_wallace as aw
from src import phenotype_association as pa
import pandas as pd


def get_args():
    '''
    
    '''
    parser = argparse.ArgumentParser(description='CAPS is a subtyping comparison command-line tool')

    subtype = parser.add_argument_group('Subtyping method comparison')
    subtype.add_argument('--input', '-i', type=str, required=True, help='Path to the input file (CSV or TSV) containing subtyping data and phenotype data for association analysis')
    subtype.add_argument('--output', '-o', type=str, required=True, help='Path to the output file where results will be saved')
    subtype.add_argument('--comparator', '-c', default=None, help='name of column (subtyping method) to use as comparator. ONLY NEEDED IF NUMBER OF SUBTYPING METHODS > 2' )

    phenotype = parser.add_argument_group('Phenotype association analysis')
    phenotype.add_argument('--phenotype', '-p', type=str, help='name of the column containing the phenotype data to be used for association analysis')

    return parser.parse_args()

def main():
    args = get_args()

    df, phenotype = iop.DataProcessor(args.input, args.phenotype)

    DI, DI_low, DI_high = di.DiscriminationIndex(df)

    if args.comparator is not None:
        if args.comparator not in df.columns:
            raise ValueError(f"Comparator column '{args.comparator}' not found in input file")
        comparator_col = df[args.comparator]

        AW_ab, AW_ba = aw.AdjustedWallace(df.iloc[:, 0], comparator_col)
    else:
        AW_ab, AW_ba = aw.AdjustedWallace(df.iloc[: 0], df.iloc[:, 1])


    for i in df.columns:
        print('Testing for perfect seperation in data')
        seperation = pa.SeperationInvestigator(df[i], phenotype)

        if seperation == True:
            print('SEPERATION DETECTED!! Will use Firth logistic regression to compute odds ratios')

        else:
            print('NO SEPERATION DETECTED!! Will compute odds ratios directly from the 2x2 contingency table')

        print(f"Testing for associations between {i} and {args.phenotype}")
        assoc_df = pa.AssociationTest(df[i], phenotype, seperation=seperation)
        iop.OutputProcessor(assoc_df, args.output + f"{i}_association_")


    results = pd.DataFrame({
        'DI': [DI],
        'CI_low': [DI_low],
        'CI_high': [DI_high],
        'AWC_AvsB' : [AW_ab],
        'AWC_BvsA' : [AW_ba]
    })

    iop.OutputProcessor(results, args.output + "DI_AWC_")

if __name__ == "__main__":
    main()

