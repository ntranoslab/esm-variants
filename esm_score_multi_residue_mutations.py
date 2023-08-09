# Python standard libraries
import argparse

# Third-party libraries
import pandas as pd
import torch

# Local applications/libraries
from esm_variants_utils import load_esm_model, get_PLLR

def main(args):
    """
    Execute the main script logic.
    
    Parameters:
    args (Namespace): Arguments parsed from command line input.
    """
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print('Using {}.'.format('GPU' if device == 'cuda' else 'CPU (this may be much slower)'))

    input_df = pd.read_csv(args.input_csv_file)

    print('Loading the model ({})...'.format(args.model_name))
    model, alphabet, batch_converter, _ = load_esm_model(args.model_name, device)

    print('Invoking the model...')
    PLLRs = []
    for _, row in input_df.iterrows():
        wt_seq = row['wt_seq']
        mut_seq = row['mut_seq']
        start_pos = int(row['start_pos'])
        PLLR = get_PLLR(wt_seq, mut_seq, start_pos, model, alphabet, batch_converter, device=device)
        PLLRs.append(PLLR)

    print('Saving results...')
    results_df = input_df.copy()
    results_df['esm_score'] = PLLRs
    results_df.to_csv(args.output_csv_file, index=False)

    print('Done.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute ESM effect scores for specified multi-residue variants in a set of protein sequences.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--input-csv-file', dest='input_csv_file', required=True, metavar='/path/to/input_mutations.csv', help='Path to the input CSV file with the protein mutations to calculate ESM scores for.')
    parser.add_argument('--output-csv-file', dest='output_csv_file', required=True, metavar='./esm_multi_residue_effect_scores.csv', help='Path to save the output CSV file.')
    parser.add_argument('--model-name', dest='model_name', default='esm1b_t33_650M_UR50S', metavar='esm1b_t33_650M_UR50S', help='Name of the ESM model to use. See list of options here: https://github.com/facebookresearch/esm#available')

    args = parser.parse_args()

    main(args)
