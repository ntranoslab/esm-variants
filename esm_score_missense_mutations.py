# Python standard libraries
import argparse

# Third-party libraries
import pandas as pd
from Bio import SeqIO
import torch

# Local applications/libraries
from esm_variants_utils import load_esm_model, get_wt_LLR

def fasta_to_dataframe(fasta_file):
    """
    Convert a FASTA file into a Pandas DataFrame.
    
    Parameters:
    fasta_file (str): Path of the FASTA file.
    
    Returns:
    df (DataFrame): Pandas DataFrame containing id, gene, seq, and length columns.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    data = {
        "id": [record.id for record in records],
        "gene": [record.description.split()[1] if len(record.description.split()) > 1 else record.id for record in records],
        "seq": [str(record.seq) for record in records],
        "length": [len(record.seq) for record in records]
    }

    df = pd.DataFrame(data)

    return df


def main(args):
    """
    Execute the main script logic.
    
    Parameters:
    args (Namespace): Arguments parsed from command line input.
    """
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print('Using {}.'.format('GPU' if device == 'cuda' else 'CPU (this may be much slower)'))

    input_df = fasta_to_dataframe(args.input_fasta_file)

    print('Loading the model ({})...'.format(args.model_name))
    model, alphabet, batch_converter, repr_layer = load_esm_model(args.model_name, device)

    print('Invoking the model...')
    input_df_ids, LLRs = get_wt_LLR(input_df, model=model, alphabet=alphabet, batch_converter=batch_converter, device=device)

    print('Saving results...')
    results = []
    for seq_id, LLR in zip(input_df_ids, LLRs):
        raw_seq_results = LLR.transpose().stack().reset_index().rename(columns = {'level_0': 'wt_aa_and_pos', 'level_1': 'mut_aa', 0: 'esm_score'})
        seq_results = pd.DataFrame({'seq_id': seq_id, 'mut_name': raw_seq_results['wt_aa_and_pos'].str.replace(' ', '') + raw_seq_results['mut_aa'], 'esm_score': raw_seq_results['esm_score']})
        results.append(seq_results)
        
    results = pd.concat(results).reset_index(drop=True)
    results.to_csv(args.output_csv_file, index=False)

    print('Done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute ESM effect scores for all possible missense variants in a set of protein sequences.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--input-fasta-file', dest='input_fasta_file', required=True, metavar='/path/to/input_aa_seqs.fasta', help='Path to the input FASTA file with the amino-acid sequences to calculate ESM scores for.')
    parser.add_argument('--output-csv-file', dest='output_csv_file', required=True, metavar='./esm_missense_effect_scores.csv', help='Path to save the output CSV file.')
    parser.add_argument('--model-name', dest='model_name', default='esm1b_t33_650M_UR50S', metavar='esm1b_t33_650M_UR50S', help='Name of the ESM model to use. See list of options here: https://github.com/facebookresearch/esm#available')

    args = parser.parse_args()

    main(args)