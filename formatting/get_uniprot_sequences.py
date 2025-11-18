"""
Pull sequences for a list of uniprot IDs

Author: Serena G. Lotreck
"""
import argparse
from os.path import abspath
import pandas as pd
from tqdm import tqdm
import requests
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


def format_fasta(protein_sequences, fasta_name):
    """
    Format sequences into fasta and save.

    parameters:
        protein_sequences, dict: keys are uniprot IDs, values are sequences
        fasta_name, str: path to save fasta file

    returns: None
    """
    # Make seq record objects
    seq_list = []
    for prot, sequence in protein_sequences.items():
        record = SeqRecord(
            Seq(sequence),
            id=prot)
        seq_list.append(record)

    # Write to disk
    with open(fasta_name, "w") as output_handle:
        SeqIO.write(seq_list, output_handle, "fasta")
    print(f'Fasta file has been saved to {fasta_name}')
    

def pull_sequences(uniprot_id_list):
    """
    Make API requests to UniProt.

    parameters:
        uniprot_id_list, list of str: uniprot IDs to pull

    returns:
        protein_sequences, dict: keys are uniprot IDs, values are sequences
    """
    protein_sequences = {}

    for up_id in tqdm(uniprot_id_list):
        params = {
          "fields": [
            "sequence"
          ]
        }
        headers = {
          "accept": "application/json"
        }
        base_url = f"https://rest.uniprot.org/uniprotkb/{up_id}"
        
        response = requests.get(base_url, headers=headers, params=params)
        if not response.ok:
          response.raise_for_status()
          sys.exit()
        
        data = response.json()
        protein_sequences[data["primaryAccession"]] = data["sequence"]["value"]

    return protein_sequences


def main(uniprot_ids, fasta_name):

    print('\nReading in IDs...')
    uniprot_id_list = pd.read_csv(uniprot_ids).dropna().uniprot.tolist()

    print('\nQuerying UniProt API for sequences...')
    protein_sequences = pull_sequences(uniprot_id_list)

    print('\nFormatting fasta...')
    format_fasta(protein_sequences, fasta_name)

    print('\nDone!')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Pull sequences for a list of uniprot IDs')

    parser.add_argument('uniprot_ids', type=str,
                       help='Path to csv with uniprot IDs, column with IDs must be named '
                       '"uniprot"')
    parser.add_argument('fasta_name', type=str,
                       help='Path to a fasta file to save output')


    args = parser.parse_args()

    args.uniprot_ids = abspath(args.uniprot_ids)
    args.fasta_name = abspath(args.fasta_name)

    main(args.uniprot_ids, args.fasta_name)
    
