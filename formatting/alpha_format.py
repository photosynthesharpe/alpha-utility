"""
Process fasta files with ligands and cofactors for input into AlphaFold3.

If more than one cofactor or ligand is passed, the program assumes that a json
should be generated for each ligand/cofactor/protein-protein combination. This
kind of combinatorics can get unwieldy very quickly, so use with caution! The
same goes for how many proteins you provide for the many_v_many or one_v_many
options -- the computational complexity of combinations can get prohibitively
high relatively quickly.

Authors: Serena & Luke
"""
import argparse
from os.path import abspath
import json
from tqdm import tqdm
from itertools import product, combinations
from Bio import SeqIO


def make_folding_only_jsons(protein_sequences, ligand_list=None, cofactor_list=None):
    """
    Makes a set of jsons for AlphaFold where there is one protein and optionally
    a cofactor and/or ligand.

    parameters:
        protein_sequences, SeqIO object: protein sequences to use
        ligand_list, list of str or None: ligands to use
        cofactor_list, list of str or None: cofactors to use

    returns:
        jsons, list of dict: formatted jsons
    """
    for protein


    
    for protein in protein_sequences:
        proteins[fasta.id] = str(fasta.seq)
        
        inputs = {}
        for pair in product(proteins.keys(), ligand_list):
    
            pair_dict = {}
            pair_dict['name'] = '_'.join(pair)
            pair_dict['sequences'] = [
                {
                    'protein': {
                    'id': 'P', # This could also be wrong but idk lol
                    'sequence': proteins[pair[0]]
                    }
                },
                {
                  'ligand': {
                    "id": 'L', ## this could be wrong
                    "ccdCodes": [pair[1]]
                    }
                }
            ]
    
            inputs['_'.join(pair)] = dict(pair_dict)
    
            for cof in cofactor_list:
                pair_dict['sequences'].append(
                    {
                        'ligand': {
                            'id': 'C',
                            'ccdCodes': [cof]
                        }
                    }
                )
                inputs['_'.join(pair) + f'_{cof}'] = dict(pair_dict)
    
        return inputs

def generateJSONs(fasta_sequences, ligand_list=None, cofactor_list=None, anchor_sequences=None, protein_comparison_type='folding_only'):
    """
    Make paired json files for ligands and proteins.

    parameters:
        protein_fasta, str: path to a fasta file containing protein sequences to compare, file path
        ligand_list, list of str: ligand CCD codes
        cofactor_list, list of str, optional: cofactor CCD codes

    returns:
        inputs, dict of dict: keys are pair names, values are formatted json inputs
    """
    # Define empty list if no cofactors or ligands
    if cofactor_list is None:
        cofactor_list = []
    if ligand_list is None:
        ligand_list = []

    # Determine what kind of output we're generating and call the corresponding
    # function. Calling separate functions here allows us to just make a new one
    # if we want a new functionality, as opposed to having to figure out where
    # to put it in a larger body of code with multiple purposes
    if protein_comparison_type == 'folding_only':
        return make_folding_only_jsons(fasta_sequences, ligand_list, cofactor_list)
    elif protein_comparison_type == 'many_v_many':
        return make_many_v_many_jsons(fasta_sequences, ligand_list, cofactor_list)
    elif protein_comparison_type == 'one_v_many':
        return make_one_v_many_jsons(fasta_sequences, ligand_list, cofactor_list, anchor_sequences=anchor_fasta_sequences)


def main(fasta, out_loc, outprefix, ligands_file, cofactor_file, anchor_fasta, protein_comparison_type):

    # Read in the provided input files
    print('\nReading inputs...')
    
    fasta_sequences = SeqIO.parse(open(fasta),'fasta')
    print(f'There are {len(fasta_sequences)} fasta sequences.')

    if ligands_file is not None:
        with open(ligands_file) as f:
            # Make sure there are no trailing whitespaces
            ligands = [l.strip() for l in f.readlines()]
        print(f'There are {len(ligands)} ligands.')

    if cofactor_file is not None:
        with open(cofactor_file) as f:
            # Make sure there are no trailing whitespaces
            cofactors = [l.strip() for l in f.readlines()]
        print(f'There are {len(cofactors)} cofactors.')
        
    if anchor_fasta is not None:
        anchor_fasta_sequences = SeqIO.parse(open(anchor_fasta),'fasta')
        print(f'There are {len(anchor_fasta_sequences)} anchor fasta sequences.')

    # Call the function to generate the jsons
    print('\nGenerating jsons...')
    inputs = generateJSONs(fasta_sequences, ligands, cofactors, anchor_fasta_sequences, protein_comparison_type)

    # Save the output
    print('\nSaving jsons...')
    for specific_name, inp in inputs.items():
        with open(f'{out_loc}/{out_prefix}_{specific_name}.json', 'w') as f:
            json.dump(inp, f)

    print('\nDone!')


if __name__ == "__main__":

    # Define our parser
    parser = argparse.ArgumentParser(description='Format data for AlphaFold3')

    # Define the command line args that the parser will accept
    parser.add_argument(
        'fasta',
        type=str,
        help='Path to fasta file with protein sequences'
    )
    parser.add_argument(
        'out_loc',
        type=str,
        help='Destination directory for the output JSON files'
    )
    parser.add_argument(
        'outprefix',
        type=str,
        help='The experiment being preformed, will prepend every output filename'
    )
    paraser.add_argument(
        '-ligands_file',
        default=None,
        type=str,
        help='Path to to file with the CCD ligand codes. Should be a .txt file '
        'with one code per line'  
    )
    parser.add_argument(
        '-cofactor_file',
        default=None,
        type=str,
        help='Path to file with the relevant cofactors. Should be a .txt file '
        'with one code per line'
    )
    parser.add_argument(
        '-anchor_fasta',
        default=None,
        type=str,
        help='Path to a secondary fasta file for one- or few-vs-many comparisons'
    )
    parser.add_argument(
        '-protein_comparison_type',
        default='folding_only',
        type=str,
        help='Options are: '
        '    many_v_many, gets every possible pair of proteins. Only requires '
        '        the single required fasta input'
        '    one_v_many, gets every comparison of the protein(s) in anchor_fasta '
        '        with the proteins in the required fasta input. Anchor fasta is '
        '        allowed to contain more than one protein'
        '    folding_only, doesn\'t pair any proteins, only cofactors and ligands'
    )

    # Get the absolute paths of files and directories
    args = parser.parse_args()
    args.fasta = abspath(args.fasta)
    args.out_loc = abspath(args.out_loc)
    
    # We want to leave the optional path args as None if not specified,
    # Using abspath would throw an error
    if args.ligand_file is not None:
        args.ligands_file = abspath(args.ligands_file)
    if args.cofactor_file is not None:
        args.cofactor_file = abspath(args.cofactor_file)
    if args.anchor_fasta is not None:
        args.anchor_fasta = abspath(args.anchor_fasta)

    main(args.fasta, args.out_loc, args.outprefix, args.ligands_file, args.cofactor_file, args.anchor_fasta, args.protein_comparison_type)