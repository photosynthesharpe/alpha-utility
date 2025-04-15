"""
Format jsons for alphafold.

Authors: Serena & Luke
"""
import argparse
from os.path import abspath
import json
from itertools import product
from Bio import SeqIO
# Some fasta shit


def generateJSONs(protein_fasta, ligand_list, cofactor_list=None):
    """
    Make paired json files for ligands and proteins.

    parameters:
        protein_fasta, str: path to a fasta file containing protein sequences to compare, file path
        ligand_list, list of str: ligand CCD codes
        cofactor_list, list of str, optional: cofactor CCD codes

    returns:
        inputs, dict of dict: keys are pair names, values are formatted json inputs
    """
    # Define empty list if no cofactors
    if cofactor_list is None:
        cofactor_list = []
    
    # Read the protein fasta
    ## do something --> proteins
    # proteins  is a dict of format {'proteinID': 'sequence'}
    # implement a fasta reader that 
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
proteins = {}
for fasta in fasta_sequences:
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

def main(fasta, ligands_file, cofactor_file, out_loc, outprefix):

    print('\nReading inputs...')
    # Read the files we need
    ## ligands
    ## cofactors

    print('\nGenerating jsons...')
    inputs = generateJSONs(protein_fasta, ligand_list, cofactor_list)

    print('\nSaving jsons...')
    for specific_name, inp in inputs.items():
        with open(f'{out_loc}/{out_prefix}_{specific_name}.json', 'w') as f:
            json.dump(inp, f)

    print('\nDone!')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Format for AlphaFold3')

    parser.add_argument(
        'fasta',
        type=str,
        help='Path to fasta file with protein sequences'
    )
    paraser.add_argument(
        'ligands_file',
        type=str,
        help='Path to to file with the CCD ligand codes'
    )
    parser.add_argument(
        'cofactor_file',
        type=str,
        help='Path to file with the relevant cofactors, can be used or not'
    )
    parser.add_argument(
        'out_loc',
        type=str,
        help='Destination of the output JSON files'
    )
    parser.add_argument(
        'outprefix',
        type=str,
        help='The experiment being preformed, will lead every output JSON file'

    # blah blah blah

    args = parser.parse_args()

    args.fasta = abspath(args.fasta)
    args.ligands_file = abspath(args.ligands_file)
    args.cofactor_file = abspath(args.cofactor_file)
    args.out_loc = abspath(args.out_loc)
    # args.outprefix = abspath(args.outprefix)
    # do this for anyt arg that's a path
    # i did this for everything except the output prefix, i believe that it is the only non-path element

    main(args.fasta, args.ligands_file, args.cofactor_file, args.out_loc)