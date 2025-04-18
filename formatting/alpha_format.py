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


def generate_seq_dicts(group, fasta_sequences, anchor_sequences=None):
    """
    Generate sequence dicts from all desired combinations.

    parameters:
        group, tuple: group of proteins, ligands and cofactors
        fasta_sequences, SeqIO index object: protein sequences to compare
        anchor_sequences, SeqIO index object, optional: anchor protein
            sequences for one_v_many comparisons

    returns:
        instance_dict, dict: the formatted json for this entry
    """
    # I mainly split this off for aesthetic reasons, as with all the comments the
    # function was getting a little unwieldy in terms of lines, but it also allows
    # us to unit test this section separately from the section that gets the basic
    # elements we use to pair and generates the combinations.
    
    instance_dict = {}
    instance_dict['name'] = '_'.join(
        group)  # Could use some finessing for aesthetics

    # We just have to go over every element and make a dict for it. Proteins
    # are a little different so we'll use a full loop as opposed to a
    # (typically way prettier) list comp
    instance_dict['sequences'] = []
    for item in group:
        # There could be other underscores in an ID, so we need to safeguard
        # when splitting the string to get the prepended string back.
        id_str = item.split('_')[0]
        name = '_'.join(item.split('_')
                        [1:])  # Rather than just the second item of the split
        # If protein, we need the 'sequence' key rather than 'ccdCodes'
        if id_str == 'P':
            # Look for the sequence in the protein list
            try:
                seq_dict = {
                    'protein': {
                        'id':
                        id_str,  ## TODO figure out what is actually supposed to go here
                        'sequence':str(fasta_sequences[name].seq)
                    }
                }
            # If using one_v_many, it could fail because the protein is in
            # the anchor sequences, so we catch that error and try the other
            except KeyError:
                seq_dict = {
                    'protein': {
                        'id': id_str, ## TODO figure out what is actually supposed to go here
                        'sequence': str(anchor_sequences[name].seq)
                    }
                }
        # Both cofactor and ligand use the keys 'ligand' and 'ccdCodes'
        else:
            seq_dict = {'ligand': {"id": id_str, "ccdCodes": [name]}}
        # Add it to the list of sequences for this dict
        instance_dict['sequences'].append(seq_dict)

        # Add additional required arguments
        ## TODO I hardcoded these because we're all going to be using AlphaFold3
        ## but you could easily take them as an argument with argparse and pass
        ## them through to here
        instance_dict['dialect'] = 'alphafold3'
        instance_dict['version'] = 2
        # Re: this next param, no clude how many structures people want, def want
        # to update this to be passable
        instance_dict['modelSeeds'] = [1855]

    return instance_dict


def generateJSONs(fasta_sequences,
                  ligand_list=None,
                  cofactor_list=None,
                  anchor_sequences=None,
                  protein_comparison_type='folding_only',
                  make_self_matches=False):
    """
    Make paired json files for ligands and proteins.

    parameters:
        fasta_sequences, SeqIO index object: protein sequences to compare
        ligand_list, list of str, optional: ligand CCD codes
        cofactor_list, list of str, optional: cofactor CCD codes
        anchor_sequences, SeqIO index object, optional: anchor protein
            sequences for one_v_many comparisons
        protein_comparison_type, str: what type of comparison is being made
        make_self_matches, bool: whether or not to include self matches in
            the many_v_many case

    returns:
        inputs, dict of dict: keys are pair names, values are formatted json
            inputs
    """
    # The challenge of this code is to try and be as efficient as possible with
    # making a bunch of 3- or 4-way combinations. The most semantically
    # straightforward way to do this would be a bunch of nested lists, where you
    # go over the proteins, then the ligands and then the cofactors. However, a
    # cleaner (syntactically, I am actually not sure if it's more efficient
    # computationally speaking) way to do this would be to use itertools.product
    # over a group of lists to make tuples of (protein1, protein2 (optional),
    # ligand, cofactor). The only challenge there is, how do we keep track of
    # which CCD codes are ligands vs. cofactors when we do that? Here, we'll pre-
    # pend each ligand with 'L_' and each cofactor with 'C_', as L and C are the
    # id strings we need anyway for the json. We will also put 'P_' on the
    # beginning of proteins, as they require an id of P as well.
    protein_list = ['P_' + prot for prot in fasta_sequences.keys()]
    if cofactor_list is not None:
        cofactor_list = ['C_' + cof for cof in cofactor_list]
    if ligand_list is not None:
        ligand_list = ['L_' + lig for lig in ligand_list]
    if anchor_sequences is not None:
        anchor_list = ['P_' + prot for prot in anchor_sequences.keys()]
    else:
        anchor_list = None

    # The only thing preventing product from being a truly pretty implementation
    # is that passing an empty list as one of the iterables to product makes it
    # return an empty list. Rather than making different calls with if statements,
    # I'll make a list of lists that just excludes anything empty, and pass that.
    product_components = [
        l for l in [protein_list, ligand_list, cofactor_list, anchor_list]
        if (l is not None) and l != []
    ]
    if protein_comparison_type == 'many_v_many':
        # Add the protein list again so we can get all the protein combinations.
        # Note how we're inserting the protein list at the front, next to the
        # other protein list (rather than appending on the end) -- this is on
        # purpose, to simplify dealing with repeated protein combinations 
        product_components.insert(0, protein_list)

    # Get the combinations and iterate to make jsons. We don't make this a list
    # first because the iterator only holds one item in memory at a time and is
    # more efficient. tqdm will show progress but not percentage of total, because
    # it can't know without storing the whole thing in memory
    inputs = {}
    protein_combo_tracker = []
    for combo in tqdm(product(*product_components)): # Have to unpack the list of lists with *
        # The one problem with this approach is that when we provide the list of
        # protein sequences twice for many_v_many, we will get both orders of
        # every combination; i.e. (proteinA, proteinB) and (proteinB, proteinA)
        # will both appear in the results. We don't actually want this behavior,
        # because from AlphaFold's perspective those are the same pair, and it
        # will double the size of what you end up running. So here, we will check
        # if the reverse of the protein pair is already in the results, and if
        # so, skip it. The one problem with this method is that if there are
        # ligands or cofactors, we need to ignore them; this is partially solved
        # by having put the protein lists together in the input to product, as
        # we can then just look at the first two things, but we still need to
        # include the ligand and cofactor or else risk tossing things we actually
        # still need
        if protein_comparison_type == 'many_v_many':
            if (combo[1], combo[0],) + combo[2:] in protein_combo_tracker:
                continue
            # Also use this opportunity to skip self matches if needed
            if not make_self_matches:
                if combo[0] == combo[1]:
                    continue
        inp = generate_seq_dicts(combo, fasta_sequences, anchor_sequences)
        inputs['_'.join(combo)] = inp
        protein_combo_tracker.append(combo)

    return inputs


def main(fasta, out_loc, outprefix, ligands_file, cofactor_file, anchor_fasta,
         protein_comparison_type, make_self_matches):

    # Read in the provided input files
    print('\nReading inputs...')

    # Read the proteins in with index instead of parse, because it allows us to
    # get sequences using the protein ID as a key
    fasta_sequences = SeqIO.index(fasta, 'fasta')
    print(f'There are {len(fasta_sequences)} fasta sequences.')

    if ligands_file is not None:
        with open(ligands_file) as f:
            # Make sure there are no trailing whitespaces
            ligands = [l.strip() for l in f.readlines()]
        print(f'There are {len(ligands)} ligands.')
    else:
        ligands = None

    if cofactor_file is not None:
        with open(cofactor_file) as f:
            # Make sure there are no trailing whitespaces
            cofactors = [l.strip() for l in f.readlines()]
        print(f'There are {len(cofactors)} cofactors.')
    else:
        cofactors = None

    if anchor_fasta is not None:
        anchor_fasta_sequences = SeqIO.index(anchor_fasta, 'fasta')
        print(
            f'There are {len(anchor_fasta_sequences)} anchor fasta sequences.')
    else:
        anchor_fasta_sequences = None

    # Call the function to generate the jsons
    print('\nGenerating jsons...')
    inputs = generateJSONs(fasta_sequences, ligands, cofactors,
                           anchor_fasta_sequences, protein_comparison_type,
                          make_self_matches)

    # Save the output
    print('\nSaving jsons...')
    for specific_name, inp in inputs.items():
        with open(f'{out_loc}/{outprefix}_{specific_name}.json', 'w') as f:
            json.dump(inp, f)

    print('\nDone!')


if __name__ == "__main__":

    # Define our parser
    parser = argparse.ArgumentParser(description='Format data for AlphaFold3')

    # Define the command line args that the parser will accept
    parser.add_argument('fasta',
                        type=str,
                        help='Path to fasta file with protein sequences')
    parser.add_argument('out_loc',
                        type=str,
                        help='Destination directory for the output JSON files')
    parser.add_argument(
        'outprefix',
        type=str,
        help=
        'The experiment being preformed, will prepend every output filename')
    parser.add_argument(
        '-ligands_file',
        default=None,
        type=str,
        help='Path to to file with the CCD ligand codes. Should be a .txt file '
        'with one code per line')
    parser.add_argument(
        '-cofactor_file',
        default=None,
        type=str,
        help='Path to file with the relevant cofactors. Should be a .txt file '
        'with one code per line')
    parser.add_argument(
        '-anchor_fasta',
        default=None,
        type=str,
        help=
        'Path to a secondary fasta file for one- or few-vs-many comparisons')
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
    parser.add_argument(
        '--make_self_matches',
        action='store_true',
        help='Whether or not to include each protein against itself in the '
        'many_v_many scenario' 
    )

    # Get the absolute paths of files and directories
    args = parser.parse_args()
    args.fasta = abspath(args.fasta)
    args.out_loc = abspath(args.out_loc)

    # We want to leave the optional path args as None if not specified,
    # Using abspath would throw an error
    if args.ligands_file is not None:
        args.ligands_file = abspath(args.ligands_file)
    if args.cofactor_file is not None:
        args.cofactor_file = abspath(args.cofactor_file)
    if args.anchor_fasta is not None:
        args.anchor_fasta = abspath(args.anchor_fasta)

    main(args.fasta, args.out_loc, args.outprefix, args.ligands_file,
         args.cofactor_file, args.anchor_fasta, args.protein_comparison_type,
        args.make_self_matches)
