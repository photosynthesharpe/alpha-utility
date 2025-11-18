"""
Process fasta files with ligands and cofactors for input into AlphaFold3.

If more than one cofactor or ligand is passed, the program assumes that a json
should be generated for each ligand/cofactor/protein-protein combination. This
kind of combinatorics can get unwieldy very quickly, so use with caution! The
same goes for how many proteins you provide for the many_v_many or one_v_many
options -- the computational complexity of combinations can get prohibitively
high relatively quickly.

This attempts to generalize the includion of multi-subunit proteins; however,
I have not exhaustively tested all the ways they could appear in the data, so
use with caution.

Authors: Serena & Luke
"""
import argparse
from os.path import abspath
import json
from tqdm import tqdm
from itertools import product, combinations
from Bio import SeqIO
from string import ascii_uppercase


def format_protein_seqs(name,
                        fasta_sequences,
                        anchor_sequences=None,
                        multi_subunit_proteins=None,
                        in_complex=None,
                        prev_letters=0,
                       remove_trailing_asterisks=False):
    """
    Format the sequence dicts for proteins, including multi-subunit proteins.

    parameters:
        name, str: protein ID
        fasta_sequences, SeqIO index object: protein sequences to compare
        anchor_sequences, SeqIO index object, optional: anchor protein
            sequences for one_v_many comparisons
        multi_subunit_proteins, dict, optional: keys are ID's for proteins with
            multiple subunits, values are dicts with fasta ID: number of
            occurrences pairs
        in_complex, dict: keys are ID's that are in a complex, values are their
            complex ID
        prev_letters, int, default 0: how many sequences are already in the
            dict -- can't repeat letters of the alphabet across subunits or
            pairs of sequences
        remove_trailing_asterisks, bool: whether or not to remove trailing
            asterisks from protein sequences

    returns:
        protein_seqs, list of dict: protein sequences
        total_letters
    """
    # So we can check use if statement later
    if in_complex is None:
        in_complex = {}

    # Turn alphabet string into a list
    alphabet = [letter for letter in ascii_uppercase]
    
    # Figure out which dict has the protein --
    # We're going to assume that all the subunits appear in one or the
    # other of the two fasta sequence files
    if name in fasta_sequences:
        has_seqs = fasta_sequences
    else:
        has_seqs = anchor_sequences

    # Check if it has multiple subunits
    if name in in_complex:
        complex_name = in_complex[name]
        # Now we need to fish out all the subunits
        seqs = []
        total_letters = prev_letters  # So we can do the ID letters without repeats
        for sub, num in multi_subunit_proteins[complex_name].items():
            id_letters = alphabet[total_letters:total_letters + num]
            # Need to start doubling letters if we get past 26
            alphabet_iter = total_letters // 26
            # I am only implementing the ability to have as many sequences as double
            # letters can allow for, because that is a lot, I can't imagine somebody
            # doing more! But if I'm wrong, you'll have to change this implementation
            assert alphabet_iter < 26, 'Too many proteins! Ran out of letters!'
            if alphabet_iter > 0:
                id_letters = [letter + alphabet[alphabet_iter - 1] for letter in id_letters]
            total_letters += num
            if remove_trailing_asterisks:
                seq = str(has_seqs[sub].seq)
                try:
                    aster_idx = seq.index('*')
                    if aster_idx == len(seq) - 1:
                        seq = seq[:-1]
                except ValueError:
                    seq = seq
            else:
                seq = str(has_seqs[sub].seq)
            sub_seq = {
                'protein': {
                    'id': id_letters,
                    'sequence': seq
                }
            }
            seqs.append(sub_seq)

    # If not, just make one sequence, but put in a list for consistency
    else:
        seq = str(has_seqs[name].seq)
        if remove_trailing_asterisks:   
            try:
                aster_idx = seq.index('*')
                if aster_idx == len(seq) - 1:
                    seq = seq[:-1]
            except ValueError:
                seq = seq
        seqs = [{
            'protein': {
                'id': [alphabet[prev_letters]],
                'sequence': seq
            }
        }]
        total_letters = 1
    
    return seqs, total_letters


def generate_seq_dicts(group,
                       fasta_sequences,
                       anchor_sequences=None,
                       multi_subunit_proteins=None,
                      remove_trailing_asterisks=False):
    """
    Generate sequence dicts from all desired combinations.

    parameters:
        group, tuple: group of proteins, ligands and cofactors
        fasta_sequences, SeqIO index object: protein sequences to compare
        anchor_sequences, SeqIO index object, optional: anchor protein
            sequences for one_v_many comparisons
        multi_subunit_proteins, dict, optional: keys are ID's for proteins with
            multiple subunits, values are dicts with fasta ID: number of
            occurrences pairs
        remove_trailing_asterisks, bool: whether or not to remove trailing
            asterisks from protein sequences

    returns:
        instance_dict, dict: the formatted json for this entry
    """
    # I mainly split this off for aesthetic reasons, as with all the comments the
    # function was getting a little unwieldy in terms of lines, but it also allows
    # us to unit test this section separately from the section that gets the basic
    # elements we use to pair and generates the combinations.

    # Reindex multi subunit proteins for future reference
    if multi_subunit_proteins is not None:
        in_complex = {
            f_id: comp_id
            for comp_id, members in multi_subunit_proteins.items()
            for f_id, num in members.items()
        }
    else:
        in_complex = {}

    instance_dict = {}
    instance_name_parts = []

    # We just have to go over every element and make a dict for it. Proteins
    # are a little different so we'll use a full loop as opposed to a
    # (typically way prettier) list comp
    instance_dict['sequences'] = []
    total_letters = 0
    for item in group:
        # There could be other underscores in an ID, so we need to safeguard
        # when splitting the string to get the prepended string back.
        seq_type = item.split('_')[0]
        name = '_'.join(item.split('_')
                        [1:])  # Rather than just the second item of the split
        # If protein, we need the 'sequence' key rather than 'ccdCodes'
        if seq_type == 'P':
            if name in in_complex:
                instance_name_parts.append(in_complex[name])
            else:
                instance_name_parts.append(name)
            seqs, update_total_letters = format_protein_seqs(name, fasta_sequences,
                                                      anchor_sequences,
                                                      multi_subunit_proteins,
                                                      in_complex,
                                                      total_letters,
                                                      remove_trailing_asterisks)
            total_letters += update_total_letters
        # Both cofactor and ligand use the keys 'ligand' and 'ccdCodes'
        else:
            # Check if we're putting multiple ligands/cofactors in the same fold
            if len(name.split('_')) > 1:
                all_ligands = name.split('_')
            else:
                all_ligands = [name]
            seqs = [{'ligand': {"id": [ascii_uppercase[total_letters]], "ccdCodes": [n]}} for n in all_ligands]
            total_letters += 1
            instance_name_parts.extend(all_ligands)
        # Add it to the list of sequences for this dict
        instance_dict['sequences'].extend(seqs)

        # Add additional required arguments
        instance_dict['name'] = '_'.join(instance_name_parts)
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
                  group_ligands=False,
                  cofactor_list=None,
                  group_cofactors=False,
                  anchor_sequences=None,
                  protein_comparison_type='folding_only',
                  multi_subunit_proteins=None,
                  make_self_matches=False,
                 remove_trailing_asterisks=False):
    """
    Make paired json files for ligands and proteins.

    parameters:
        fasta_sequences, SeqIO index object: protein sequences to compare
        ligand_list, list of str, optional: ligand CCD codes
        cofactor_list, list of str, optional: cofactor CCD codes
        anchor_sequences, SeqIO index object, optional: anchor protein
            sequences for one_v_many comparisons
        protein_comparison_type, str: what type of comparison is being made
        multi_subunit_proteins, dict, optional: keys are ID's for proteins with
            multiple subunits, values are dicts with fasta ID: number of
            occurrences pairs
        make_self_matches, bool: whether or not to include self matches in
            the many_v_many case
        remove_trailing_asterisks, bool: whether or not to remove trailing
            asterisks from protein sequences

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
    # pend each ligand with 'L_' and each cofactor with 'C_'. We will also put
    # 'P_' on the beginning of proteins.

    # An additional consideration for wanting to put multiple ligands or cofactors
    # in the same fold is that we need to treat all of them like one item, but then
    # unpack them later. In order to preserve compatibility with the way I identify
    # the components in the generate_seq_dicts function (splitting on '_'), I'll use
    # a delimiter between the different cofactors. CCD codes are 1-3 alphanumeric
    # characters, so we should beable to use any reasonable seperator; however, I'll
    # check for the delimiter character first, because https://xkcd.com/327/
    protein_list = ['P_' + prot for prot in fasta_sequences.keys()]
    if cofactor_list is not None:
        if group_cofactors:
            for cof in cofactor_list:
                assert '_' in cof, 'Please remove _ characters from cofactors'
            cofactor_list = ['C_' + '_'.join(cofactor_list)]
        else:
            cofactor_list = ['C_' + cof for cof in cofactor_list]
    if ligand_list is not None:
        if group_ligands:
            for lig in ligand_list:
                assert '_' in lig, 'Please remove _ characters from ligands'
            ligand_list = ['L_' + '_'.join(ligand_list)]
        else:
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
    for combo in tqdm(product(
            *product_components)):  # Have to unpack the list of lists with *
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
        # still need. We also have to check for the same order because each subunit
        # will end up with its own comparison
        if protein_comparison_type == 'many_v_many':
            if (
                    combo[1],
                    combo[0],
            ) + combo[
                    2:] in protein_combo_tracker or combo in protein_combo_tracker:
                continue
            # Also use this opportunity to skip self matches if needed
            if not make_self_matches:
                if combo[0] == combo[1]:
                    continue
        inp = generate_seq_dicts(combo, fasta_sequences, anchor_sequences,
                                 multi_subunit_proteins, remove_trailing_asterisks)
        inputs[inp['name']] = inp
        protein_combo_tracker.append(combo)

    return inputs


def main(fasta, out_loc, outprefix, ligands_file, group_ligands, cofactor_file, group_cofactors, anchor_fasta,
         protein_comparison_type, multi_subunit_proteins, make_self_matches,
        remove_trailing_asterisks):

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

    if multi_subunit_proteins is not None:
        with open(multi_subunit_proteins) as f:
            multi_subunit_proteins = json.load(f)

    # Call the function to generate the jsons
    print('\nGenerating jsons...')
    inputs = generateJSONs(fasta_sequences, ligands, group_ligands, cofactors, group_cofactors,
                           anchor_fasta_sequences, protein_comparison_type,
                           multi_subunit_proteins, make_self_matches,
                          remove_trailing_asterisks)

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
    parser.add_argument('--group_ligands',
        action='store_true',
        help='Provides all ligands to every fold, rather than one at a time')
    parser.add_argument(
        '-cofactor_file',
        default=None,
        type=str,
        help='Path to file with the relevant cofactors. Should be a .txt file '
        'with one code per line')
    parser.add_argument('--group_cofactors',
        action='store_true',
        help='Provides all cofactors to every fold, rather than one at a time')
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
        help='Options are:'
        '    many_v_many, gets every possible pair of proteins. Only requires '
        '        the single required fasta input'
        '    one_v_many, gets every comparison of the protein(s) in anchor_fasta '
        '        with the proteins in the required fasta input. Anchor fasta is '
        '        allowed to contain more than one protein'
        '    folding_only, doesn\'t pair any proteins, only cofactors and ligands'
    )
    parser.add_argument(
        '-multi_subunit_proteins',
        default=None,
        type=str,
        help='Path to a json that identifies multi-subunit proteins. The '
        'structure of the json for a single protein with 2 kinds of subunits '
        'is the following: '
        '{"complex_name": {"subunit_1": 8, "subunit_2": 8}}, where complex_name '
        'can be any arbitrary string, but subunit_1 and subunit_2 need to be '
        'valid identifiers in at least one provided fasta, and the values in '
        'the internal dict are the number of times that subunit appears in the '
        'complex. If only one subunit is found, a warning will be raised, but '
        'will not cause an exception.')
    parser.add_argument(
        '--make_self_matches',
        action='store_true',
        help='Whether or not to include each protein against itself in the '
        'many_v_many scenario')
    parser.add_argument(
        '--remove_trailing_asterisks', ##TODO add a unit test for this
        action='store_true',
        help='Whether or not to remove asterisks on the end of protein sequences. '
        '* are stop codons in a fasta, and AlphaFold does not allow them. However, '
        'as caution should be exercised removing asterisks not at the end of '
        'sequences, only trailing asterisks will be removed with this option.'
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
    if args.multi_subunit_proteins is not None:
        args.multi_subunit_proteins = abspath(args.multi_subunit_proteins)

    main(args.fasta, args.out_loc, args.outprefix, args.ligands_file, args.group_ligands,
         args.cofactor_file, args.group_cofactors, args.anchor_fasta, args.protein_comparison_type,
         args.multi_subunit_proteins, args.make_self_matches,
         args.remove_trailing_asterisks)
