"""
Get the RMSD matrix for a set of AlphaFold outputs.

Adapted from Josh Vermaas' PLastidCut analysis.

Author: Serena G. Lotreck
"""
import argparse
from os.path import abspath, dirname
import glob
import pandas as pd
import numpy as np
from vmd import Molecule, atomsel
from tqdm import tqdm


def getRMSDmatrix(filelist):
    """
    Get the RMSD matrix for a given fold.

    Unmodified from PlastidCut analysis.

    parameters:
        filelist, list of str: paths to .cif files

    returns:
        rmsd_matrix, np array: RMSD matrix
    """
    mol = Molecule.Molecule()
    for f in filelist:
        mol.load(f, 'pdbx')
    seli = atomsel("type CA") ## NOTE: this is for py-VMD, you need "name CA"
    selj = atomsel("type CA") ## for a full installation of VMD
    rmsdlist = []
    for i in range(mol.numFrames()):
        seli.frame = i
        for j in range(i + 1, mol.numFrames()):
            selj.frame = j
            rmsd = seli.rmsdQCP(selj)
            rmsdlist.append(rmsd)
    return np.array(rmsdlist)


def main(alphafold_outputs, outpath, out_prefix):

    # Get a list of the files we want
    print('\nIdentifying files to analyze...')
    rankscores = sorted(glob.glob(f'{alphafold_outputs}/*/ranking_scores.csv'))

    # Go through outputs
    print('\nCalculating RMSD sums...')
    rmsd_data = {'name': [], 'meanconfidence': [], 'rmsdsum': []}
    for fold in tqdm(rankscores):
        # Get the fold name from the directory containing these results
        name = fold.split('/')[-2]
        # Load the ranking scores
        data = pd.read_csv(fold)
        # Get the cif files
        ciffiles = sorted(glob.glob(f'{dirname(fold)}/*/model.cif'))
        # Pass to the heavy lifter function
        rmsdsum = np.sum(getRMSDmatrix(ciffiles))
        # Append to dataframe input
        rmsd_data['name'].append(name)
        rmsd_data['meanconfidence'].append(data['ranking_score'].mean())
        rmsd_data['rmsdsum'].append(rmsdsum)

    # Make the dataframe
    print('\nMaking dataframe...')
    rmsd_output = pd.DataFrame(rmsd_data)
    print(f'\nSnapshot of dataframe:\n{rmsd_output.head()}')

    # Save
    print('\nSaving...')
    rmsd_output.to_csv(f'{outpath}/{out_prefix}_rmsd_analysis.csv', index=False)
    print(f'Saved output as {outpath}/{out_prefix}_rmsd_analysis.csv')

    print('\nDone!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get AlphaFold RMSD')

    parser.add_argument('alphafold_outputs',
                        type=str,
                        help='Path to a directory containing the AlphaFold '
                        'outputs. Should contain one directory per fold.')
    parser.add_argument('outpath', type=str, help='Path to save output file')
    parser.add_argument('out_prefix',
                        type=str,
                        help='String to prepend to output filename')

    args = parser.parse_args()

    args.alphafold_outputs = abspath(args.alphafold_outputs)
    args.outpath = abspath(args.outpath)

    main(args.alphafold_outputs, args.outpath, args.out_prefix)
