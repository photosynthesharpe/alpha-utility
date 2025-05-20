"""
Split a set of AlphaFold inputs into 1000 fold sets.

Assumes file_list has already been generated

Author: Serena G. Lotreck
"""
import argparse
from os.path import abspath


def main(f_list):

    # Read in file list
    with open(f_list) as f:
        files = [l.strip() for l in f.readlines()]

    # Split into 1000s
    subsets = [files[i:i + 1000] for i in range(0, len(files), 1000)]

    # Write out
    for i, sub in enumerate(subsets):
        with open(f'{f_list}_{i}', 'w') as f:
            f.write('\n'.join(sub))
    print(f'Wrote out {i+1} subsets as {f_list}_#')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Subset folds')

    parser.add_argument('f_list', type=str,
            help='Patht o file containing a list of either files '
            'or directories for an array job')

    args = parser.parse_args()

    args.json_dir = abspath(args.f_list)

    main(args.f_list)
