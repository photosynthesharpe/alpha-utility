"""
Split a set of AlphaFold inputs into 1000 fold sets.

Assumes file_list has already been generated

Author: Serena G. Lotreck
"""
import argparse
from os.path import abspath


def main(json_dir):

    # Read in file list
    with open(f'{json_dir}/file_list') as f:
        files = [l.strip() for l in f.readlines()]

    # Split into 1000s
    subsets = [files[i:i + 1000] for i in range(0, len(files), 1000)]

    # Write out
    for i, sub in enumerate(subsets):
        with open(f'{json_dir}/file_list_{i}', 'w') as f:
            f.write('\n'.join(sub))
    print(f'Wrote out {i} subsets as {json_dir}/file_list_#')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Subset folds')

    parser.add_argument('json_dir')

    args = parser.parse_args()

    args.json_dir = abspath(args.json_dir)

    main(args.json_dir)