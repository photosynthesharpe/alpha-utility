"""
Get the top X interaction candidates based on the ipTM scores for each fold.

Assumes the protein of interest is the first chain; in a monomer-monomer pair
this doesn't matter, but for a multimer, does assume the interactor comes first.

Comparison to a CoIP dataset from C. merolae suggests an optimal X is 5; however,
more characterization is needed to see how best to use this metric.

Author: Serena Lotreck
"""
import argparse
from os import listdir
from os.path import abspath, isdir
import json
import pandas as pd
import matplotlib.pyplot as plt


def main(results_dir, out_loc, out_prefix, num_top):

    print('\nReading in metrics...')
    summary_confidence_dicts = {}
    for elt in listdir(results_dir):
        if isdir(f'{results_dir}/{elt}'):
            try:
                with open(f'{results_dir}/{elt}/{elt}_summary_confidences.json') as f:
                    summary_confidence_dicts[elt] = json.load(f)
            except FileNotFoundError:
                print('Failed run: ', list(listdir(f'{results_dir}/{elt}')))

    print('\nPlotting score distribution...')
    iptm = {k: v['chain_iptm'][0] for k, v in summary_confidence_dicts.items()}
    plt.hist(iptm.values(), color='turquoise', alpha=0.5, label='chain_iptm')
    plt.title('CoIP-precipitated proteins')
    plt.savefig(f'{out_loc}/{out_prefix}_ipTM_histogram.png', format='png', dpi=300, bbox_inches='tight')
 
    print('\nGetting top scorers...')
    top5 = sorted(iptm.items(), key=lambda item: item[1], reverse=True)[:5]
    top5_df = pd.DataFrame(top5, columns=['Comparison', 'ipTM'])
    top5_df.to_csv(f'{out_loc}/{out_prefix}_top{num_top}_pairs.csv', index=False)

    print('\nDone!')

    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Get top X proteins')

    parser.add_argument('alpha_output_topdir', type=str,
                       help='Path to where alphafold outputs are')
    parser.add_argument('out_loc', type=str,
                       help='Where to save outputs')
    parser.add_argument('out_prefix', type=str,
                       help='String to prepend to output files')
    parser.add_argument('-num_top', type=int, default=5,
                       help='Number of top ranked interactors to report')


    args = parser.parse_args()

    args.alpha_output_topdir = abspath(args.alpha_output_topdir)
    args.out_loc = abspath(args.out_loc)
    
    main(args.alpha_output_topdir, args.out_loc, args.out_prefix, args.num_top)