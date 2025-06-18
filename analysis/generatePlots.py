"""
Make basic graphs of AlphaFold outputs.

Author: Josh Vermaas, adapted by Serena G. Lotreck
"""
import argparse
from os.path import abspath
import pandas as pd
import matplotlib
import numpy as np
import mpl_scatter_density
import matplotlib.pyplot as plt
import networkx as nx
from tqdm import tqdm
import warnings
try:
    from astropy.convolution import Gaussian2DKernel, convolve
    astro_smooth = True
except ImportError as IE:
    astro_smooth = False


def make_contour_plot(data, outpath, outprefix):
    """
    Make the contour plot of log2 RMSD sum vs. mean confidence.

    parameters:
        data, pandas df: filtered data

    retunrs: None
    """
    x = data['meanconfidence']
    y = np.log10(data['rmsdsum'])
    H, xedges, yedges = np.histogram2d(x,
                                       y,
                                       bins=(20, 20),
                                       range=[[0, 1], [0.5, 3]])
    xmesh, ymesh = np.meshgrid(xedges[:-1], yedges[:-1])
    H = np.log10(H)
    # Smooth the contours (if astropy is installed)
    if astro_smooth:
        kernel = Gaussian2DKernel(0.1)
        H = convolve(H, kernel)

    fig, ax = plt.subplots(1, figsize=(7, 6))
    clevels = ax.imshow(H.T,
                        extent=[0, 1, 0.5, 3],
                        origin='lower',
                        cmap='spring',
                        alpha=0.5,
                        zorder=90,
                        aspect='auto')
    ax.scatter(x, y, color='k', s=0.1, zorder=2)
    ax.set_xlabel("AlphaFold Confidence")
    ax.set_ylabel("log10 RMSD Sum ($\\AA$)")

    fig.colorbar(clevels, label="Log10 (Counts)")
    fig.savefig(f'{outpath}/{outprefix}_confidence_contour.png',
                format='png',
                dpi=600)


def main(simple_analysis, outpath, outprefix, network_graph,
         include_all_pairs):

    # Initialize graph
    print('\nInitializing graph...')
    G = nx.Graph()

    # Load data
    print('\nLoading data...')
    simple_data = pd.read_csv(simple_analysis)
    print(f'Snapshot of the data: {simple_data.head()}')

    # Filter data by those that have a confidence of higher than 0.5 or an
    # RMSD sum that is less than 10^1.5
    if not include_all_pairs:
        print('\nFiltering data...')
        filtered_data = simple_data[(simple_data['meanconfidence'] > 0.5) | (
            simple_data['rmsdsum'] < 10**1.5)].set_index('name')
        print(
            f'{len(simple_data) - len(filtered_data)} pairs were removed upon filtering.'
        )
    else:
        filtered_data = simple_data.set_index('name')

    # Iterate over the pairs and add them to the graph
    print('\nAdding nodes and edges to the network...')
    for pair in tqdm(filtered_data.index):
        try:
            n1, n2 = pair.split("_")
        except ValueError:
            print(f'Too many underscores in name {pair}, skipping')
            continue
        G.add_node(n1)
        G.add_node(n2)
        G.add_edge(n1,
                   n2,
                   confidence=filtered_data.loc[pair, 'meanconfidence'],
                   RMSD=filtered_data.loc[pair, 'rmsdsum'])
    gexf_path = f'{outpath}/{outprefix}_network.gexf'
    print(f'Saving as a gexf to: {gexf_path}')
    nx.write_gexf(G, gexf_path)

    # Make a graph visual
    if network_graph:
        print('\nGenerating a graph visual...')
        pos = nx.forceatlas2_layout(G,
                                    scaling_ratio=20,
                                    linlog=True,
                                    max_iter=2000)
        nx.draw_networkx(G, pos=pos)
        # Set margins for the axes so that nodes aren't clipped
        ax = plt.gca()
        ax.margins(0.20)
        plt.axis("off")
        plt.savefig(f'{outpath}/{outprefix}_graph.png', format='png', dpi=600)

    # Make the heatmap
    print('\nGenerating confidence heatmap...')
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    with warnings.catch_warnings(
    ):  # See https://github.com/astrofrog/mpl-scatter-density/issues/35
        warnings.simplefilter('ignore')
        ax.scatter_density(filtered_data['meanconfidence'],
                           np.log10(filtered_data['rmsdsum']))
    ax.set_xlabel("AlphaFold Confidence")
    ax.set_ylabel("log10 RMSD Sum ($\\AA$)")
    fig.savefig(f'{outpath}/{outprefix}_confidence_heatmap.png',
                format='png',
                dpi=600)

    # Make the countour plot
    print('\nGenerating contour plot...')
    make_contour_plot(filtered_data, outpath, outprefix)

    print('\nDone!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Make network of proteins')

    parser.add_argument(
        'simple_analysis',
        type=str,
        help='Path to file that is output from getRMSDmatrix.py')
    parser.add_argument('outpath', type=str, help='Directory to save outputs')
    parser.add_argument('outprefix',
                        type=str,
                        help='string to prepend to output filenames')
    parser.add_argument('--network_graph',
                        action='store_true',
                        help='Whether or not to save a plot of the network '
                        'interactions. This uses the ForceAtlas2 layout, '
                        'which may fail on certain types of networks. The '
                        '.gexf file can be loaded into Gephi for interactive '
                        'visualization')
    parser.add_argument('--include_all_pairs',
                        action='store_true',
                        help='Whether or not to filter by confidence')

    args = parser.parse_args()

    args.simple_analysis = abspath(args.simple_analysis)
    args.outpath = abspath(args.outpath)

    main(args.simple_analysis, args.outpath, args.outprefix,
         args.network_graph, args.include_all_pairs)
