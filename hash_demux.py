# -*- coding: utf-8 -*-
"""
mission bio single-cell pipeline code
written by ben 10.8.2020

"""

# modules
import numpy as np
import pandas as pd
import os
import shutil
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.signal import find_peaks
import itertools
import seaborn as sns
from collections import Counter

# plotting
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['axes.unicode_minus'] = False

# load external files
import resources

def demux_cells(ab_clr_count_file, hash_csv, hashing_folder):
    # associate cells with a sample based on their hash

    # create hashing output folder if it doesn't exist
    if not os.path.exists(hashing_folder):
        os.mkdir(hashing_folder)
    else:
        shutil.rmtree(hashing_folder)
        os.mkdir(hashing_folder)

    # extract sample name
    sample_name = os.path.basename(ab_clr_count_file).split('.umi_counts.clr.cells.tsv')[0]

    # load antibody hash descriptions
    hashes = resources.load_barcodes(hash_csv, 1, False)
    hash_names, sample_labels = zip(*hashes.items())

    # load clr counts of each ab
    ab_counts_clr = pd.read_csv(ab_clr_count_file, sep='\t', header=0, index_col=0)

    # determine thresholds for each hash and plot fitted distributions
    plt.figure(figsize=(7, 5))
    thresholds = []
    colors = []
    for h in hash_names:

        # fit double gaussian distributions to the hash count data
        # the distributions should fit to the negative and doublet populations
        try:
            data = ab_counts_clr[h]
        except KeyError:
            print('Hashing Ab not in list of sample Abs. Ignoring this hash...')
            continue

        color = plt.get_cmap('tab10')(hash_names.index(h))
        colors.append(color)
        y, x, _ = plt.hist(data, 100, alpha=.3, facecolor=color, label=h)
        x = (x[1:] + x[:-1]) / 2  # for len(x)==len(y)

        def gauss(x, mu, sigma, A):
            return A * np.exp(-(x - mu) ** 2 / 2 / sigma ** 2)

        def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
            return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)

        # estimate the mean of the distributions using the max of the kde
        kde_x, kde_y = sns.kdeplot(data, bw_method=0.2, lw=0).lines[0].get_data()
        peaks, _ = find_peaks(kde_y, prominence=1e-3)
        if len(peaks) == 2:
            expected = (kde_x[peaks[0]], 1, 50, kde_x[peaks[1]], 1, 50)
        elif len(peaks) == 1:
            expected = (kde_x[peaks[0]], 1, 50, 0.5, 1, 50)
        else:
            expected = (-2, 1, 50, 0.5, 1, 50)

        # fit the double gaussian to the data
        params, cov = curve_fit(bimodal, x, y, expected)

        plt.plot(x, bimodal(x, *params), lw=2, c=color)
        plt.legend()
        plt.xlabel('CLR')
        plt.ylabel('# Cells')

        # set the threshold at the 99th percentile of the negative distribution
        thresh = norm(params[0], params[1]).ppf(0.99)
        thresholds.append(thresh)
        plt.axvline(x=thresh, linestyle=':', lw=2, c=color)
        print(h + ': ' + str(thresh))

    plt.savefig(hashing_folder + sample_name + '.curve_fit.png', dpi=300)

    # assign cells to samples based on thresholds
    # create dataframe indicating whether a cell is positive for each hash
    hash_df = pd.DataFrame()
    for i in range(len(hash_names)):
        hash_df[hash_names[i]] = (ab_counts_clr[hash_names[i]] > thresholds[i])
    hash_df['num_hashes'] = hash_df.sum(axis=1)
    hash_df['sample'] = 'NO HASH'
    hash_df.loc[hash_df['num_hashes'] > 1, 'sample'] = 'MULTIPLET'
    for i in range(len(hash_names)):
        hash_df.loc[(hash_df['num_hashes'] == 1) & (hash_df[hash_names[i]]), 'sample'] = sample_labels[i]

    # count hashes per cell and create summary
    hash_counts = Counter(hash_df['num_hashes'])
    num_cells = sum(hash_counts.values())
    summary = ''
    labels, values = zip(*hash_counts.items())
    labels, values = (list(t) for t in zip(*sorted(zip(labels, values))))

    for n in labels:
        desc = '%d hashes: %d cells (%0.1f%% of cells)' % (n, hash_counts[n], hash_counts[n] / num_cells * 100)
        summary += desc + '\n'

    # plot bar chart of hashes per cell
    plt.figure()
    indexes = np.arange(len(labels))
    width = 0.8
    plt.bar(indexes, values, width)
    plt.xticks(indexes, labels)
    plt.xlabel('# Hashes per Cell', fontsize=12)
    plt.ylabel('# Cells', fontsize=12)
    plt.gcf().text(1.0, 0.4, summary, fontsize=12)
    plt.tight_layout()
    plt.savefig(hashing_folder + sample_name + '.hashes_per_cell.png',
                bbox_inches='tight',
                dpi=300)

    # plot a similar bar chart of cells per sample
    sample_counts = Counter(hash_df['sample'])
    summary = ''
    labels, values = zip(*sample_counts.items())

    for s in labels:
        desc = '%s: %d cells (%0.1f%% of cells)' % (s, sample_counts[s], sample_counts[s] / num_cells * 100)
        summary += desc + '\n'

    plt.figure()
    indexes = np.arange(len(labels))
    width = 0.8
    plt.bar(indexes, values, width)
    plt.xticks(indexes, labels)
    plt.xlabel('Sample', fontsize=12)
    plt.ylabel('# Cells', fontsize=12)
    plt.gcf().text(1.0, 0.4, summary, fontsize=12)
    plt.tight_layout()
    plt.savefig(hashing_folder + sample_name + '.cells_per_sample.png',
                bbox_inches='tight',
                dpi=300)

    # plot all pairs of hashes as scatter plots
    hash_combos = [list(x) for x in itertools.combinations(hash_names, 2)]
    color_dict = dict(zip(sample_labels, colors))
    color_dict['NO HASH'] = 'grey'
    color_dict['MULTIPLET'] = 'k'
    for h in hash_combos:
        plt.figure(figsize=(5, 5))
        ax = sns.scatterplot(x=ab_counts_clr[h[0]],
                             y=ab_counts_clr[h[1]],
                             s=4,
                             hue=hash_df['sample'],
                             palette=color_dict,
                             edgecolor=None)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles, labels=labels)
        plt.axvline(x=thresholds[hash_names.index(h[0])], c='k', linestyle='--')
        plt.axhline(y=thresholds[hash_names.index(h[1])], c='k', linestyle='--')
        plt.xlabel('CLR ' + h[0] + ' (' + hashes[h[0]] + ')')
        plt.ylabel('CLR ' + h[1] + ' (' + hashes[h[1]] + ')')
        plt.savefig(hashing_folder + sample_name + '.' + h[0] + '.' + h[1] + '.scatter.png',
                    bbox_inches='tight',
                    dpi=300)

    # save hash truth table and sample assignments to file
    hash_df.index = hash_df.index + '-' + sample_name
    hash_df.to_csv(path_or_buf=hashing_folder + sample_name + '.sample_hashes.tsv', sep='\t')