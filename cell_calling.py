# -*- coding: utf-8 -*-
"""
mission bio single-cell pipeline code
written by ben 11.25.2018

"""

# modules
from __future__ import division
import os
import os.path
import collections
import csv
from itertools import product, combinations, izip, groupby
from multiprocessing import Process, Queue
import json
from collections import Counter
from itertools import izip
import subprocess
import sys
import copy
import numpy as np
import pandas as pd

# peak calling
from scipy.stats import binned_statistic
from scipy.interpolate import UnivariateSpline
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

# plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['axes.unicode_minus'] = False

def call(output_folder, experiment_name, method='second_derivative', threshold=True):
    # call cells using selected method. returns list of valid cell barcodes.
    # if threshold is on, further refine to cells with 'min_coverage' in 'min_fraction' intervals.

    # set parameters for coverage minimums (applied when threshold is True)

    # minimum coverage for an interval
    min_coverage = 8

    # minimum fraction of intervals with min_coverage
    min_fraction = 0.6

    # load alignment df from all barcodes into dataframe
    all_tsv = output_folder + experiment_name + '.all.tsv'  # alignment counts for all barcodes (previously saved)
    all_df = pd.read_csv(all_tsv, sep='\t', header=0, index_col=0)

    # extract barcodes and read totals
    barcodes = list(all_df.index)
    reads_per_cell = [int(i) for i in list(all_df.sum(axis=1))]
    reads_per_cell, barcodes = (list(t) for t in zip(*sorted(zip(reads_per_cell, barcodes), reverse=True)))

    if method == 'second_derivative':
        # this method uses the second derivative (inflection point) of the knee plot to identify cells

        # 1 - first derivative of cell rank plot

        # exclude barcodes with low numbers of reads
        rpc_thresh = [x for x in reads_per_cell if x >= 100]

        x = np.log10(range(1, len(rpc_thresh) + 1))
        y = np.log10(np.array(rpc_thresh))

        dy = np.zeros(y.shape, np.float)
        dy[0:-1] = np.diff(y) / np.diff(x)
        dy[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])

        dy = -dy  # invert for positive graph

        # 2 - smooth with binning + filter

        dy = pd.Series(dy)  # convert to series

        # bin the data
        bin_means, bin_edges, binnumber = binned_statistic(range(len(dy)), dy, statistic='mean', bins=1000)
        bin_width = (bin_edges[1] - bin_edges[0])
        bin_centers = bin_edges[1:] - bin_width / 2

        # smooth the data by filtering and call first peak
        yhat = savgol_filter(bin_means, 151, 3)  # window size, polynomial order
        first_peak = find_peaks(yhat, prominence=0.1)[0][0]

        # first n cell barcodes are valid
        n_cells = int(bin_centers[first_peak])
        cell_barcodes = barcodes[:n_cells]

    elif method == 'simple_minimum':
        # this method uses a simple minimum to call cells

        num_targets = len(all_df.columns) - 1
        min_reads_per_cell = min_coverage * min_fraction * num_targets
        cell_barcodes = [b for b in barcodes if barcodes[b] >= min_reads_per_cell]

    else:
        print 'Invalid method selected. Exiting...'
        raise SystemExit

    # if thresholding is on, further filter cells according to coverage
    if threshold:
        valid_cells = []
        for c in cell_barcodes:
            amplicon_counts = list(all_df.loc[c, : ])
            if len([i for i in range(len(amplicon_counts)) if amplicon_counts[i] >= min_coverage])\
                    / len(amplicon_counts) >= min_fraction:
                valid_cells.append(c)

    else:
        valid_cells = cell_barcodes

    # save cell barcode file
    cell_df = all_df[all_df.index.isin(valid_cells)]
    cell_tsv = output_folder + experiment_name + '.cells.tsv'
    cell_df.to_csv(path_or_buf=cell_tsv, sep='\t')

    # generate knee plot and amplicon boxplots
    kneeplot_path = output_folder + experiment_name + '.kneeplot.png'
    boxplot_path = output_folder + experiment_name + '.amplicon_boxplot.png'

    knee_plot(all_df, cell_df, kneeplot_path)
    amplicon_boxplot(cell_df, boxplot_path)

    return valid_cells

def knee_plot(all_df, cell_df, output_file):
    # generate knee plot with auc colored according to cell association

    reads_per_cell = [int(i) for i in list(all_df.sum(axis=1))]
    n_cells = len(cell_df.index)

    # plot log-log reads per cell vs cells
    reads_per_cell.sort(reverse=True)
    plt.figure(figsize=(7, 7))
    ax = plt.axes()
    ax.grid()
    plt.loglog(range(len(reads_per_cell)), reads_per_cell, color='k', linewidth=1.5)
    plt.axvline(x=n_cells, color='k', linestyle='--', label='Cell Threshold')
    plt.xlabel('Cell Barcode #', fontsize=18, labelpad=12)
    plt.ylabel('Reads per Cell Barcode', fontsize=18, labelpad=12)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    t1 = plt.annotate('N' + r'$_{cells}$' + ' = %d' % n_cells, xy=(2, 15), fontsize=18)
    t1.set_bbox(dict(facecolor=None, edgecolor=None, alpha=0))

    per_valid = cell_df.to_numpy().sum() / all_df.to_numpy().sum() * 100
    t2 = plt.annotate('%% reads assigned\nto valid cells = %0.1f%%' % per_valid, xy=(2, 3), fontsize=18)
    t2.set_bbox(dict(facecolor=None, edgecolor=None, alpha=0))

    # fill auc
    plt.fill_between(range(1, n_cells), reads_per_cell[1:n_cells], [1] * (n_cells - 1), color='g', alpha=0.2,
                     label='Cell Reads')
    total_barcodes = len(reads_per_cell)
    plt.fill_between(range(n_cells, total_barcodes), reads_per_cell[n_cells:], [1] * (total_barcodes - n_cells),
                     color='grey', alpha=0.2, label='Non-Cell Reads')
    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)

def amplicon_boxplot(cell_df, output_file):
    # generate amplicon box plot showing number of reads per amplicon per cell

    # plot reads per cell for each amplicon in called cells
    plt.figure(figsize=(15, 7))

    df = cell_df.reindex(cell_df.mean().sort_values(ascending=False).index, axis=1)  # sort columns by mean value
    amplicons = list(df.columns)  # amplicon names

    plt.boxplot(np.transpose(df), 0, '', boxprops=dict(facecolor='gainsboro', color='k'), patch_artist=True)

    plt.xticks(np.arange(1, len(amplicons) + 1.2), amplicons, rotation='vertical', fontsize=8)
    plt.yticks(fontsize=12)
    plt.xlim([0, len(amplicons) + 1])
    plt.ylabel('Number of Reads per Cell (N = %d)' % len(cell_df.index), fontsize=14, labelpad=20)
    plt.xlabel('Amplicon', fontsize=14, labelpad=20)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)








