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

# plotting
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
    sample_name = os.path.basename(ab_clr_count_file).split('.clr.cells.tsv')[0]

    # load antibody hash descriptions
    hashes = resources.load_barcodes(hash_csv, 1, False)
    hash_names = hashes.keys()

    # load clr counts of each ab
    ab_counts_clr = pd.read_csv(ab_clr_count_file, sep='\t', header=0, index_col=0)

    # determine thresholds for each hash and plot fitted distributions
    plt.figure(figsize=(7, 5))
    thresholds = []
    for h in hash_names:

        # fit double gaussian distributions to the hash count data
        # the distributions should fit to the negative and doublet populations
        try:
            data = ab_counts_clr[h]
        except KeyError:
            print('Hashing Ab not in list of sample Abs. Ignoring this hash...')
            continue

        color = plt.get_cmap('tab10')(hash_names.index(h))
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

    # plot all pairs of hashes as scatter plots
    hash_combos = [list(x) for x in itertools.combinations(hash_names, 2)]
    for h in hash_combos:
        plt.figure(figsize=(5, 5))
        sns.scatterplot(x=ab_counts_clr[h[0]], y=ab_counts_clr[h[1]], s=3)
        plt.axvline(x=thresholds[hash_names.index(h[0])])
        plt.plot(np.linspace(-5, 5, 100), [thresholds[hash_names.index(h[1])]] * 100)
        plt.savefig

    # assign cells to samples based on thresholds

