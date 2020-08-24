'''
simple module for merging abs (custom and totalseq)
ben demaree 8.21.2020
'''

import os
import sys
import json
import math
import subprocess
import re
import argparse
import shutil

if __name__ == "__main__":

    input_folder = '/drive3/smith_lab/cohorts/hashing_pilot_1/all_tubes/fastq/ab_fastq/'

    samples = list(set([f.split('totalseq_')[-1].split('cust_')[-1][:5]
                        for f in os.listdir(input_folder) if '.fastq.gz' in f]))

    for j in range(len(samples)):

        print('Merging sample %s' % samples[j])

        sample_fastq = [f for f in os.listdir(input_folder) if samples[j] in f and '.fastq.gz' in f]

        # merge R1
        r1_files = re.compile('.*_R1_.*')
        to_merge_r1 = filter(r1_files.match, sample_fastq)
        to_merge_r1.sort()
        r1_new = re.sub(r'.*_S', input_folder + samples[j] + r'_S', to_merge_r1[0])
        to_merge_r1 = [input_folder + f for f in to_merge_r1]

        print("     Merging files %s into %s..." % (', '.join(to_merge_r1), r1_new))

        subprocess.call('cat %s > %s' % (' '.join(to_merge_r1), r1_new), shell=True)

        # merge R2 - if R2 files exist
        try:
            r2_files = re.compile('.*_R2_.*')
            to_merge_r2 = filter(r2_files.match, sample_fastq)
            to_merge_r2.sort()
            r2_new = re.sub(r'.*_S', input_folder + samples[j] + r'_S', to_merge_r2[0])
            to_merge_r2 = [input_folder + f for f in to_merge_r2]

            print("     Merging files %s into %s..." % (', '.join(to_merge_r2), r2_new))

            subprocess.call('cat %s > %s' % (' '.join(to_merge_r2), r2_new), shell=True)

        except IndexError:
            pass