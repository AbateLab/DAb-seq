'''
simple lane merger module
ben demaree 8.2.2019
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

    base_dir = '<base_folder>'

    sample_dirs = os.listdir(base_dir)

    for i in range(len(sample_dirs)):

        # sample_files = [f for f in os.listdir(base_dir + sample_dirs[i]) if f[-9:] == '.fastq.gz']
        samples = list(set([f.split('_')[0] for f in os.listdir(base_dir + sample_dirs[i]) if '.fastq.gz' in f]))
        output_folder = base_dir + sample_dirs[i] + '/merged'

        try:
            os.mkdir(output_folder)
        except OSError:
            pass

        for j in range(len(samples)):

            print 'Merging sample %s' % samples[j]

            sample_fastq = [f for f in os.listdir(base_dir + sample_dirs[i]) if f.split('_')[0] == samples[j]]

            # merge R1
            r1_files = re.compile('.*_R1_.*')
            to_merge_r1 = filter(r1_files.match, sample_fastq)
            to_merge_r1.sort()
            r1_new = re.sub(r'_L\d\d\d_R1_', r'_L000_R1_', to_merge_r1[0])
            r1_new = os.path.join(output_folder, os.path.basename(r1_new))

            to_merge_r1 = [base_dir + sample_dirs[i] + '/' + f for f in to_merge_r1]

            print "     Merging files %s into %s..." % (', '.join(to_merge_r1), r1_new)
            subprocess.call('cat %s > %s' % (' '.join(to_merge_r1), r1_new), shell=True)

            # merge R2 - if R2 files exist
            try:
                r2_files = re.compile('.*_R2_.*')
                to_merge_r2 = filter(r2_files.match, sample_fastq)
                to_merge_r2.sort()
                r2_new = re.sub(r'_L\d\d\d_R2_', r'_L000_R2_', to_merge_r2[0])
                r2_new = os.path.join(output_folder, os.path.basename(r2_new))

                to_merge_r2 = [base_dir + sample_dirs[i] + '/' + f for f in to_merge_r2]

                print "     Merging files %s into %s..." % (', '.join(to_merge_r2), r2_new)
                subprocess.call('cat %s > %s' % (' '.join(to_merge_r2), r2_new), shell=True)

            except IndexError:
                pass


