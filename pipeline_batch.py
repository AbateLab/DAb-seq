'''

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 8.8.2019

a simple script for batch running the dabseq pipeline in parallel
'''

import os
import subprocess
import sys
import argparse
import copy
import time
from multiprocessing import Process

def wait(processes):
    # waits for processes to finish
    return [process.communicate() for process in processes]

def run_pipeline(sample_basename, config_file):
    # initiates the mb pipeline

    pipeline_cmd = 'python /home/bdemaree/code/dab-seq/mb_pipeline_v2.py %s %s' % (sample_basename, config_file)

    process = subprocess.Popen(pipeline_cmd, shell=True)

    return process


if __name__ == "__main__":

    sample_basenames = ['abseq' + str(i) for i in range(13, 20)]

    config_file = '/drive3/dabseq/output/dabseq_13-19.cfg'

    # preprocess samples in chunks (hardware-limited)
    samples_per_chunk = 2

    # split cell list into chunks
    sample_chunks = [sample_basenames[i:i + samples_per_chunk]
                     for i in xrange(0, len(sample_basenames), samples_per_chunk)]

    for chunk in sample_chunks:
        # call the dabseq script
        run_pipeline = [run_pipeline(sample, config_file) for sample in chunk]

        # wait for all samples in chunk to finish before continuing
        wait(run_pipeline)