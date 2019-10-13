'''

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 8.8.2019

a simple script for batch running the dabseq pipeline in parallel

'''

import multiprocessing
import subprocess

def pipeline(sample_basename, config_file):
    # initiates the mb pipeline

    pipeline_cmd = 'python /home/bdemaree/code/dab-seq/mb_pipeline_v2.py --gvcf-only %s %s' % (sample_basename, config_file)
    subprocess.call(pipeline_cmd, shell=True)

if __name__ == "__main__":

    config_file = '/home/bdemaree/code/dab-seq/dabseq.cfg'

    # sample_basenames = ['abseq' + str(i) for i in [6, 8, 9, 10, 11] + range(13, 20)]
    sample_basenames = ['abseq' + str(i) for i in [11, 14, 19, 20, 21, 8, 18, 10, 13, 17, 15, 16, 6, 9]]

    # number of samples to process in parallel at a time
    # limit to 2 for dabseq samples on greyhound
    n_parallel = 2

    # cap number of concurrent jobs to total number of samples
    if n_parallel > len(sample_basenames):
        n_parallel = len(sample_basenames)

    # create pool of workers and run through all samples
    pool = multiprocessing.Pool(processes=n_parallel)

    for i in range(len(sample_basenames)):
        pool.apply_async(pipeline, args=(sample_basenames[i], config_file,))

    pool.close()
    pool.join()