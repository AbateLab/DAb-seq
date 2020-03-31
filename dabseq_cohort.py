'''

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 8.8.2019

a simple script for batch running the dabseq pipeline in parallel

'''

import multiprocessing
import subprocess
import os

def barcode_sample(sample_name, config_file, chem):
    # initiates the dabseq pipeline

    barcode_cmd = 'python /home/bdemaree/code/dab-seq/dabseq_pipeline.py ' \
                  '%s barcode %s --sample-name %s --chem %s' % (cohort_name, config_file, sample_name, chem)
    subprocess.call(barcode_cmd, shell=True)

if __name__ == "__main__":

    # settings for this cohort
    # other flags may be specified in the 'barcode_sample' function

    cohort_name = 'cohort1'
    cohort_dir = '/drive3/testing/' + cohort_name + '/'
    chem = 'V1'
    config_file = '/home/bdemaree/code/dab-seq/cfg/dabseq.cfg'

    # get sample names from cohort directory
    sample_names = [s for s in os.listdir(cohort_dir) if 'GENOTYPING' not in s]

    # number of samples to process in parallel at a time
    # limit to 2 for dabseq samples on greyhound
    n_parallel = 2

    # cap number of concurrent jobs to total number of samples
    if n_parallel > len(sample_names):
        n_parallel = len(sample_names)

    # # create pool of workers and barcode all samples
    # pool = multiprocessing.Pool(processes=n_parallel)
    #
    # for i in range(len(sample_names)):
    #     pool.apply_async(barcode_sample, args=(sample_names[i], config_file, chem,))
    #
    # pool.close()
    # pool.join()

    # genotype all samples
    genotype_cmd = 'python /home/bdemaree/code/dab-seq/dabseq_pipeline.py ' \
                   '%s genotype %s' % (cohort_name, config_file)

    subprocess.call(genotype_cmd, shell=True)