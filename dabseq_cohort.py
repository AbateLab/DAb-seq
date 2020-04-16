'''

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 8.8.2019

a simple script for batch running the dabseq pipeline in parallel

'''

import multiprocessing
import subprocess
import os

def barcode_sample(sample_name, config_file, chem, slack_token):
    # initiates the dabseq pipeline

    barcode_cmd = 'python /home/bdemaree/code/dab-seq/dabseq_pipeline.py ' \
                  '%s barcode %s --sample-name %s --chem %s --slack-token %s' % (cohort_name,
                                                                                 config_file,
                                                                                 sample_name,
                                                                                 chem,
                                                                                 slack_token)
    subprocess.call(barcode_cmd, shell=True)

if __name__ == "__main__":

    # settings for this cohort
    # other flags may be specified in the 'barcode_sample' function

    base_dir = '/drive3/hiv/'
    cohort_name = 'hiv2'
    cohort_dir = base_dir + cohort_name + '/'
    chem = 'V2'
    config_file = '/home/bdemaree/code/dab-seq/cfg/dabseq.hiv2.cfg'
    slack_token_file = '/home/bdemaree/.slack_token'
    with open(slack_token_file, 'r') as f:
        slack_token = f.readline().strip()

    # get sample names from cohort directory
    sample_names = [s for s in os.listdir(cohort_dir) if 'GENOTYPING' not in s]

    # number of samples to process in parallel at a time
    # limit to 2 for dabseq samples on greyhound
    n_parallel = 2

    # cap number of concurrent jobs to total number of samples
    if n_parallel > len(sample_names):
        n_parallel = len(sample_names)

    # create pool of workers and barcode all samples
    pool = multiprocessing.Pool(processes=n_parallel)

    for i in range(len(sample_names)):
        pool.apply_async(barcode_sample, args=(sample_names[i], config_file, chem, slack_token,))

    pool.close()
    pool.join()

    # genotype all samples
    genotype_cmd = 'python /home/bdemaree/code/dab-seq/dabseq_pipeline.py ' \
                   '%s genotype %s ' \
                   '--slack-token %s' % (cohort_name, config_file, slack_token)

    subprocess.call(genotype_cmd, shell=True)