'''

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 8.8.2019

'''

import subprocess
import os

if __name__ == "__main__":

    # input folder
    input_dir = '/drive4/harish_novaseq/'
    output_dir = '/drive4/harish_novaseq/repaired/'
    fastq_files = [f for f in os.listdir(input_dir) if '.fastq.gz' in f]
    fastq_files.sort()

    print(fastq_files)

    output_files = [output_dir + f.replace('.fastq.gz', '.repaired.fastq.gz') for f in fastq_files]

    for i in range(0, len(output_files), 2):

        in1 = fastq_files[i]
        in2 = fastq_files[i + 1]

        out1 = output_files[i]
        out2 = output_files[i + 1]

        repair_cmd = 'repair.sh in=%s in2=%s out=%s out2=%s' % (in1, in2, out1, out2)
        subprocess.call(repair_cmd, shell=True)