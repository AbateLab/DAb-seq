'''

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 7.9.2019

script for performing joint genotyping of a cohort (including longitudinal data)
the main processing script must be run beforehand
caution: this script uses A LOT of memory

'''

import os
import subprocess
from slackclient import SlackClient
import time

# import functions from external files
import resources_v2

def slack_message(message):
    # for posting a notification to the server-alerts slack channel

    channel = 'server-alerts'
    token = 'xoxp-7171342752-7171794564-486340412737-91fd92781cde6307b077f30f9ea1b700'
    sc = SlackClient(token)

    sc.api_call('chat.postMessage', channel=channel,
                text=message, username='pipelines',
                icon_emoji=':adam:')

def wait(processes):
    # waits for processes to finish
    return [process.communicate() for process in processes]

def db_import(db_path, sample_map_path, interval):
    # import single-cell gvcfs for one interval

    # make call to gatk for genomics db import
    db_import_cmd = 'gatk --java-options "-Xmx4g" GenomicsDBImport ' \
                    '--genomicsdb-workspace-path %s ' \
                    '--batch-size 50 ' \
                    '--reader-threads 2 ' \
                    '--validate-sample-name-map true ' \
                    '-L %s ' \
                    '--sample-name-map %s'% (db_path,
                                             interval,
                                             sample_map_path)

    process = subprocess.Popen(db_import_cmd, shell=True)

    return process

def joint_genotype(db_path, fasta, interval, dbSNP_file, output_vcf):
    # perform joint genotyping across a cohort using data from a genomicsdb store

    # make call to gatk for genotyping
    genotype_cmd = 'gatk --java-options "-Xmx4g" GenotypeGVCFs ' \
                   '-V %s ' \
                   '-R %s ' \
                   '-L %s ' \
                   '-D %s ' \
                   '-O %s ' \
                   '--include-non-variant-sites' \
                   % (db_path,
                      fasta,
                      interval,
                      dbSNP_file,
                      output_vcf)

    process = subprocess.Popen(genotype_cmd, shell=True)

    return process

if __name__ == "__main__":

    ####################################################################################
    # set variables and filepaths for this sample cohort
    ####################################################################################

    # single sample
    output_dir = '/drive3/dabseq/cohorts/aml/abseq21_only/'
    cohort_name = 'abseq21_only'

    # list of sample names
    sample_names = ['abseq21']

    # # pediatric aml
    # output_dir = '/drive3/dabseq/cohorts/aml/abseq8_only_test_all/'
    # cohort_name = 'ped_aml_abseq8'
    #
    # # list of sample names
    # sample_names = ['abseq8',
    #                 'abseq18',
    #                 'abseq10']

    # # patient 655
    # output_dir = '/drive3/dabseq/cohorts/aml/patient_655/'
    # cohort_name = 'patient_655'
    #
    # # list of sample names
    # sample_names = ['abseq11',
    #                 'abseq14',
    #                 'abseq19',
    #                 'abseq21']

    # # patient 577
    # output_dir = '/drive3/dabseq/cohorts/aml/patient_577/'
    # cohort_name = 'patient_577'
    #
    # # list of sample names
    # sample_names = ['abseq13',
    #                 'abseq17',
    #                 'abseq15',
    #                 'abseq16']

    # list of sample ids
    sample_ids = [s.split('seq')[1] for s in sample_names]

    # list of directories containing single-sample gvcfs
    gvcf_paths = ['/drive3/dabseq/output/' + s + '/cells/gvcf/' for s in sample_names]

    # interval file path (EXCLUDING primer coordinates)
    interval_file = '/home/bdemaree/abseq/panel_files/AML/AML.bed'

    # human reference genome fasta file path
    human_fasta = '/drive2/igenomes/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'

    # snpeff config file path
    snpeff_config = '/usr/local/bin/snpEff/snpEff.config'

    # snpeff summary file path
    snpeff_summary = output_dir + cohort_name + '.snpeff_summary.html'

    # genotyped vcf file
    genotyped_vcf = output_dir + cohort_name + '.genotyped.vcf'

    # split vcf file
    split_vcf = output_dir + cohort_name + '.split.vcf'

    # snpeff annotated vcf
    snpeff_annot_vcf = output_dir + cohort_name + '.snpeff.annotated.vcf'

    # final annotated vcf
    annot_vcf = output_dir + cohort_name + '.annotated.vcf'

    # genotype hdf5 file
    geno_hdf5 = output_dir + cohort_name + '.genotypes.hdf5'

    # variant information table
    variants_tsv = output_dir + cohort_name + '.variants.tsv'

    # dbsnp file path
    dbsnp_file = '/drive2/snp_dbs/dbSNP/common_all_20180423.vcf.gz'

    # clinvar vcf file path
    clinvar_vcf = '/drive2/snp_dbs/clinvar_20190805_hg19/clinvar_20190805.hg19.vcf.gz'

    # cosmic vcf file path
    cosmic_vcf = '/drive2/snp_dbs/cosmic_v89_hg19/CosmicAll_v89.hg19.vcf.gz'

    ####################################################################################
    # prepare gvcf and interval files for genotyping
    ####################################################################################

    # send slack notification
    start_time = time.time()
    start_time_fmt = str(time.strftime('%m-%d-%Y %H:%M:%S', time.localtime(start_time)))
    slack_message('Joint genotyping started for cohort %s [%s] at %s.' % (cohort_name,
                                                                                  ', '.join(sample_names),
                                                                                  start_time_fmt))

    # base path for genomics db
    db_dir = output_dir + 'dbs/'
    if not os.path.exists(db_dir):
        os.mkdir(db_dir)

    # base path for single-interval vcfs
    vcf_dir = output_dir + 'vcfs/'
    if not os.path.exists(vcf_dir):
        os.mkdir(vcf_dir)

    # list of lists of gvcf files for each sample
    gvcf_files = []
    for path in gvcf_paths:
        gvcf_files.append([f for f in os.listdir(path) if f.endswith('.g.vcf')])

    # create sample map file
    sample_map_path = output_dir + cohort_name + '.sample_map.tsv'
    with open(sample_map_path, 'w') as f:
        for i in range(len(sample_names)):
            for g in gvcf_files[i]:
                f.write(g.split('.')[0] + '-' + sample_ids[i] + '\t' + gvcf_paths[i] + g + '\n')

    # extract intervals from bed file
    # exclude RUNX1_4 (used for antibodies)
    intervals = {}
    with open(interval_file, 'r') as f:
        for line in f:
            if 'RUNX1_4' not in line:
                fields = line.strip().split('\t')
                intervals[fields[3]] = fields[0] + ':' + fields[1] + '-' + fields[2]

    # create a genomics db and output vcf for each interval
    db_paths = {}
    output_vcfs = {}
    for L in intervals:
        db_paths[L] = db_dir + L + '.genomics.db'
        output_vcfs[L] = vcf_dir + L + '.genotyped.vcf'

    ####################################################################################
    # import gvcfs into genomics DB
    ####################################################################################

    # parallelize by starting a process for each interval
    import_processes = []

    for L in intervals:
        import_processes.append(db_import(db_paths[L],
                                          sample_map_path,
                                          intervals[L]))

    # wait for processes to finish
    wait(import_processes)

    ####################################################################################
    # perform joint genotyping
    ####################################################################################

    # parallelize by starting a process for each interval
    genotype_processes = []

    for L in intervals:
        genotype_processes.append(joint_genotype('gendb://' + db_paths[L],
                                                 human_fasta,
                                                 intervals[L],
                                                 dbsnp_file,
                                                 output_vcfs[L]))

    # wait for processes to finish
    wait(genotype_processes)

    ####################################################################################
    # merge single-interval vcfs
    ####################################################################################

    # call gatk to perform vcf merging
    vcfs_to_merge = ['-I ' + v for v in output_vcfs.values()]
    merge_cmd = 'gatk MergeVcfs %s -O %s' % (' '.join(vcfs_to_merge), genotyped_vcf)

    subprocess.call(merge_cmd, shell=True)

    ####################################################################################
    # split multiallelic sites and annotate vcf
    ####################################################################################

    # split multiallelics, left-align, and trim
    resources_v2.left_align_trim(human_fasta, genotyped_vcf, split_vcf)

    # annotate vcf with snpeff (functional predictions)
    resources_v2.snpeff_annotate(snpeff_summary, snpeff_config, split_vcf, snpeff_annot_vcf)

    # annotate with bcftools
    # use clinvar database
    resources_v2.bcftools_annotate(clinvar_vcf, snpeff_annot_vcf, '-c INFO', annot_vcf)

    # convert vcf to variant matrix in hdf5 format
    resources_v2.vcf_to_tables(annot_vcf, geno_hdf5, variants_tsv)

    print 'Pipeline complete!'

    # send slack notification
    elapsed_time = time.time() - start_time
    elapsed_time_fmt = str(time.strftime('%Hh %Mm %Ss', time.gmtime(elapsed_time)))
    slack_message('Pipeline complete for cohort %s! Total elapsed time is %s.' % (cohort_name, elapsed_time_fmt))