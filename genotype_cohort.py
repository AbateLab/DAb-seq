'''

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 7.9.2019

script for performing joint genotyping of a cohort (multiple samples and/or longitudinal data)
the main processing script must be run beforehand
caution: this script uses A LOT of memory - run one or two samples at a time

'''

import os
import subprocess
import time
import sys
from collections import OrderedDict

# import functions from external files
import resources

# option to enable slack messages
slack_enabled = True
slack_token_file = '/home/bdemaree/.slack_token'

def slack_message(message, enabled=False, slack_token_file=None):
    # for posting a notification to a slack channel

    if enabled:
        from slackclient import SlackClient

        with open(slack_token_file) as f:
            token = f.readline().strip()

        channel = 'server-alerts'
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

def joint_genotype(db_path, fasta, interval, output_vcf):
    # perform joint genotyping across a cohort using data from a genomicsdb store

    # make call to gatk for genotyping
    genotype_cmd = 'gatk --java-options "-Xmx4g" GenotypeGVCFs ' \
                   '-V %s ' \
                   '-R %s ' \
                   '-L %s ' \
                   '-O %s ' \
                   '--include-non-variant-sites' \
                   % (db_path,
                      fasta,
                      interval,
                      output_vcf)

    process = subprocess.Popen(genotype_cmd, shell=True)

    return process

if __name__ == "__main__":

    ####################################################################################
    # set variables and filepaths for sample cohorts
    ####################################################################################

    # dict of cohort names and corresponding abseq samples

    cohorts = OrderedDict()
    cohorts['abseq6'] = ['abseq6']
    cohorts['abseq9'] = ['abseq9']
    cohorts['ped_aml'] = ['abseq8', 'abseq18', 'abseq10']
    cohorts['patient_577'] = ['abseq13', 'abseq17', 'abseq15', 'abseq16']
    cohorts['patient_655'] = ['abseq11', 'abseq14', 'abseq19', 'abseq21']

    base_output_dir = '/drive3/dabseq/cohorts/aml/'

    for cohort_name in cohorts:

        sample_names = cohorts[cohort_name]
        output_dir = base_output_dir + cohort_name + '/'

        # create cohort folder if it does not exist
        # note that GATK will not overwrite existing genomicsdbs
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # option to analyze flt3-itd variants
        flt3_itd = True

        # minimum read depth and single-cell vaf for all flt3 variants
        # this helps filter false positives
        min_itd_dp = 15
        min_itd_vaf = 0.15

        # list of sample ids
        sample_ids = [s.split('seq')[1] for s in sample_names]

        # sample map for gatk
        sample_map = output_dir + cohort_name + '.sample_map.tsv'

        # list of directories containing single-sample gvcfs
        gvcf_paths = ['/drive3/dabseq/output/' + s + '/cells/gvcf/' for s in sample_names]

        # list of directories containing flt3-itd vcf files
        itd_paths = ['/drive3/dabseq/output/' + s + '/cells/flt3/' for s in sample_names]

        # interval file path (EXCLUDING primer coordinates)
        interval_file = sys.path[0] + '/panel_files/AML/AML.bed'

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

        # combined flt3 vcf
        flt3_vcf = output_dir + cohort_name + '.flt3.vcf'

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
                                                                              start_time_fmt),
                      slack_enabled,
                      slack_token_file)

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

        # optional: list of lists of flt3 vcf files for each sample
        if flt3_itd:
            itd_files = []
            for path in itd_paths:
                itd_files.append([path + f for f in os.listdir(path) if f.endswith('.vcf')])

        # create sample map file
        with open(sample_map, 'w') as f:
            for i in range(len(sample_names)):
                for g in gvcf_files[i]:
                    f.write(g.split('.')[0] + '-' + sample_ids[i] + '\t' + gvcf_paths[i] + g + '\n')

        # extract intervals from bed file
        # option to exclude RUNX1_4 (used for antibodies)
        exclude_RUNX1_4 = True
        intervals = {}
        with open(interval_file, 'r') as f:
            for line in f:
                if exclude_RUNX1_4 and 'RUNX1_4' not in line:
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
                                              sample_map,
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
        resources.left_align_trim(human_fasta, genotyped_vcf, split_vcf)

        # annotate vcf with snpeff (functional predictions)
        resources.snpeff_annotate(snpeff_summary, snpeff_config, split_vcf, snpeff_annot_vcf)

        # annotate with bcftools
        # use clinvar database
        resources.bcftools_annotate(clinvar_vcf, snpeff_annot_vcf, '-c INFO', annot_vcf)

        # convert vcf to variant matrix in hdf5 format
        # optional: include itd calls

        if flt3_itd:
            # combine flt3 itd vcf files from all samples
            # only include variants with sufficient read depth

            with open(flt3_vcf, 'w') as v:
                header = True
                for i in range(len(sample_names)):
                    for itd in itd_files[i]:
                        with open(itd, 'r') as f:
                            for line in f:
                                if header and line[:2] == '##':
                                    v.write(line)
                                    continue

                                elif header and line[:2] == '#C':
                                    v.write(line)

                                header = False
                                cell_barcode = itd.split('.')[0].split('/')[-1] + '-' + sample_ids[i]

                                if line[0] != '#':
                                    depth = int(line.split('\t')[5])
                                    vaf = float(line.strip().split('=')[-1])
                                    if depth >= min_itd_dp and vaf >= min_itd_vaf:
                                        # add cell barcode to id column
                                        vcf_record = '\t'.join(line.split('\t')[:2] + [cell_barcode] + line.split('\t')[3:])
                                        v.write(vcf_record)
                                    else:
                                        break

            resources.vcf_to_tables(annot_vcf, geno_hdf5, variants_tsv, flt3_vcf)

        else:
            resources.vcf_to_tables(annot_vcf, geno_hdf5, variants_tsv)

        # compress the final annotated vcf to reduce size
        subprocess.call('gzip %s' % annot_vcf, shell=True)

        ####################################################################################
        # clean up temporary files
        ####################################################################################

        try:
            os.remove(sample_map)
            os.remove(genotyped_vcf)
            os.remove(genotyped_vcf + '.idx')
            os.remove(split_vcf)
            os.remove(snpeff_annot_vcf + '.gz')
            os.remove(snpeff_annot_vcf + '.gz.tbi')

        except OSError:
            pass

        ####################################################################################

        print 'Pipeline complete!'

        # send slack notification
        elapsed_time = time.time() - start_time
        elapsed_time_fmt = str(time.strftime('%Hh %Mm %Ss', time.gmtime(elapsed_time)))
        slack_message('Pipeline complete for cohort %s! Total elapsed time is %s.' % (cohort_name, elapsed_time_fmt),
                      slack_enabled,
                      slack_token_file)
