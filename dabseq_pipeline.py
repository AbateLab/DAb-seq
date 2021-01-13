'''

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 2020

the main script for pipeline execution

'''

import os
import shutil
import subprocess
import argparse
import copy
import sys
import time
from multiprocessing.pool import ThreadPool
from multiprocessing import Process

# import functions from external files
import resources
import cell_calling
import hash_demux

def slack_message(message, enabled=False, slack_token=None):
    # for posting a notification to a slack channel

    if enabled:
        client = WebClient(token=slack_token)

        try:
            response = client.chat_postMessage(
                channel='#server-alerts',
                text=message,
                username='pipelines',
                icon_emoji=':adam:')

        except SlackApiError as e:
            # You will get a SlackApiError if "ok" is False
            assert e.response["ok"] is False
            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
            print(f"Got an error: {e.response['error']}")

    else:
        pass

def wait(processes):
    # waits for processes to finish
    return [process.communicate() for process in processes]

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='''
    
    dab-seq: single-cell dna genotyping and antibody sequencing
    ben demaree 2020
    
    input requirements:
    -config file defining file paths and variables (e.g. dabseq.hg19.cfg)
    -raw fastq files (targeted sequencing panel and/or antibody tags)
    -cell and ab barcode csvs
    -panel interval files
    
    requires the following programs in path:
    -gatk
    -bowtie2
    -itdseek (flt3-calling only)
    -samtools
    -bedtools
    -bcftools
    -cutadapt
    -bbmap
    -snpeff
    ''', formatter_class=argparse.RawTextHelpFormatter)

    # required arguments
    parser.add_argument('cohort_name', type=str,
                        help='cohort name')
    parser.add_argument('mode', type=str, choices=['barcode', 'genotype', 'both'],
                        help='pipeline mode ("barcode", "genotype", "both")')
    parser.add_argument('cfg_file', type=str,
                        help='config filename')

    # optional arguments
    parser.add_argument('--sample-name', default=None, type=str,
                        help='sample name (required when mode is "barcoding" or "both")')
    parser.add_argument('--chem', type=str, default='V2', choices=['V1', 'V2'],
                        help='chemistry version (V1 or V2) (default: V2)')
    parser.add_argument('--slack-token', type=str, default=None,
                        help='option to provide slack token string')
    parser.add_argument('--dna-only', action='store_true', default=False,
                        help='option to run dna panel pipeline only')
    parser.add_argument('--dna-panel', type=str, default=None,
                        help='basename (including path) of DNA panel BED files (exclude ".bed" and ".amplicons" '
                             'extensions) - only required for DNA runs')
    parser.add_argument('--ab-only', action='store_true', default=False,
                        help='option to run ab panel pipeline only')
    parser.add_argument('--ab-csv', type=str, default=None,
                        help='ab barcode CSV file (required when processing ab reads)')
    parser.add_argument('--ab-method', type=str, default='all', choices=['unique', 'all'],
                        help='ab clustering method ("unique" or "all") (default: "all")')
    parser.add_argument('--hash-csv', type=str, default=None,
                        help='csv file containing hashing barcodes - if provided, cells will be demultiplexed based on hash')
    parser.add_argument('--r1-dna-min-len', type=int, default=50,
                        help='minimum dna panel read 1 length after barcode and adapter trimming. may need to be'
                             ' modified for shorter sequencing runs (default: 50 bp; >30 bp is recommended')
    parser.add_argument('--r2-dna-min-len', type=int, default=50,
                        help='minimum dna panel read 2 length after adapter trimming. may need to be modified for'
                             ' shorter sequencing runs (default: 50 bp; >30 bp is recommended')
    parser.add_argument('--skip-flt3', action='store_true', default=False,
                        help='option to skip FLT3-ITD calling')
    parser.add_argument('--non-human', action='store_true', default=False,
                        help='option to skip steps requiring a human reference (variant prediction and annotation)')
    parser.add_argument('--ignore-panel-uniformity', action='store_true', default=False,
                        help='option to ignore panel uniformity threshold when calling cells (default cell threshold '
                             'is 60% of amplicons with >= 10X coverage)')
    parser.add_argument('--min-reads', type=int, default=None,
                        help='minimum number of reads for cell calling. when set, will override all other cell calling'
                        ' criteria and only use this read threshold.')
    parser.add_argument('--ploidy', type=int, default=2, choices=[1, 2],
                        help='organism ploidy (default: 2) [haploid (1) or diploid (2) supported only]')

    # parse arguments
    args = parser.parse_args()

    cohort_name = args.cohort_name
    pipeline_mode = args.mode
    cfg_f = args.cfg_file
    sample_name = args.sample_name
    chem = args.chem
    slack_token = args.slack_token
    clustering_method = args.ab_method
    dna_only = args.dna_only
    dna_panel = args.dna_panel
    ab_only = args.ab_only
    ab_barcode_csv = args.ab_csv
    hash_csv = args.hash_csv
    r1_min_len_panel = args.r1_dna_min_len
    r2_min_len_panel = args.r2_dna_min_len
    skip_flt3 = args.skip_flt3
    non_human = args.non_human
    ignore_panel_uniformity = args.ignore_panel_uniformity
    min_reads = args.min_reads
    ploidy = args.ploidy

    # require a sample name when in barcoding mode (or both)
    if ((pipeline_mode == 'barcode' or pipeline_mode == 'both') and sample_name is None) or sample_name == 'GENOTYPING':
        print('Please specify a valid sample name when barcoding.')
        raise SystemExit

    # cohort and sample names cannot include dashes
    if '-' in cohort_name:
        print('Cohort name cannot contain a dash (-).')
        raise SystemExit

    if pipeline_mode == 'barcode' or pipeline_mode == 'both':
        if '-' in sample_name:
            print('Sample name cannot contain a dash (-).')
            raise SystemExit

    # check for config file
    if not os.path.isfile(cfg_f):
        print('Config file not found! Please check the file name and path.')
        raise SystemExit

    # cannot select both dna_only and ab_only
    if dna_only and ab_only:
        print('DNA-only and Ab-only options cannot be used simultaneously.')
        raise SystemExit

    # require ab csv when processing ab reads
    if not dna_only and pipeline_mode != 'genotype':
        if ab_barcode_csv is None:
            print('An antibody CSV barcode file is required when processing Ab reads.')
            raise SystemExit
        elif not os.path.isfile(ab_barcode_csv):
            print('Antibody CSV barcode file not found! Please check the file name and path.')
            raise SystemExit

    # check for ab hashing csv (if specified)
    if hash_csv is not None:
        if not os.path.isfile(hash_csv):
            print('Hashing CSV file not found! Please check the file name and path.')
            raise SystemExit

    # cannot perform genotyping in ab_only mode
    if (pipeline_mode == 'genotyping' or pipeline_mode == 'both') and ab_only:
        print('Cannot perform genotyping using Ab-only setting.')
        raise SystemExit

    # skip flt3 calling when using non-human references
    if non_human:
        skip_flt3 = True

    # ensure min_reads, if set, is a positive number
    if min_reads is not None:
        if min_reads < 1:
            print('Minimum number of reads must be positive.')
            raise SystemExit

    # load config file variables
    # be careful about using exec - never run this pipeline as root
    with open(cfg_f, 'r') as cfg:
        for line in cfg:

            if (line[0] == '#') or (line in ['\n', '\r\n']) or (line[0] == ' '):
                continue

            elif line[0] == '[':

                if '[Barcoding]' in line and pipeline_mode == 'genotype':
                    line = next(cfg)
                    while line[0] != '[':
                        try:
                            line = next(cfg)
                        except StopIteration:
                            break

                if '[DNA]' in line and ab_only:
                    line = next(cfg)
                    while line[0] != '[':
                        try:
                            line = next(cfg)
                        except StopIteration:
                            break

                if '[Antibodies]' in line and (dna_only or pipeline_mode == 'genotype'):
                    line = next(cfg)
                    while line[0] != '[':
                        try:
                            line = next(cfg)
                        except StopIteration:
                            break

                if '[Genotyping]' in line and (ab_only or pipeline_mode == 'barcode'):
                    line = next(cfg)
                    while line[0] != '[':
                        try:
                            line = next(cfg)
                        except StopIteration:
                            break
            else:
                # save a configuration parameter to a variable
                cfg_variable = line.split('=')[0]

                if 'panel_fastq' in cfg_variable and pipeline_mode == 'genotype':
                    pass
                elif 'human' in cfg_variable and non_human:
                    pass
                else:
                    var = line.split("#", 1)[0].strip()  # to remove inline comments
                    exec (var)

    # check all files exist
    all_vars = copy.copy(globals())
    input_files = [all_vars[f] for f in all_vars if '_file' in f and f != '__file__']
    if not ab_only and (pipeline_mode == 'barcode' or pipeline_mode == 'both'):
        input_files.append(bt2_ref + '.1.bt2')
    missing_files = []
    for f in input_files:
        if not os.path.exists(f):
            missing_files.append(f)

    # also check for dna bed files
    if not ab_only:
        # interval file path (EXCLUDING primer coordinates)
        try:
            interval_file = dna_panel + '.bed'
            if not os.path.exists(interval_file):
                missing_files.append(interval_file)
        except TypeError:
            missing_files.append('interval_file')

        # amplicon file path (INCLUDING primer coordinates)
        try:
            amplicon_file = dna_panel + '.amplicons'
            if not os.path.exists(amplicon_file):
                missing_files.append(amplicon_file)
        except TypeError:
            missing_files.append('amplicons_file')

    # print missing files, if any, and exit
    if len(missing_files) > 0:
        print('\nThe following input files could not be found:')
        for f in missing_files:
            print(f)
        print('\nExiting...\n')
        raise SystemExit

    # create base directories for this run
    # if it already exists, ignore and continue
    dirs = [all_vars[d] for d in all_vars if '_dir' in d]
    dirs.sort()
    for d in dirs:
        if not os.path.exists(d):
            os.mkdir(d)

    # option to enable slack messages
    if slack_token is not None:
        slack_enabled = True
    else:
        slack_enabled = False
    if slack_enabled:
        from slack import WebClient
        from slack.errors import SlackApiError

########################################################################################################################
#                                                   BARCODING
########################################################################################################################

    if pipeline_mode == 'barcode' or pipeline_mode == 'both':

        # check that the input fastq directories and create others
        if not ab_only:

            if not os.path.exists(panel_fastq):
                print('FASTQ input directory (%s) does not exist! Exiting...\n' % panel_fastq)
                raise SystemExit

            elif os.listdir(panel_fastq) == []:
                print('FASTQ input directory (%s) is empty! Exiting...\n' % panel_fastq)
                raise SystemExit

            # the following files and directories will be automatically created

            by_cell_dir = sample_dir + 'cells/'             # directory for single-cell files
            by_cell_fastq_dir = by_cell_dir + 'fastq/'      # directory for fastq by cell
            by_cell_bam_dir = by_cell_dir + 'bam/'          # directory for bam (and bai) by cell
            by_cell_gvcf_dir = by_cell_dir + 'gvcf/'        # directory for gvcf by cell
            by_cell_flt3_dir = by_cell_dir + 'flt3/'        # directory for flt3 itd calling by cell

        if not dna_only:

            if not os.path.exists(ab_fastq):
                print('FASTQ input directory (%s) does not exist! Exiting...\n' % ab_fastq)
                raise SystemExit

            elif os.listdir(ab_fastq) == []:
                print('FASTQ input directory (%s) is empty! Exiting...\n' % ab_fastq)
                raise SystemExit

            # the following files and directories will be automatically created
            ab_dir = sample_dir + 'abs/'  # directory for storing antibody data
            umi_counts_merged = ab_dir + sample_name + '.umi_counts.tsv'  # sample umi counts file

        # create other directories for this run
        # if it already exists, ignore and continue
        all_vars = copy.copy(globals())
        dirs = [all_vars[d] for d in all_vars if '_dir' in d]
        dirs.sort()
        for d in dirs:
            if not os.path.exists(d):
                os.mkdir(d)

        # setup barcode extraction parameters based on chemistry
        if chem == 'V1':
            cell_barcode_csv_file = sys.path[0] + '/barcodes/mb_cell_barcodes_v1.csv'

            r1_start = 'CGATGACG'
            r1_end = 'CTGTCTCTTATACACATCT'
            r2_end = 'CGTCATCG'

            bar_ind_1, bar_ind_2 = list(range(8)), list(range(-8, 0))

        elif chem == 'V2':
            cell_barcode_csv_file = sys.path[0] + '/barcodes/mb_cell_barcodes_v2.csv'

            r1_start = 'GTACTCGCAGTAGTC'
            r1_end = 'CTGTCTCTTATACACATCT'
            r2_end = 'GACTACTGCGAGTAC'

            bar_ind_1, bar_ind_2 = list(range(9)), list(range(-9, 0))

        # set minimum ab sequence lengths after trimming
        # should be acceptable for most tags
        r1_min_len_ab = 0
        r2_min_len_ab = 20

        # send slack notification
        start_time = time.time()
        start_time_fmt = str(time.strftime('%m-%d-%Y %H:%M:%S', time.localtime(start_time)))
        slack_message('Barcoding pipeline started for cohort %s, sample %s at %s.' % (cohort_name,
                                                                                      sample_name,
                                                                                      start_time_fmt),
                      slack_enabled,
                      slack_token)

        # get input file names and store in TapestriTube objects
        R1_files, R2_files = [], []

        if not ab_only:
            # get fastq files for dna panel
            R1_files_panel, R2_files_panel = resources.get_fastq_names(panel_fastq)

            R1_files += R1_files_panel
            R2_files += R2_files_panel

        if not dna_only:

            # get fastq files for antibodies
            R1_files_ab, R2_files_ab = resources.get_fastq_names(ab_fastq)

            R1_files += R1_files_ab
            R2_files += R2_files_ab

        R1_files.sort()
        R2_files.sort()

        # store sample info in TapestriTube objects
        if dna_only:
            ab_dir = ''
        tubes = resources.generate_sample(R1_files, R2_files, dna_only, ab_only, sample_name, temp_dir, ab_dir)

        # display and write sample summary file
        resources.file_summary(tubes, cfg_f, summary, header_file=sys.path[0]+'/squidward.txt')

        print('''
####################################################################################
# filter reads for cell barcode and perform error correction
####################################################################################
''')

        # load mission bio barcode csv file
        barcodes = resources.load_barcodes(cell_barcode_csv_file, 1, False)

        # generate hamming dictionary for error correction
        barcodes = resources.generate_hamming_dict(barcodes)

        print('%d unique barcode sequences loaded into dictionary.\n' % len(barcodes))

        # for panel reads, filter reads with valid barcode structure and export to new fastq
        print('Extracting barcodes from raw fastq files...\n')

        # process reads with valid barcode
        process_barcodes = []

        for tube in tubes:
            if not ab_only:
                # panel files
                p = Process(
                    target=tube.barcode_reads,
                    args=(r1_start,
                          r1_end,
                          r2_end,
                          barcodes,
                          bar_ind_1,
                          bar_ind_2,
                          r1_min_len_panel,
                          r2_min_len_panel,
                          'panel'))
                process_barcodes.append(p)
                p.start()

            if not dna_only:
                # ab files
                p = Process(
                    target=tube.barcode_reads,
                    args=(r1_start,
                          r1_end,
                          r2_end,
                          barcodes,
                          bar_ind_1,
                          bar_ind_2,
                          r1_min_len_ab,
                          r2_min_len_ab,
                          'ab'))
                process_barcodes.append(p)
                p.start()

        # wait for processes to finish
        for p in process_barcodes:
            p.join()

        if not dna_only:
            print('''
####################################################################################
# import ab reads, error correct, and count
####################################################################################
''')

            # load ab barcode csv file (with descriptions)
            ab_barcodes = resources.load_barcodes(ab_barcode_csv, 1, False)
            barcode_descriptions = copy.deepcopy(ab_barcodes)

            # generate hamming dictionary for error correction
            ab_barcodes = resources.generate_hamming_dict(ab_barcodes)

            # process ab reads and look for barcode
            ab_process = []
            for tube in tubes:
                # extract ab reads from files
                p = Process(
                    target=tube.process_abs,
                    args=(ab_barcodes,
                          barcode_descriptions,
                          ab_handles,
                          ab_bar_coord,
                          ab_umi_coord,
                          min_umi_qual))
                ab_process.append(p)
                p.start()

            # wait for processes to finish
            for p in ab_process:
                p.join()

            # sort ab reads files
            sort_files = []
            for tube in tubes:
                p = subprocess.Popen('sort --parallel=2 %s -o %s' % (tube.ab_reads, tube.ab_reads), shell=True)
                sort_files.append(p)

            # wait for processes to finish
            wait(sort_files)

            # count the ab read umis
            ab_counts = []
            for tube in tubes:
                p = Process(
                    target=tube.count_umis,
                    args=(clustering_method, ))
                ab_counts.append(p)
                p.start()

            # wait for processes to finish
            for p in ab_counts:
                p.join()

            # merge umi counts files
            subprocess.call('head -1 %s > %s' % (tubes[0].umi_counts,
                                                 umi_counts_merged), shell=True)

            subprocess.call('tail -n +2 -q %s >> %s' % (' '.join([t.umi_counts for t in tubes]),
                                                        umi_counts_merged), shell=True)

            # remove old umi count files
            umi_file_cleanup = [os.remove(t.umi_counts) for t in tubes]

            # create table of ab counts for all barcodes
            resources.umi_counts_by_cell(umi_counts_merged,
                                         ab_barcode_csv,
                                         ab_dir,
                                         None)

        if not ab_only:
            print('''
###################################################################################
# count read alignments to inserts
###################################################################################
''')

            # get r1 reads for all panel tubes
            panel_r1_files = [t.panel_r1_temp for t in tubes]

            # align r1 reads to inserts to obtain read counts across all barcodes
            all_tsv = barcode_dir + sample_name + '.all.tsv'  # alignment counts for all barcodes
            aln_stats_file = barcode_dir + sample_name + '.alignment_stats.txt'
            resources.count_alignments(panel_r1_files,
                                       amplicon_file,
                                       ref_fasta_file,
                                       all_tsv,
                                       aln_stats_file,
                                       temp_dir)

            print('''
###################################################################################
# call valid cells using selected method
###################################################################################
''')

            # call valid cells using cell_caller function
            valid_cells = cell_calling.call(barcode_dir,
                                            sample_name,
                                            min_reads,
                                            'second_derivative',
                                            threshold=(not ignore_panel_uniformity))

            # create SingleCell objects for each valid cell
            cells = [resources.SingleCell(barcode,
                                          by_cell_fastq_dir,
                                          by_cell_bam_dir,
                                          by_cell_gvcf_dir,
                                          by_cell_flt3_dir)
                     for barcode in valid_cells]

            print('%s valid cells found!' % len(cells))

            if not dna_only:
                # create table of per-cell ab counts for valid cells only
                # also calculate clr for valid cells
                ab_counts_unique = resources.umi_counts_by_cell(umi_counts_merged,
                                                                ab_barcode_csv,
                                                                ab_dir,
                                                                cells)

                # if cell hashing csv is specified, demultiplex cells
                if hash_csv is not None:
                    hash_demux.demux_cells(ab_counts_unique,
                                           hash_csv,
                                           ab_dir + 'hashes/',
                                           clr_abs='hashing_only')

            print('''
###################################################################################
# write valid cells from panel reads to separate fastq files
###################################################################################
''')

            # check that all barcodes have the same length for demultiplexing
            bar_length = list(set([len(k) for k in valid_cells]))
            assert len(bar_length) == 1, 'All barcodes must have the same length!'

            # split cell fastq files for panel reads associated with valid cells
            barcode_names = temp_dir + 'barcodes.txt'
            with open(barcode_names, 'w') as f:
                for b in valid_cells:
                    f.write(b + '\n')

            # split files by cell barcode using bbmap demuxbyname.sh
            split_files = []
            for tube in tubes:
                demux_cmd = 'demuxbyname.sh prefixmode=f -Xmx10g length=%d in1=%s in2=%s out=%s names=%s' % (
                    bar_length[0],
                    tube.panel_r1_temp,
                    tube.panel_r2_temp,
                    by_cell_fastq_dir + '%.fastq',
                    barcode_names
                )
                p = subprocess.Popen(demux_cmd, shell=True)
                split_files.append(p)

            # wait for processes to finish
            wait(split_files)

            print('''
####################################################################################
# align, convert, sort, index panel reads (optional: call FLT3-ITDs)
####################################################################################
''')

            # limit number of cells to preprocess at a time (based on hardware limitations)
            n_preprocess = 100  #100/300 standard/greyhound
            if n_preprocess > len(cells):
                n_preprocess = len(cells)

            # create pool of workers and run through all cells
            preprocess_pool = ThreadPool(processes=n_preprocess)

            # align and index cells
            for c in cells:
                preprocess_pool.apply_async(resources.SingleCell.align_and_index, args=(c, bt2_ref,))

            preprocess_pool.close()
            preprocess_pool.join()

            # optionally, call FLT3-ITDs using ITDseek
            if not skip_flt3:

                # verify flt3 folder is empty before continuing (remove and recreate folder)
                shutil.rmtree(by_cell_flt3_dir)
                os.mkdir(by_cell_flt3_dir)

                preprocess_pool = ThreadPool(processes=n_preprocess)

                for c in cells:
                    preprocess_pool.apply_async(resources.SingleCell.call_flt3, args=(c, ref_fasta_file,))

                preprocess_pool.close()
                preprocess_pool.join()

            else:
                os.rmdir(by_cell_flt3_dir)

            print('''
####################################################################################
# perform variant calling for all cells
####################################################################################
''')
            # delete by-cell fastq directory to save space
            assert 2 * len(os.listdir(by_cell_fastq_dir)) == len(os.listdir(by_cell_bam_dir)),\
                'Number of single-cell files in FASTQ directory is inconsistent with the number in the BAM directory.'
            shutil.rmtree(by_cell_fastq_dir)

            # verify gvcf folder is empty before continuing (remove and recreate folder)
            shutil.rmtree(by_cell_gvcf_dir)
            os.mkdir(by_cell_gvcf_dir)

            # limit number of cells to call variants at a time (based on hardware limitations)
            n_call_variants = 30   #30/70 standard/greyhound
            if n_call_variants > len(cells):
                n_call_variants = len(cells)

            # create pool of workers and run through all cells
            call_variants_pool = ThreadPool(processes=n_call_variants)

            for c in cells:
                call_variants_pool.apply_async(resources.SingleCell.call_variants, args=(c,
                                                                                         ref_fasta_file,
                                                                                         interval_file,
                                                                                         ploidy,))

            call_variants_pool.close()
            call_variants_pool.join()

        # send slack notification
        elapsed_time = time.time() - start_time
        elapsed_time_fmt = str(time.strftime('%Hh %Mm %Ss', time.gmtime(elapsed_time)))
        slack_message('Barcoding complete for cohort %s, sample %s! Total elapsed time is %s.' % (cohort_name,
                                                                                                  sample_name,
                                                                                                  elapsed_time_fmt),
                      slack_enabled,
                      slack_token)

########################################################################################################################
#                                                   GENOTYPING
########################################################################################################################

    if pipeline_mode == 'genotype' or pipeline_mode == 'both':

        print('''
####################################################################################
# prepare gvcf and interval files for genotyping
####################################################################################
''')
        # get list of all samples in this cohort
        sample_names = [f for f in os.listdir(cohort_dir) if f != 'GENOTYPING']

        # the following files and directories will be automatically created

        sample_map = cohort_genotyping_dir + cohort_name + '.sample_map.tsv'                # sample map for gatk
        genotpying_summary = cohort_genotyping_dir + cohort_name + '.cell_counts.tsv'       # number of cells per sample
        snpeff_summary = cohort_genotyping_dir + cohort_name + '.snpeff_summary.html'       # snpeff summary file path
        genotyped_vcf = cohort_genotyping_dir + cohort_name + '.genotyped.vcf'              # genotyped vcf file
        split_vcf = cohort_genotyping_dir + cohort_name + '.split.vcf'                      # split vcf file
        snpeff_annot_vcf = cohort_genotyping_dir + cohort_name + '.snpeff.annotated.vcf'    # snpeff annotated vcf
        annot_vcf = cohort_genotyping_dir + cohort_name + '.annotated.vcf'                  # final annotated vcf
        geno_hdf5 = cohort_genotyping_dir + cohort_name + '.genotypes.hdf5'                 # genotype hdf5 file
        variants_tsv = cohort_genotyping_dir + cohort_name + '.variants.tsv'                # variant information table
        flt3_vcf = cohort_genotyping_dir + cohort_name + '.flt3.vcf'                        # combined flt3 vcf

        # find single-cell files for genotyping
        # check that barcodes match in single-cell files and cell tsv for each sample

        gvcf_files = []
        itd_files = []
        num_cells = []

        for s in sample_names:

            # get tsv barcodes
            s_cell_tsv = cohort_dir + '%s/barcodes/%s.cells.tsv' % (s, s)
            tsv_barcodes = []
            first_line = True
            with open(s_cell_tsv, 'r') as c:
                for line in c:
                    if first_line:
                        first_line = False
                    else:
                        tsv_barcodes.append(line.split('\t')[0].encode('utf8'))

            # get and check gvcf files
            gvcf_folder = cohort_dir + s + '/cells/gvcf/'
            s_gvcf = [f[:-6].encode('utf8') for f in os.listdir(gvcf_folder) if f.endswith('.g.vcf')]

            assert sorted(s_gvcf) == sorted(tsv_barcodes), 'GVCF and TSV barcodes do not match!'

            gvcf_files.append([gvcf_folder + b.decode('utf-8') + '.g.vcf' for b in s_gvcf])
            num_cells.append(len(tsv_barcodes))

            if not skip_flt3:

                # get and check flt3-itd files
                itd_folder = cohort_dir + s + '/cells/flt3/'
                s_itd = [f[:-12].encode('utf8') for f in os.listdir(itd_folder) if f.endswith('.vcf')]

                assert sorted(s_itd) == sorted(tsv_barcodes), 'FLT3 and TSV barcodes do not match!'

                itd_files.append([itd_folder + b.decode('utf-8') + '.flt3itd.vcf' for b in s_itd])

        # base path for genomics db
        db_dir = cohort_genotyping_dir + 'dbs/'
        if os.path.exists(db_dir):
            print('The genomics DB directory %s already exists. Please delete or rename it.\n' % db_dir)
            raise SystemExit
        os.mkdir(db_dir)

        # base path for single-interval vcfs
        vcf_dir = cohort_genotyping_dir + 'vcfs/'
        if not os.path.exists(vcf_dir):
            os.mkdir(vcf_dir)

        # write a short cell count file
        with open(genotpying_summary, 'w') as f:
            f.write('sample\tnumber_of_cells\n')
            for i in range(len(sample_names)):
                f.write(sample_names[i] + '\t' + str(num_cells[i]) + '\n')

        # create sample map file
        with open(sample_map, 'w') as f:
            for i in range(len(sample_names)):
                for g in gvcf_files[i]:
                    f.write(os.path.basename(g).split('.')[0] + '-' + sample_names[i] + '\t' + g + '\n')

        # extract intervals from bed file
        intervals = {}
        with open(interval_file, 'r') as f:
            for line in f:
                if exclude_RUNX1_4 and 'RUNX1_4' in line:
                    continue
                else:
                    fields = line.strip().split('\t')
                    intervals[fields[3]] = fields[0] + ':' + fields[1] + '-' + fields[2]

        # create a genomics db and output vcf for each interval
        db_paths = {}
        output_vcfs = {}
        for L in intervals:
            db_paths[L] = db_dir + L + '.genomics.db'
            output_vcfs[L] = vcf_dir + L + '.genotyped.vcf'

        # send slack notification
        start_time = time.time()
        start_time_fmt = str(time.strftime('%m-%d-%Y %H:%M:%S', time.localtime(start_time)))
        slack_message('Genotyping pipeline started for cohort %s at %s.' % (cohort_name, start_time_fmt),
                      slack_enabled,
                      slack_token)

        print('''
####################################################################################
# import gvcfs into genomics DB
####################################################################################
''')

        # limit number of intervals to import at a time (based on hardware limitations)
        n_import = 30       #30/60 standard/greyhound
        if n_import > len(intervals):
            n_import = len(intervals)

        # create pool of workers and run through all samples
        import_pool = ThreadPool(processes=n_import)

        # import intervals
        for L in intervals:
            import_pool.apply_async(resources.db_import, args=(db_paths[L],
                                                               sample_map,
                                                               intervals[L],))

        import_pool.close()
        import_pool.join()

        print('''
####################################################################################
# perform joint genotyping
####################################################################################
''')

        # limit number of intervals to import at a time (based on hardware limitations)
        n_genotype = 30     #30/60 standard/greyhound
        if n_genotype > len(intervals):
            n_genotype = len(intervals)

        # create pool of workers and run through all samples
        genotype_pool = ThreadPool(processes=n_genotype)

        # joint genotype intervals
        for L in intervals:
            genotype_pool.apply_async(resources.joint_genotype, args=('gendb://' + db_paths[L],
                                                                      ref_fasta_file,
                                                                      intervals[L],
                                                                      output_vcfs[L],
                                                                      ploidy,))

        genotype_pool.close()
        genotype_pool.join()

        print('''
####################################################################################
# merge single-interval vcfs
####################################################################################
''')

        # call gatk to perform vcf merging
        vcfs_to_merge = ['-I ' + v for v in output_vcfs.values()]
        merge_cmd = 'gatk MergeVcfs %s -O %s' % (' '.join(vcfs_to_merge), genotyped_vcf)

        subprocess.call(merge_cmd, shell=True)

        print('''
####################################################################################
# split multiallelic sites and annotate vcf (human reference only)
####################################################################################
''')

        # split multiallelics, left-align, and trim
        resources.left_align_trim(ref_fasta_file, genotyped_vcf, split_vcf)

        if not non_human:

            # annotate vcf with snpeff (functional predictions)
            resources.snpeff_annotate(snpeff_summary, snpeff_human_config_file, split_vcf, snpeff_annot_vcf)

            # annotate with bcftools
            # use clinvar database
            resources.bcftools_annotate(clinvar_human_vcf_file, snpeff_annot_vcf, '-c INFO', annot_vcf)

        print('''
####################################################################################
# export genotyping calls to compressed HDF5 file
####################################################################################
''')

        # convert vcf to variant matrix in hdf5 format

        # combine flt3 itd vcf files from all samples
        # only include variants with sufficient read depth
        if not skip_flt3 and not non_human:
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
                                cell_barcode = os.path.basename(itd).split('.')[0] + '-' + sample_names[i]

                                if line[0] != '#':
                                    depth = int(line.split('\t')[5])
                                    vaf = float(line.strip().split('=')[-1])
                                    if depth >= min_itd_dp and vaf >= min_itd_vaf:
                                        # add cell barcode to id column
                                        vcf_record = '\t'.join(line.split('\t')[:2] + [cell_barcode] + line.split('\t')[3:])
                                        v.write(vcf_record)
                                    else:
                                        break

            resources.vcf_to_tables(annot_vcf, geno_hdf5, variants_tsv, ploidy, flt3_vcf, non_human=False)

        # no flt3 itd calling with human samples
        elif not non_human:
            resources.vcf_to_tables(annot_vcf, geno_hdf5, variants_tsv, ploidy, itd_vcf_file=False, non_human=False)

        # non-human samples
        elif non_human:
            resources.vcf_to_tables(split_vcf, geno_hdf5, variants_tsv, ploidy, itd_vcf_file=False, non_human=True)

        # add amplicon count tables (from valid cells) to hdf5 file
        amplicon_counts_tsvs = [cohort_dir + s + '/barcodes/' + s + '.cells.tsv' for s in sample_names]
        resources.add_amplicon_counts(geno_hdf5, sample_names, amplicon_counts_tsvs)

        # add hashing tables, if available, to hdf5 file
        hash_tables = [cohort_dir + s + '/abs/hashes/' + s + '.sample_hashes.tsv' for s in sample_names]
        hash_tables = [f for f in hash_tables if os.path.isfile(f)]
        if len(hash_tables) > 0:
            resources.add_hashes(geno_hdf5, hash_tables)

        # if abs in experiment, add counts to hdf5 file
        # also perform antibody correction using GLM
        # note: clr is calculated independently for each cell, whereas the corrected counts regress across all cells
        if not dna_only:
            ab_count_dirs = [cohort_dir + s + '/abs/by_method.CELLS/' for s in sample_names]
            resources.add_hdf5_ab_counts(geno_hdf5, sample_names, ab_count_dirs)
            resources.GLM_regression(geno_hdf5, ['IgG1', 'amplicon_total', 'ab_total_raw', 'ab_total_umi'])

        # compress the final annotated vcf to reduce size
        if non_human:
            subprocess.call('pigz -p 8 -f %s' % split_vcf, shell=True)
        else:
            subprocess.call('pigz -p 8 -f %s' % annot_vcf, shell=True)

        # send slack notification
        elapsed_time = time.time() - start_time
        elapsed_time_fmt = str(time.strftime('%Hh %Mm %Ss', time.gmtime(elapsed_time)))
        slack_message('Genotyping pipeline complete for cohort %s! Total elapsed time is %s.' % (cohort_name,
                                                                                                 elapsed_time_fmt),
                      slack_enabled,
                      slack_token)

########################################################################################################################

    # delete temporary files, if selected
    if file_cleanup:

        # remove temp folder from barcoding mode
        if pipeline_mode == 'barcode' or pipeline_mode == 'both':
            try:
                shutil.rmtree(temp_dir)
            except OSError:
                pass

        # remove temp files from genotyping mode
        if pipeline_mode == 'genotype' or pipeline_mode == 'both':

            # remove limited files for non-human samples
            if non_human:
                try:
                    os.remove(sample_map)
                    os.remove(genotyped_vcf)
                    os.remove(genotyped_vcf + '.idx')

                except OSError:
                    pass

            # human samples (additional files to remove)
            else:
                try:
                    os.remove(sample_map)
                    os.remove(genotyped_vcf)
                    os.remove(genotyped_vcf + '.idx')
                    os.remove(split_vcf)
                    os.remove(snpeff_annot_vcf + '.gz')
                    os.remove(snpeff_annot_vcf + '.gz.tbi')

                except OSError:
                    pass

            # remove single-interval vcfs and dbs
            try:
                shutil.rmtree(db_dir)
                shutil.rmtree(vcf_dir)
            except OSError:
                pass

    print('Pipeline complete!\n')