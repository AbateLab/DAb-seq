'''

dab-seq: single-cell abseq + genotyping
ben demaree 7.9.2019

input requirements:
-raw fastq files (panel and abs)
-cell and ab barcode csvs
-panel bed file

requires the following programs in path:
-gatk
-bowtie2
-samtools
-bedtools
-cutadapt
-bbmap

umis counted using methods described in Genome Research paper by Smith et al:
https://genome.cshlp.org/content/early/2017/01/18/gr.209601.116

'''

import os
import subprocess
import sys
import argparse
import copy
from slackclient import SlackClient
import time
from multiprocessing import Process

# add umi-tools folder to path
sys.path.insert(0, 'umi-tools')
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

def file_summary(samples, filename='file_summary.txt', to_file=True):
    # displays sample files and class variables

    if to_file:
        f = open(filename, 'w')

    for sample in samples:

        s = vars(sample)

        for item in s:
            print item, ': ', s[item]
        print '\n'

        if to_file:
            for item in s:
                if type(s[item]) is list:
                    f.write(item + ': ' + ', '.join(s[item]) + '\n')

                else:
                    f.write(item + ': ' + str(s[item]) + '\n')
            f.write('\n')

    if to_file:
        f.close()

def get_fastq_names(path_to_fastq, paired=True):
    # gets fastq filenames in a given directory and runs some simple checks
    # assumes files are compressed with .gz extensions

    R1_files = []
    R2_files = []

    for file in os.listdir(path_to_fastq):

        if file.endswith('.fastq.gz'):
            R = file.split('_')[-2]  # R1 or R2

            if R == 'R1':
                R1_files += [path_to_fastq + file]

            elif R == 'R2' and paired:
                R2_files += [path_to_fastq + file]

            else:
                print 'Unexpected filename structure. Exiting...'
                raise SystemExit

    if len(R1_files) != len(R2_files) and paired:
        print 'Unpaired FASTQ files exist! Check input files.'
        raise SystemExit

    R1_files.sort()
    R2_files.sort()

    return R1_files, R2_files

def generate_samples(R1_files, R2_files):
    # store sample filenames in Sample objects

    samples = []        # list of sample objects

    for i in range(0, len(R1_files) / 2):
        # assign filenames to sample types
        # note: using alphabetization pattern which may not exist in future

        ab_r1 = R1_files[i]
        ab_r2 = R2_files[i]
        panel_r1 = R1_files[i + len(R1_files) / 2]
        panel_r2 = R2_files[i + len(R1_files) / 2]

        sample_num = i + 1
        sample_name = sample_basename + '-' + str(sample_num)

        # set file locations and append to sample object

        panel_r1_temp = temp_dir + panel_r1.split('.fastq.gz')[0].split('/')[-1] + '_temp.fastq'
        panel_r2_temp = temp_dir + panel_r2.split('.fastq.gz')[0].split('/')[-1] + '_temp.fastq'

        ab_r1_temp = temp_dir + ab_r1.split('.fastq.gz')[0].split('/')[-1] + '_temp.fastq'
        ab_r2_temp = temp_dir + ab_r2.split('.fastq.gz')[0].split('/')[-1] + '_temp.fastq'

        panel_barcodes = barcode_dir + sample_name + '_barcodes_panel.json'
        ab_barcodes = barcode_dir + sample_name + '_barcodes_abs.json'

        ab_reads = ab_dir + sample_name + '_ab_reads.tsv'

        samples.append(resources_v2.TapestriSample(sample_num,
                                                   panel_r1,
                                                   panel_r2,
                                                   panel_r1_temp,
                                                   panel_r2_temp,
                                                   ab_r1,
                                                   ab_r2,
                                                   ab_r1_temp,
                                                   ab_r2_temp,
                                                   panel_barcodes,
                                                   ab_barcodes,
                                                   ab_reads))

    return samples

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='''
    dab-seq: single-cell abseq + genotyping
    ben demaree 2019
    
    input requirements:
    -raw fastq files (panel and abs)
    -cell and ab barcode csvs
    -panel bed file
    
    software requirements:
    -gatk 4
    -bowtie2
    -samtools
    -cutadapt
    -bedtools
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('cfg_file', type=str, help='config filename')

    args = parser.parse_args()  # parse arguments

    cfg_f = args.cfg_file  # config filename

    print '''
    
                                                                                                                                              
                                                bbbbbbbb                                                                                      
DDDDDDDDDDDDD                  AAA              b::::::b                                                                                      
D::::::::::::DDD              A:::A             b::::::b                                                                                      
D:::::::::::::::DD           A:::::A            b::::::b                                                                                      
DDD:::::DDDDD:::::D         A:::::::A            b:::::b                                                                                      
  D:::::D    D:::::D       A:::::::::A           b:::::bbbbbbbbb                         ssssssssss       eeeeeeeeeeee       qqqqqqqqq   qqqqq
  D:::::D     D:::::D     A:::::A:::::A          b::::::::::::::bb                     ss::::::::::s    ee::::::::::::ee    q:::::::::qqq::::q
  D:::::D     D:::::D    A:::::A A:::::A         b::::::::::::::::b                  ss:::::::::::::s  e::::::eeeee:::::ee q:::::::::::::::::q
  D:::::D     D:::::D   A:::::A   A:::::A        b:::::bbbbb:::::::b --------------- s::::::ssss:::::se::::::e     e:::::eq::::::qqqqq::::::qq
  D:::::D     D:::::D  A:::::A     A:::::A       b:::::b    b::::::b -:::::::::::::-  s:::::s  ssssss e:::::::eeeee::::::eq:::::q     q:::::q 
  D:::::D     D:::::D A:::::AAAAAAAAA:::::A      b:::::b     b:::::b ---------------    s::::::s      e:::::::::::::::::e q:::::q     q:::::q 
  D:::::D     D:::::DA:::::::::::::::::::::A     b:::::b     b:::::b                       s::::::s   e::::::eeeeeeeeeee  q:::::q     q:::::q 
  D:::::D    D:::::DA:::::AAAAAAAAAAAAA:::::A    b:::::b     b:::::b                 ssssss   s:::::s e:::::::e           q::::::q    q:::::q 
DDD:::::DDDDD:::::DA:::::A             A:::::A   b:::::bbbbbb::::::b                 s:::::ssss::::::se::::::::e          q:::::::qqqqq:::::q 
D:::::::::::::::DDA:::::A               A:::::A  b::::::::::::::::b                  s::::::::::::::s  e::::::::eeeeeeee   q::::::::::::::::q 
D::::::::::::DDD A:::::A                 A:::::A b:::::::::::::::b                    s:::::::::::ss    ee:::::::::::::e    qq::::::::::::::q 
DDDDDDDDDDDDD   AAAAAAA                   AAAAAAAbbbbbbbbbbbbbbbb                      sssssssssss        eeeeeeeeeeeeee      qqqqqqqq::::::q 
                                                                                                                                      q:::::q 
                                                                                                                                      q:::::q 
                                                                                                                                     q:::::::q
                                                                                                                                     q:::::::q
                                                                                                                                     q:::::::q
                                                                                                                                     qqqqqqqqq
                                                                                                                                              
Initializing pipeline...

    '''

    # load config file variables
    # be careful about using exec - it can run bad things
    if not os.path.isfile(cfg_f):
        print 'Config file not found! Please check the file name and path.'
        raise SystemExit

    else:
        with open(cfg_f, 'r') as cfg:
            for line in cfg:
                if line[0] == '[' or line[0] == '#' or line[0] == ' ':
                    continue
                else:
                    var = line.split("#", 1)[0].strip()  # to remove inline comments
                    exec(var)

    # check all files exist
    all_vars = copy.copy(globals())
    input_files = [all_vars[f] for f in all_vars if '_file' in f and f != '__file__']
    input_files.append(bt2_ref + '.1.bt2')
    missing_files = []
    for f in input_files:
        if not os.path.exists(f):
            missing_files.append(f)

    # print missing files, if any, and exit
    if missing_files != []:
        print 'The following input files could not be found:'
        for f in missing_files:
            print f
        print '\nExiting...\n'
        raise SystemExit

    # check that the input fastq directory exists and is not empty
    if not os.path.exists(fastq_dir):
        print 'FASTQ input directory (%s) does not exist! Exiting...\n' % fastq_dir
        raise SystemExit

    elif os.listdir(fastq_dir) == []:
        print 'FASTQ input directory (%s) is empty! Exiting...\n' % fastq_dir
        raise SystemExit

    # create all other directories for this run
    # if it already exists, ignore and continue
    dirs = [all_vars[d] for d in all_vars if '_dir' in d]
    dirs.sort()
    for d in dirs:
        if not os.path.exists(d):
            os.mkdir(d)

    # send slack notification
    start_time = time.time()
    start_time_fmt = str(time.strftime('%m-%d-%Y %H:%M:%S', time.localtime(start_time)))
    slack_message('Pipeline started at %s.' % start_time_fmt)


    print '''
####################################################################################
# Step 1: get input file names and store in TapestriSample objects
####################################################################################
    '''

    # get all fastq filenames
    R1_files, R2_files = get_fastq_names(fastq_dir)

    # store sample info in Sample objects
    samples = generate_samples(R1_files, R2_files)

    # display and write sample summary file
    file_summary(samples, summary, to_file=True)

    print '''
####################################################################################
# Step 2: filter reads for cell barcode and perform error correction
####################################################################################
'''

    # load mission bio barcode csv file
    barcodes = resources_v2.load_barcodes(cell_barcode_csv_file, 1, False)

    # generate hamming dictionary for error correction
    barcodes = resources_v2.generate_hamming_dict(barcodes)

    print '%d unique barcode sequences loaded into dictionary.\n' % len(barcodes)

    # for panel reads, filter reads with valid barcode structure and export to new fastq
    print 'Extracting barcodes from raw fastq files...\n'

    # first, identify reads with valid barcode

    process_barcodes = []

    for sample in samples:
        # panel files
        p = Process(
            target=sample.filter_valid_reads,
            args=(r1_start,
                  barcodes,
                  bar_ind_1,
                  bar_ind_2,
                  'panel'))
        process_barcodes.append(p)
        p.start()

        # ab files
        p = Process(
            target=sample.filter_valid_reads,
            args=(r1_start,
                  barcodes,
                  bar_ind_1,
                  bar_ind_2,
                  'ab'))
        process_barcodes.append(p)
        p.start()

    # wait for processes to finish
    for p in process_barcodes:
        p.join()

    # second, cut adapters from reads and add barcodes to header

    cut_adapters = []

    for sample in samples:
        # panel files
        p = Process(
            target=sample.barcode_reads,
            args=(r1_start,
                  r1_end,
                  r2_end,
                  r1_min_len_panel,
                  r2_min_len_panel,
                  'panel'))
        cut_adapters.append(p)
        p.start()

        # ab files
        p = Process(
            target=sample.barcode_reads,
            args=(r1_start,
                  r1_end,
                  r2_end,
                  r1_min_len_ab,
                  r2_min_len_ab,
                  'ab'))
        cut_adapters.append(p)
        p.start()

    # wait for processes to finish
    for p in cut_adapters:
        p.join()

    print '''
# ###################################################################################
# Step 3: import ab reads, error correct, and count
# ###################################################################################
'''

    # load ab barcode csv file (with descriptions)
    barcodes = resources_v2.load_barcodes(ab_barcode_csv_file, 1, False)
    barcode_descriptions = copy.deepcopy(barcodes)

    # generate hamming dictionary for error correction
    barcodes = resources_v2.generate_hamming_dict(barcodes)

    # process ab reads and look for barcode
    ab_process = []
    for sample in samples:
        # ab files
        p = Process(
            target=sample.process_abs,
            args=(barcodes,
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

    # merge ab reads files into a single file and remove old files
    ab_files = [s.ab_reads for s in samples]

    # concatenate files and remove old ones
    wait([subprocess.Popen('cat %s > %s' % (' '.join(ab_files), ab_reads_merged), shell=True)])
    for f in ab_files:
        try:
            os.remove(f)
        except OSError:
            pass

    # count the ab read umis
    resources_v2.count_umis(ab_reads_merged, umi_counts_merged)

    print '''
###################################################################################
Step 4: count read alignments to inserts
###################################################################################
'''

    # get r1 reads for all panel samples
    panel_r1_files = [s.panel_r1_temp for s in samples]

    # align r1 reads to inserts to obtain read counts across all barcodes
    resources_v2.count_alignments(panel_r1_files, amplicon_file, human_fasta_file, human_genome_index, all_tsv, temp_dir)

    print '''
###################################################################################
Step 5: call valid cells using selected method
###################################################################################
'''

    # call valid cells using cell_caller function
    valid_cells = resources_v2.cell_caller(all_tsv, min_coverage, min_fraction)

    # create SingleCell objects for each valid cell
    cells = [resources_v2.SingleCell(barcode,
                                     by_cell_fastq_dir,
                                     by_cell_bam_dir,
                                     by_cell_gvcf_dir)
             for barcode in valid_cells]

    print '''
###################################################################################
Step 6: write valid cells from panel reads to separate fastq files
###################################################################################
'''

    # check that all barcodes have the same length
    bar_length = list(set([len(k) for k in valid_cells]))
    assert len(bar_length) == 1, 'All barcodes must have the same length!'

    # split cell fastq files for panel reads associated with valid cells
    barcode_names = temp_dir + 'barcodes.txt'
    with open(barcode_names, 'w') as f:
        for b in valid_cells:
            f.write(b + '\n')

    # split files by cell barcode using bbmap demuxbyname.sh
    split_files = []
    for sample in samples:
        demux_cmd = 'demuxbyname.sh prefixmode=f length=%d in1=%s in2=%s out=%s names=%s' % (
            bar_length[0],
            sample.panel_r1_temp,
            sample.panel_r2_temp,
            by_cell_fastq_dir + '%.fastq',
            barcode_names
        )
        p = subprocess.Popen(demux_cmd, shell=True)
        split_files.append(p)

    # wait for processes to finish
    wait(split_files)

    print '''
####################################################################################
# Step 7: align, convert, sort, index panel reads
####################################################################################
'''

    # preprocess cells in chunks (hardware-limited)
    # size of chunks for preprocessing batching (based on hardware limitations)
    samples_per_chunk = 300
    cells_processed = 0

    # split cell list into chunks
    sample_chunks = [cells[i:i + samples_per_chunk] for i in xrange(0, len(cells), samples_per_chunk)]

    for chunk in sample_chunks:
        # preprocess cells for variant calling
        preprocess_cells = [resources_v2.SingleCell.align_sample(cell, bt2_ref) for cell in chunk]
        # wait for all processes to finish before continuing
        wait(preprocess_cells)

        # index bam files
        index_bam = [resources_v2.SingleCell.index_bam(cell) for cell in chunk]
        # wait for all processes to finish before continuing
        wait(index_bam)

        cells_processed += samples_per_chunk

        print '\n%d of %d cells pre-processed.' % (cells_processed, len(cells))

    print '''
####################################################################################
# Step 8: perform variant calling for all cells
####################################################################################
'''

    # call variants on cells in chunks (hardware-limited)
    # size of chunks for variant call batching (based on hardware limitations)
    samples_per_chunk = 50
    cells_processed = 0

    # split cell list into chunks
    sample_chunks = [cells[i:i + samples_per_chunk] for i in xrange(0, len(cells), samples_per_chunk)]

    for chunk in sample_chunks:
        # perform variant calling on cells
        variant_calling = [resources_v2.SingleCell.call_variants(cell,
                                                                 human_fasta_file,
                                                                 interval_file,
                                                                 dbsnp_file)
                           for cell in chunk]

        # wait for all processes to finish before continuing
        wait(variant_calling)

        cells_processed += samples_per_chunk

        print '\n%d of %d cells genotyped.\n' % (cells_processed, len(cells))

    print '''
####################################################################################
# Step 9: combine gvcf files
####################################################################################
'''

    # combine gvcfs in chunks (hardware-limited)
    # size of chunks for gvcf merger batching (based on hardware limitations)
    samples_per_chunk = 200
    cells_processed = 0

    # split cell list into chunks
    sample_chunks = [cells[i:i + samples_per_chunk] for i in xrange(0, len(cells), samples_per_chunk)]
    chunk_number = 0

    # first stage merger
    for chunk in sample_chunks:
        # combine gvcfs in batches
        combine_gvcfs = [resources_v2.SingleCell.combine_gvcfs(chunk,
                                                               chunk_number,
                                                               human_fasta_file,
                                                               interval_file,
                                                               dbsnp_file,
                                                               merged_gvcf_dir,
                                                               genotyping_dir)]

        # wait for all processes to finish before continuing
        wait(combine_gvcfs)

        chunk_number += 1
        cells_processed += samples_per_chunk

        print '\n%d of %d gvcfs combined (first stage).\n' % (cells_processed, len(cells))

    # second stage merger
    # find gvcf files from first merger stage

    gvcfs_to_merge = [merged_gvcf_dir + f for f in os.listdir(merged_gvcf_dir) if f.split('.')[-1] == 'vcf']

    # perform final merger

    print 'Performing final merger on %d GVCF files of %d cells each.' % (len(gvcfs_to_merge), samples_per_chunk)

    combine_gvcfs = [resources_v2.SingleCell.combine_gvcfs(gvcfs_to_merge,
                                                           merged_gvcf,
                                                           human_fasta_file,
                                                           interval_file,
                                                           dbsnp_file,
                                                           genotyping_dir,
                                                           genotyping_dir,
                                                           multi_sample=True)]

    # wait for all processes to finish before continuing
    wait(combine_gvcfs)

    print '''
####################################################################################
# Step 10: genotype merged gvcfs
####################################################################################
'''
    genotype_gvcfs = [resources_v2.SingleCell.genotype_gvcfs(human_fasta_file,
                                                             dbsnp_file,
                                                             merged_gvcf,
                                                             geno_gvcf)]

    # wait for all processes to finish before continuing
    wait(genotype_gvcfs)

    print 'Pipeline complete!'

    # send slack notification
    elapsed_time = time.time() - start_time
    elapsed_time_fmt = str(time.strftime('%Hh %Mm %Ss', time.gmtime(elapsed_time)))
    slack_message('Pipeline complete! Total elapsed time is %s.' % elapsed_time_fmt)