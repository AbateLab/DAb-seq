"""

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 7.9.2019

functions required for the processing pipeline

two classes are defined:
* TapestriTube: variables and methods that operate on the level of a single Tapestri tube library
* SingleCell: variables and methods that operate on the level of one cell barcode

"""

import os
import os.path
import shutil
import csv
from itertools import product, combinations, izip
from collections import Counter
import numpy as np
import subprocess
import sys
import copy
import allel
import pandas as pd
import h5py

# add the modified umi_tools directory to python path
sys.path.append(os.path.join(sys.path[0], 'umi_tools'))
import Utilities as U
import network as network

# define classes

class TapestriTube(object):
    # class for storing metadata for each tapestri tube

    def __init__(self,
                 tube_num,
                 panel_r1,
                 panel_r2,
                 panel_r1_temp,
                 panel_r2_temp,
                 ab_r1,
                 ab_r2,
                 ab_r1_temp,
                 ab_r2_temp,
                 ab_reads,
                 umi_counts):

        self.tube_num = tube_num            # number identifying tube

        self.panel_r1 = panel_r1                # panel R1 fastq path
        self.panel_r2 = panel_r2                # panel R2 fastq path
        self.panel_r1_temp = panel_r1_temp      # temp panel R1 fastq path
        self.panel_r2_temp = panel_r2_temp      # temp panel R2 fastq path

        self.ab_r1 = ab_r1                      # ab R1 fastq path
        self.ab_r2 = ab_r2                      # ab R2 fastq path
        self.ab_r1_temp = ab_r1_temp            # temp ab R1 fastq path
        self.ab_r2_temp = ab_r2_temp            # temp ab R2 fastq path

        self.ab_reads = ab_reads                # file containing filtered ab reads
        self.umi_counts = umi_counts            # file containing umi counts by barcode

    def barcode_reads(self,
                      r1_start,
                      r1_end,
                      r2_end,
                      mb_barcodes,
                      bar_ind_1,
                      bar_ind_2,
                      r1_min_len,
                      r2_min_len,
                      library_type):
        # for valid reads, add barcode header to fastq file and trim

        assert library_type == 'ab' or library_type == 'panel', 'Library type must be panel or ab!'

        # set filenames according to library type (panel or ab)
        if library_type == 'panel':

            r1_in = self.panel_r1
            r2_in = self.panel_r2

            r1_out = open(self.panel_r1_temp, 'w')
            r2_out = open(self.panel_r2_temp, 'w')

        elif library_type == 'ab':

            r1_in = self.ab_r1
            r2_in = self.ab_r2

            r1_out = open(self.ab_r1_temp, 'w')
            r2_out = open(self.ab_r2_temp, 'w')

        # trim 5' end of read and check barcode
        bar_cmd = 'cutadapt -a %s -O 8 -e 0.2 %s -j 8 --quiet' % (r1_start,
                                                                  r1_in)
        bar_file = subprocess.Popen(bar_cmd, stdout=subprocess.PIPE, shell=True)

        # hard trimming (55 bp off R1, 5 bp off R2) is used to ensure entire cell barcode region is removed
        # barcode bases for V1 chemistry: 47 bp
        # max barcode bases for V2 chemistry: 50 bp

        # trimmed reads output
        trim_cmd = 'cutadapt -a %s -A %s --interleaved -j 8 -u 55 -U 5 -n 2 %s %s --quiet' % (r1_end,
                                                                                              r2_end,
                                                                                              r1_in,
                                                                                              r2_in)
        trim_file = subprocess.Popen(trim_cmd, stdout=subprocess.PIPE, shell=True)

        total_reads = 0     # total count of all reads
        invalid_reads = 0   # no valid barcode found
        valid_reads = 0     # valid, barcoded read count
        too_short = 0

        # iterate through info file (barcodes) and trim file (reads)
        for bar_line, trim_line in izip(bar_file.stdout, trim_file.stdout):

            assert bar_line.strip() == trim_line.strip(), 'Cluster IDs do not match!'

            total_reads += 1

            # extract trimmed barcode region
            bar_seq = bar_file.stdout.next().strip()

            # if trimmed adapter is too long - no barcode present
            if len(bar_seq) > 52:

                invalid_reads += 1

                # advance through files
                for i in range(7):
                    trim_file.stdout.next()
                for i in range(2):
                    bar_file.stdout.next()

                continue

            # contains valid adapter
            else:
                # find barcodes and check that they are a valid MB barcode
                check = check_seq(bar_seq, bar_ind_1, bar_ind_2, mb_barcodes)

                # not a valid barcode
                if check == 'fail':

                    invalid_reads += 1

                    # advance through files
                    for i in range(7):
                        trim_file.stdout.next()
                    for i in range(2):
                        bar_file.stdout.next()

                    continue

                # valid barcode
                else:

                    # adavance through barcode file
                    for i in range(2):
                        bar_file.stdout.next()

                    barcode = check[0] + check[1] + '-' + str(self.tube_num)

                    # R1 from trimmed file
                    header_1 = trim_line.strip()
                    id_1 = header_1.split(' ')[0][1:]
                    seq_1 = trim_file.stdout.next().strip()
                    trim_file.stdout.next()
                    qual_1 = trim_file.stdout.next().strip()

                    # R2 from trimmed file
                    header_2 = trim_file.stdout.next().strip()
                    id_2 = header_2.split(' ')[0][1:]

                    assert id_1 == id_2, 'Cluster IDs in interleaved input do not match!'

                    seq_2 = trim_file.stdout.next().strip()
                    trim_file.stdout.next()
                    qual_2 = trim_file.stdout.next().strip()

                    # check reads for length
                    if len(seq_1) < r1_min_len or len(seq_2) < r2_min_len:
                        too_short += 1
                        continue

                    # add barcoded headers and reads to file
                    id = '@' + id_1 + '_' + barcode

                    # write to output fastq files
                    r1_out.write('%s\n%s\n+\n%s\n' % (id, seq_1, qual_1))
                    r2_out.write('%s\n%s\n+\n%s\n' % (id, seq_2, qual_2))

                    valid_reads += 1

                    # print counter
                    if valid_reads % 1e6 == 0:
                        print('Tube %d-%s: %d valid trimmed pairs saved to file.' % (self.tube_num,
                                                                                       library_type,
                                                                                       valid_reads))

    def process_abs(self,
                    ab_barcodes,
                    barcode_descriptions,
                    ab_handles,
                    ab_bar_coord,
                    ab_umi_coord,
                    min_umi_qual):
        # extract ab barcodes and umis from raw ab reads

        # write passed ab reads to file
        ab_reads_file = open(self.ab_reads, 'w')

        # use cutadapt to select reads with correct structure
        ab_cmd = 'cutadapt -j 8 %s -O 12 -e 0.2 -n 2 %s --quiet --discard-untrimmed' % (ab_handles,
                                                                                         self.ab_r2_temp)

        ab_process = subprocess.Popen(ab_cmd, stdout=subprocess.PIPE, shell=True)

        valid_ab_reads = 0

        # iterate through ab reads with correct adapters
        for line in ab_process.stdout:

            cell_barcode = line.strip().split('_')[1]  # extract corrected barcode from header

            # extract sequences
            seq = ab_process.stdout.next().strip()
            ab_process.stdout.next()
            qual = ab_process.stdout.next().strip()

            # check trimmed read length
            if len(seq) != len(ab_bar_coord + ab_umi_coord):
                continue

            # check ab barcode is valid
            bar = ''.join([seq[i] for i in ab_bar_coord])
            bar = correct_barcode(ab_barcodes, bar)
            if bar == 'invalid':
                continue

            # check umi quality
            umi = ''.join([seq[i] for i in ab_umi_coord])
            umi_qual = [ord(qual[i]) - 33 for i in ab_umi_coord]
            if not all(q >= min_umi_qual for q in umi_qual):
                continue

            # if a read passes all filters, write it to file
            valid_ab_reads += 1
            ab_reads_file.write(cell_barcode + '\t')
            ab_reads_file.write(barcode_descriptions[bar] + '\t')
            ab_reads_file.write(umi + '\n')

        ab_reads_file.close()

    def count_umis(self, clustering_method):
        # count umis using selected clustering methods
        # assumes ab reads file is sorted

        # extract ab umis from ab reads file
        ab_reads = open(self.ab_reads, 'r')

        # write to umi counts file
        umi_counts = open(self.umi_counts, 'w')
        umi_counts.write('cell_barcode\tab_description\traw\tunique')
        if clustering_method == 'all':
            umi_counts.write('\tadjacency\tdirectional\tpercentile\tcluster')

        # initialize with group id (cell barcode + ab) and group umis

        group_id_prev = ''
        group_umis = []
        first_line = True

        for line in ab_reads:

            group_id_curr = line.split('\t')[:2]
            umi = line.strip().split('\t')[2]

            # still in same umi group
            if group_id_curr == group_id_prev or first_line:
                group_umis.append(umi)
                first_line = False

            # new umi group
            else:
                # cluster umis from previous group
                counts = self.umi_tools_cluster(group_umis, clustering_method)

                # write umi count to file
                self.write_umis(counts, group_id_prev, umi_counts, clustering_method)

                # reset group umis
                group_umis = [umi]

            group_id_prev = group_id_curr

        # process final group
        counts = self.umi_tools_cluster(group_umis, clustering_method)
        self.write_umis(counts, group_id_prev, umi_counts, clustering_method)
        umi_counts.write('\n')
        umi_counts.close()

    @staticmethod
    def umi_tools_cluster(group_umis, clustering_method):
        # cluster a umi group with umi-tools

        umi_counter = Counter(group_umis)

        # set up UMIClusterer functor with parameters specific to specified method
        # choose method = 'all' for all available methods

        processor = network.UMIClusterer()  # initialize UMIclusterer
        clusters = processor(umi_counter.keys(),
                             umi_counter,
                             threshold=1,
                             cluster_method=clustering_method)
        counts = {k: str(len(v)) for k, v in clusters.items()}

        counts['raw'] = str(len(group_umis))

        return counts

    @staticmethod
    def write_umis(counts, group_id, umi_counts_handle, clustering_method):
        # write a umi group to file

        umi_counts_handle.write('\n' + '\t'.join(group_id) + '\t' + counts['raw'] + '\t' + counts['unique'])

        if clustering_method == 'all':
            umi_counts_handle.write('\t' + counts['adjacency']
                                    + '\t' + counts['directional']
                                    + '\t' + counts['percentile']
                                    + '\t' + counts['cluster'])

class SingleCell(object):
    # class for storing metadata for each single cell file

    def __init__(self, cell_barcode, fastq_dir, bam_dir, vcf_dir, flt3_dir):
        # initialize object by generating filenames

        self.cell_barcode = cell_barcode                    # cell barcode
        self.fastq = fastq_dir + cell_barcode + '.fastq'    # fastq file

        assert self.fastq, 'fastq file does not exist'

        self.bam = bam_dir + cell_barcode + '.bam'          # bam file
        self.bai = bam_dir + cell_barcode + '.bai'          # bam file index
        self.vcf = vcf_dir + cell_barcode + '.g.vcf'        # gvcf file
        self.flt3_vcf = flt3_dir + cell_barcode + '.flt3itd.vcf'     # flt3-itd vcf

        self.valid = False      # marker for valid cells
        self.alignments = {}    # alignment counts for each interval

    def align_and_index(self, bt2_ref):

        # align the panel to the bowtie2 index and generate sorted bam file
        # read filters: read mapped, mapq >= 3, primary alignment
        align_cmd = 'bowtie2 -x %s --mm --interleaved %s' \
                    ' --rg-id %s --rg SM:%s --rg PL:ILLUMINA --rg CN:UCSF --quiet' \
                    ' | samtools view -b -q 3 -F 4 -F 0X0100' \
                    ' | samtools sort -o %s' \
                    % (bt2_ref,
                       self.fastq,
                       self.cell_barcode,
                       self.cell_barcode,
                       self.bam)
        subprocess.call(align_cmd, shell=True)

        # index all bam files using samtools
        index_cmd = 'samtools index %s %s' \
                    % (self.bam,
                       self.bai)
        subprocess.call(index_cmd, shell=True)

    def call_flt3(self, fasta):
        # call flt3-itds using ITDseek (https://github.com/tommyau/itdseek)

        flt3_cmd = 'itdseek.sh %s %s samtools > %s' \
                       % (self.bam,
                          fasta,
                          self.flt3_vcf)

        subprocess.call(flt3_cmd, shell=True)

    def call_variants(self, fasta, interval_file, ploidy):
        # call variants using gatk
        # setting --max-reads-per-alignment-start to use downsampling

        variants_cmd = 'gatk HaplotypeCaller -R %s -I %s -O %s -L %s ' \
                       '--emit-ref-confidence GVCF ' \
                       '--verbosity ERROR ' \
                       '--native-pair-hmm-threads 1 ' \
                       '--max-alternate-alleles 2 ' \
                       '--standard-min-confidence-threshold-for-calling 0 ' \
                       '--max-reads-per-alignment-start 200 ' \
                       '--minimum-mapping-quality 3 ' \
                       '--sample-ploidy %s' \
                       % (fasta,
                          self.bam,
                          self.vcf,
                          interval_file,
                          ploidy)

        subprocess.call(variants_cmd, shell=True)

def count_alignments(r1_files, amplicon_file, fasta_file, tsv, aln_stats_file, temp_dir):
    # align and count r1 reads for all sample barcodes, and save to tsv file

    # get fasta file from genome for this interval
    insert_fasta = temp_dir + 'amplicons.fasta'
    subprocess.call('bedtools getfasta -fi %s -bed %s -fo %s -name' % (fasta_file, amplicon_file, insert_fasta), shell=True)

    # build bt2 index for this fasta
    insert_bt2 = temp_dir + 'inserts'
    subprocess.call('bowtie2-build %s %s' % (insert_fasta, insert_bt2), shell=True)

    # extract names of reference sequences
    refs = []
    get_refs = subprocess.Popen('bowtie2-inspect -n %s' % insert_bt2, stdout=subprocess.PIPE, shell=True)
    for line in get_refs.stdout:
        refs.append(line.strip())
    refs.sort()
    refs_dict = dict(zip(refs, [0] * len(refs)))

    # amplicon dict (key: cell barcode; value: list of amplicon counts)
    amplicons = {}

    # align reads with bowtie2
    # read filters: read mapped, mapq >= 3, primary alignment
    # this filter is the same as used for single-cell mapping
    bt2_cmd = 'bowtie2 -p 24 -x %s -U %s 2> %s | samtools view -q 3 -F 4 -F 0X0100' % (insert_bt2,
                                                                                       ' -U '.join(r1_files),
                                                                                       aln_stats_file)
    bt2_align = subprocess.Popen(bt2_cmd, stdout=subprocess.PIPE, shell=True)

    # iterate through all reads
    for line in bt2_align.stdout:

        # if read passes all filters, extract barcode
        record = line.split('\t')
        query_name = record[0]
        cell_barcode = query_name.split('_')[1]
        reference_name = record[2]

        # save alignment to dict
        try:
            amplicons[cell_barcode][reference_name] += 1

        except KeyError:
            amplicons[cell_barcode] = copy.deepcopy(refs_dict)
            amplicons[cell_barcode][reference_name] += 1

    # write amplicon tsv to file for this sample
    with open(tsv, 'w') as f:
        f.write('cell_barcode\t' + '\t'.join(refs) + '\n')
        for c in amplicons:
            f.write(c + '\t' + '\t'.join([str(amplicons[c][r]) for r in refs]) + '\n')

def umi_counts_by_cell(umi_counts_merged, ab_barcode_csv_file, ab_dir, cells=None):
    # converts umi group counts to counts by cell

    # create by_method directory if it doesn't exist
    if cells is None:
        ab_dir_by_method = ab_dir + 'by_method.ALL/'
    else:
        ab_dir_by_method = ab_dir + 'by_method.CELLS/'

    if not os.path.exists(ab_dir_by_method):
        os.mkdir(ab_dir_by_method)
    else:
        shutil.rmtree(ab_dir_by_method)
        os.mkdir(ab_dir_by_method)

    # load umi counts
    umi_counts = pd.read_csv(umi_counts_merged, sep='\t', header=0, index_col=0)

    # filter on valid cells
    if cells is not None:
        valid_barcodes = sorted([c.cell_barcode for c in cells])
        umi_counts = umi_counts[umi_counts.index.isin(valid_barcodes)]

    # get all ab barcode descriptions in experiment
    ab_barcodes = sorted(load_barcodes(ab_barcode_csv_file).values())

    # for each method, save in a dataframe
    count_methods = [m for m in umi_counts.columns if 'ab_desc' not in m]

    for m in count_methods:
        ab_counts_by_method = pd.pivot_table(umi_counts, values=m, index=['cell_barcode'], columns=['ab_description'],
                                             aggfunc=np.sum).fillna(0)

        # reorder the cell barcodes
        # if a cell is missing (unlikely but possible), fill with 0
        if cells is not None:
            ab_counts_by_method = ab_counts_by_method.reindex(valid_barcodes, fill_value=0)

        # if an ab is missing (also unlikely), fill with 0
        for ab in ab_barcodes:
            if ab not in ab_counts_by_method:
                ab_counts_by_method[ab] = 0

        # sort by column name and save umi counts file
        ab_counts_by_method.sort_index(axis=1, inplace=True)
        ab_counts_by_method = ab_counts_by_method.astype(int)
        if cells is not None:
            count_method_file = os.path.basename(umi_counts_merged)[:-3] + m + '.cells.tsv'
        else:
            count_method_file = os.path.basename(umi_counts_merged)[:-3] + m + '.all.tsv'
        ab_counts_by_method.to_csv(path_or_buf=ab_dir_by_method+count_method_file, sep='\t')

def load_barcodes(barcode_file, max_dist=1, check_barcodes=True):
    # loads barcodes from csv and checks all pairwise distances to ensure error correction will work
    # returns a dictionary of barcodes with their descriptions

    # load barcodes from csv
    reader = csv.reader(open(barcode_file, 'r'))
    barcodes = {}
    for barcode, desc in reader:
        barcodes[barcode] = desc

    # optional distance check for barcodes (can turn off once validated)
    if check_barcodes:

        # check all barcodes have same length
        lengths = map(len, barcodes.keys())
        if len(set(lengths)) != 1:
            print('Barcodes must all be same length! Exiting...')
            raise SystemExit

        # check pairwise hamming distances
        dist_req = 2 * max_dist + 1         # for this max_dist, need this distance between all barcodes
        pairs = list(combinations(barcodes.keys(), 2))
        for pair in pairs:
            if hd(pair[0], pair[1]) < dist_req:
                print('Error: The edit distance between barcodes %s and %s is less than %d.\n'
                      'An error correction of %d bases will not work.' % (pair[0], pair[1], dist_req, max_dist))

    return barcodes

def hd(s1, s2):
    # calculate the hamming distance between two strings of equal length

    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def generate_hamming_dict(barcodes):
    # return dictionary of strings within 1 hamming distance of barcodes

    for barcode in barcodes:
        hd_1 = sorted(hamming_circle(barcode, 1))
        barcodes[barcode] = dict(zip(hd_1, [1] * len(hd_1)))

    return barcodes

def hamming_circle(s, n):
    # generate strings over alphabet whose hamming distance from s is exactly n
    # from https://codereview.stackexchange.com/a/88919

    alphabet = 'ATCG'
    s = s.upper()

    for positions in combinations(range(len(s)), n):

        for replacements in product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)

            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]

                else:
                    cousin[p] = alphabet[r]

            yield ''.join(cousin)

def check_seq(seq, bar_ind_1, bar_ind_2, barcodes):
    # checks a sequence for valid barcodes
    # outputs a list of raw and corrected information if valid, 'fail' otherwise

    try:
        bar_1 = ''.join([seq[i] for i in bar_ind_1])
        bar_2 = ''.join([seq[i] for i in bar_ind_2])

    except IndexError:
        return 'fail'

    # try and error correct the ab barcode
    corr_barcode_1 = correct_barcode(barcodes, bar_1)
    corr_barcode_2 = correct_barcode(barcodes, bar_2)

    if (corr_barcode_1 == 'invalid') or (corr_barcode_2 == 'invalid'):
        return 'fail'

    return [corr_barcode_1, corr_barcode_2]

def correct_barcode(barcodes, raw_barcode):
    # attempts to correct a raw barcode using list of valid barcodes
    # return: raw barcode, corrected barcode

    # check if barcode is an exact match
    if raw_barcode in barcodes:
        return raw_barcode

    # attempt to error correct barcode
    else:
        for valid_barcode in barcodes:
            if raw_barcode in barcodes[valid_barcode]:
                return valid_barcode

    # if correction fails
    return 'invalid'

def db_import(db_path, sample_map_path, interval):
    # import single-cell gvcfs for one interval

    # make call to gatk for genomics db import
    # batch size is number of cells
    db_import_cmd = 'gatk --java-options "-Xmx4g" GenomicsDBImport ' \
                    '--genomicsdb-workspace-path %s ' \
                    '--batch-size 50 ' \
                    '--reader-threads 2 ' \
                    '--validate-sample-name-map true ' \
                    '-L %s ' \
                    '--sample-name-map %s' % (db_path,
                                              interval,
                                              sample_map_path)

    process = subprocess.call(db_import_cmd, shell=True)

def joint_genotype(db_path, fasta, interval, output_vcf, ploidy):
    # perform joint genotyping across a cohort using data from a genomicsdb store

    # make call to gatk for genotyping
    genotype_cmd = 'gatk --java-options "-Xmx4g" GenotypeGVCFs ' \
                   '-V %s ' \
                   '-R %s ' \
                   '-L %s ' \
                   '-O %s ' \
                   '--include-non-variant-sites ' \
                   '--sample-ploidy %s' \
                   % (db_path,
                      fasta,
                      interval,
                      output_vcf,
                      ploidy)

    process = subprocess.call(genotype_cmd, shell=True)

def left_align_trim(ref_fasta, geno_vcf, split_vcf):
    # uses bcftools to split multiallelics, left-align, and trim

    # split and left-align variants
    split_cmd = 'bcftools norm --threads 16 -f %s --check-ref w -m - %s > %s' %\
                (ref_fasta, geno_vcf, split_vcf)

    subprocess.call(split_cmd, shell=True)

def snpeff_annotate(snpeff_summary, snpeff_config, split_vcf, snpeff_annot_vcf):
    # annotate a vcf file with snpeff functional predictions

    annotate_cmd = 'snpEff ann -v -stats %s -c %s hg19 %s > %s' % (snpeff_summary,
                                                                   snpeff_config,
                                                                   split_vcf,
                                                                   snpeff_annot_vcf)
    subprocess.call(annotate_cmd, shell=True)

def bcftools_annotate(annotations_vcf, input_vcf, column_info, output_vcf):
    # uses bcftools to annotate a vcf file with annotation information from another
    # input and output vcf should both be uncompressed vcf

    # bgzip and index the input vcf
    subprocess.call('bgzip -f -@ 16 %s' % input_vcf, shell=True)
    subprocess.call('tabix -f %s' % input_vcf + '.gz', shell=True)

    # use bcftools to annotate the input
    bcf_cmd = 'bcftools annotate -a %s %s %s > %s' % (annotations_vcf,
                                                      column_info,
                                                      input_vcf + '.gz',
                                                      output_vcf)
    subprocess.call(bcf_cmd, shell=True)

def vcf_to_tables(vcf_file, genotype_file, variants_tsv, ploidy, itd_vcf_file=False, non_human=False):
    # parses a vcf file into a series of tables
    # if itd_files is given, adds flt3 itd variants to table

    # load vcf file into numpy array
    if non_human:
        # if non-human, exclude annotation info from snpeff
        vcf = allel.read_vcf(vcf_file,
                             fields=['variants/*',
                                     'calldata/GT',
                                     'calldata/AD',
                                     'calldata/GQ',
                                     'calldata/DP',
                                     'samples'])
    else:
        # if human, include annotation info from snpeff
        vcf = allel.read_vcf(vcf_file,
                             transformers=allel.ANNTransformer(),
                             fields=['variants/*',
                                     'calldata/GT',
                                     'calldata/AD',
                                     'calldata/GQ',
                                     'calldata/DP',
                                     'samples',
                                     'ANN'])

    # layers to extract:
    # GT: genotype (0: WT, 1: HET, 2: HOM, 3: no call) [diploid]
    # DP: total read depth
    # GQ: genotype quality
    # AD: alt allele depth
    # RD: ref allele depth

    GT = np.sum(vcf['calldata/GT'], axis=2)
    GT[GT == (-1 * ploidy)] = 3
    DP = np.stack(vcf['calldata/DP'], axis=0)
    GQ = np.stack(vcf['calldata/GQ'], axis=0)
    AD = np.stack(vcf['calldata/AD'][:, :, 1], axis=0)
    RD = np.stack(vcf['calldata/AD'][:, :, 0], axis=0)

    # non-human reference
    if non_human:
        # create variant names
        names = [vcf['variants/CHROM'][i] +
                 ':' + str(vcf['variants/POS'][i]) +
                 ':' + vcf['variants/REF'][i] +
                 '/' + vcf['variants/ALT'][:, 0][i]
                 for i in range((vcf['variants/REF'].shape[0]))]

    # human reference
    else:
        # create variant names
        names = [vcf['variants/ANN_Gene_Name'][i] +
                 ':' + vcf['variants/CHROM'][i] +
                 ':' + str(vcf['variants/POS'][i]) +
                 ':' + vcf['variants/REF'][i] +
                 '/' + vcf['variants/ALT'][:, 0][i]
                 for i in range((vcf['variants/REF'].shape[0]))]

        # assemble and save variant annotations to file
        variants_table = pd.DataFrame(data=names, columns=['Name'])

        # clinvar variation id
        variants_table['ClinVar_Variation_ID'] = vcf['variants/ID']

        # snpeff columns
        ANN_columns = [c for c in list(vcf) if '/ANN' in c]
        for ann in ANN_columns:
            variants_table['SnpEff_' + ann.split('/ANN_')[1]] = vcf[ann]

        # clinvar columns
        CLN_columns = [c for c in list(vcf) if '/CLN' in c]
        for cln in CLN_columns:
            variants_table['ClinVar_' + cln.split('/')[1]] = vcf[cln]

        # optional: add flt3-itd variants to table
        if itd_vcf_file:

            # make sure flt3 vcf is not empty
            empty = True
            with open(itd_vcf_file, 'r') as f:
                for line in f:
                    if line[0] != '#':
                        empty = False
                        break

            if not empty:

                itd_vcf = allel.read_vcf(itd_vcf_file, fields=['*'])

                # create itd variant names
                itd_names = ['FLT3-ITD' +
                             ':' + itd_vcf['variants/CHROM'][i] +
                             ':' + str(itd_vcf['variants/POS'][i]) +
                             ':' + itd_vcf['variants/REF'][i] +
                             '/' + itd_vcf['variants/ALT'][:, 0][i]
                             for i in range((itd_vcf['variants/REF'].shape[0]))]

                # add itd variant rows to variants table
                itd_table = pd.DataFrame(data=list(set(itd_names)), columns=['Name'])
                names += list(set(itd_names))
                variants_table = pd.concat([variants_table, itd_table], sort=True)

                # add itd variants to other layers
                # set RD = AD and GQ = 100 when itd is present
                # default for GT is 'WT' (0) - no ITD present

                # create additional array entries
                GT = np.concatenate((GT, np.zeros((itd_table.shape[0], GT.shape[1]))), axis=0)
                GQ = np.concatenate((GQ, np.zeros((itd_table.shape[0], GQ.shape[1]))), axis=0)
                DP = np.concatenate((DP, np.zeros((itd_table.shape[0], DP.shape[1]))), axis=0)
                AD = np.concatenate((AD, np.zeros((itd_table.shape[0], AD.shape[1]))), axis=0)
                RD = np.concatenate((RD, np.zeros((itd_table.shape[0], RD.shape[1]))), axis=0)

                # indices for adding entries to arrays
                var_ind = dict(zip(names, range(len(names))))
                bar_ind = dict(zip(vcf['samples'], range(len(vcf['samples']))))

                # for each cell barcode, add entry to genotyping array
                for i in range(len(itd_vcf['variants/ID'])):
                    cell_barcode = itd_vcf['variants/ID'][i]
                    alt_depth = itd_vcf['variants/QUAL'][i]
                    vaf = itd_vcf['variants/VAF'][i]
                    total_depth = int(round(np.true_divide(alt_depth, vaf)))

                    # set GT according to vaf
                    # het mut
                    if vaf < 0.9:
                        geno = 1

                    # hom mut
                    else:
                        geno = 2

                    # store entries in genotyping array
                    GT[var_ind[itd_names[i]], bar_ind[cell_barcode]] = geno
                    GQ[var_ind[itd_names[i]], bar_ind[cell_barcode]] = 100
                    DP[var_ind[itd_names[i]], bar_ind[cell_barcode]] = total_depth
                    RD[var_ind[itd_names[i]], bar_ind[cell_barcode]] = total_depth
                    AD[var_ind[itd_names[i]], bar_ind[cell_barcode]] = alt_depth

        # save variants to file
        variants_table.to_csv(path_or_buf=variants_tsv, sep='\t', index=False, encoding='utf-8')

    # encode variant names and cell barcodes
    names = [n.encode('utf8') for n in names]
    barcodes = [b.encode('utf8') for b in vcf['samples']]

    # save genotyping information to compressed hdf5 file
    # transpose to yield cells as rows, variants as columns
    with h5py.File(genotype_file, 'w') as f:
        f.create_dataset('GT', data=np.transpose(GT), dtype='i1', compression='gzip')
        f.create_dataset('GQ', data=np.transpose(GQ), dtype='i1', compression='gzip')
        f.create_dataset('DP', data=np.transpose(DP), dtype='i2', compression='gzip')
        f.create_dataset('AD', data=np.transpose(AD), dtype='i2', compression='gzip')
        f.create_dataset('RD', data=np.transpose(RD), dtype='i2', compression='gzip')
        f.create_dataset('VARIANTS', data=names, compression='gzip')
        f.create_dataset('CELL_BARCODES', data=barcodes, compression='gzip')

def add_hdf5_ab_counts(geno_hdf5, sample_names, ab_count_dirs):
    # add ab counts to compressed hdf5 file

    # get list of valid cell barcodes
    with h5py.File(geno_hdf5, 'r') as f:
        cell_barcodes = copy.deepcopy([c.decode('utf8') for c in f['CELL_BARCODES']])

    # combine abs from multiple samples into same dataframe, appending sample name to barcode
    count_methods_by_sample = {}    # key: method, value: list of dataframes

    for i in range(len(sample_names)):

        ab_filenames = [ab_count_dirs[i] + f for f in os.listdir(ab_count_dirs[i])]
        count_methods = [f.split('.')[-3] for f in ab_filenames]

        for j in range(len(count_methods)):

            sample_count_df = pd.read_csv(ab_filenames[j], sep='\t', header=0, index_col=0)
            sample_count_df.index = sample_count_df.index + '-' + sample_names[i]

            try:
                count_methods_by_sample[count_methods[j]] = \
                    count_methods_by_sample[count_methods[j]].append(sample_count_df)

            except KeyError:
                count_methods_by_sample[count_methods[j]] = sample_count_df

    # check that all methods have all cell barcodes present
    for m in count_methods_by_sample:
        assert set(list(count_methods_by_sample[m].index)) == set(cell_barcodes), 'Sample cell barcodes do not match!'

    # check that all methods have identical abs
    ab_columns = [list(count_methods_by_sample[m].columns) for m in count_methods_by_sample]
    assert all(x == ab_columns[0] for x in ab_columns), 'Ab column names do not match!'

    # save list of ab descriptions to the hdf5 file
    with h5py.File(geno_hdf5, 'a') as f:
        ab_names = [a.encode('utf8') for a in list(count_methods_by_sample['raw'].columns)]
        f.create_dataset('AB_DESCRIPTIONS', data=ab_names, compression='gzip')

    # save the ab count tables to the hdf5 file
    # reorder the ab count tables to the same as the cell barcodes list
    for m in count_methods_by_sample:
        with h5py.File(geno_hdf5, 'a') as f:
            dataset = np.asarray(count_methods_by_sample[m].reindex(cell_barcodes))
            f.create_dataset('ABS/'+m, data=dataset, dtype='i4', compression='gzip')

def file_summary(tubes, cfg_file, summary_file, header_file=False):
    # displays sample files and class variables

    input_msg = '''
####################################################################################
# INPUT FILE SUMMARY
####################################################################################

'''
    print(input_msg)

    f = open(summary_file, 'w')

    if header_file:
        with open(header_file, 'r') as h:
            for line in h:
                f.write(line)

    f.write(input_msg)

    # write file list to run summary

    for tube in tubes:

        s = vars(tube)

        for item in s:
            if type(s[item]) is list:
                to_write = item + ': ' + ', '.join(s[item]) + '\n'
                print(to_write.strip())
                f.write(to_write)

            else:
                to_write = item + ': ' + str(s[item]) + '\n'
                print(to_write.strip())
                f.write(to_write)

        print('')
        f.write('\n')

    # write config file to run summary

    config_msg = '''
####################################################################################
# CONFIG FILE
####################################################################################

'''
    print(config_msg)

    f.write(config_msg)

    with open(cfg_file, 'r') as cfg:

        for line in cfg:

            try:
                if line.split(' ')[1][0] == '=' and line[0] != '#':
                    setting = line.split(' ')[0]
                    setting_value = globals()[setting]
                    print(setting + ': ' + str(setting_value))
                    f.write(setting + ': ' + str(setting_value) + '\n')

                else:
                    print(line.strip())
                    f.write(line)
            except:
                print(line.strip())
                f.write(line)
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
                print('Unexpected filename structure. Exiting...')
                raise SystemExit

    if len(R1_files) != len(R2_files) and paired:
        print('Unpaired FASTQ files exist! Check input files.')
        raise SystemExit

    R1_files.sort()
    R2_files.sort()

    return R1_files, R2_files

def generate_sample(R1_files, R2_files, dna_only, ab_only, sample_name, temp_dir, ab_dir):
    # store filenames in TapestriTube objects

    tubes = []        # list of TapestriTube objects

    # for experiments with only DNA panels
    if dna_only:

        for i in range(0, len(R1_files)):
            # assign filenames to sample types
            # note: using alphabetization pattern which may not exist in future

            tube_num = i + 1

            # dna panel

            panel_r1 = R1_files[i]
            panel_r2 = R2_files[i]

            # set file locations

            panel_r1_temp = temp_dir + panel_r1.split('.fastq.gz')[0].split('/')[-1] + '.temp.fastq'
            panel_r2_temp = temp_dir + panel_r2.split('.fastq.gz')[0].split('/')[-1] + '.temp.fastq'

            # antibodies - set filenames to empty string

            ab_r1 = ''
            ab_r2 = ''

            ab_r1_temp = ''
            ab_r2_temp = ''

            ab_reads = ''
            umi_counts = ''

            tubes.append(TapestriTube(tube_num,
                                      panel_r1,
                                      panel_r2,
                                      panel_r1_temp,
                                      panel_r2_temp,
                                      ab_r1,
                                      ab_r2,
                                      ab_r1_temp,
                                      ab_r2_temp,
                                      ab_reads,
                                      umi_counts))

    # for experiments with only Ab panels
    elif ab_only:

        for i in range(0, len(R1_files)):
            # assign filenames to sample types
            # note: using alphabetization pattern which may not exist in future

            tube_num = i + 1

            # dna panel - set files to empty strings

            panel_r1 = ''
            panel_r2 = ''

            # set file locations

            panel_r1_temp = ''
            panel_r2_temp = ''

            # antibodies

            ab_r1 = R1_files[i]
            ab_r2 = R2_files[i]

            ab_r1_temp = temp_dir + ab_r1.split('.fastq.gz')[0].split('/')[-1] + '.temp.fastq'
            ab_r2_temp = temp_dir + ab_r2.split('.fastq.gz')[0].split('/')[-1] + '.temp.fastq'

            ab_reads = ab_dir + sample_name + '-' + str(tube_num) + '.ab_reads.tsv'
            umi_counts = ab_dir + sample_name + '-' + str(tube_num) + '.umi_counts.tsv'

            tubes.append(TapestriTube(tube_num,
                                      panel_r1,
                                      panel_r2,
                                      panel_r1_temp,
                                      panel_r2_temp,
                                      ab_r1,
                                      ab_r2,
                                      ab_r1_temp,
                                      ab_r2_temp,
                                      ab_reads,
                                      umi_counts))

    # for experiments with antibody and panel data
    else:

        for i in range(0, len(R1_files) / 2):
            # assign filenames to sample types
            # note: using alphabetization pattern which may not exist in future

            tube_num = i + 1

            # dna panel

            panel_r1 = R1_files[i + len(R1_files) / 2]
            panel_r2 = R2_files[i + len(R1_files) / 2]

            # set file locations

            panel_r1_temp = temp_dir + panel_r1.split('.fastq.gz')[0].split('/')[-1] + '.temp.fastq'
            panel_r2_temp = temp_dir + panel_r2.split('.fastq.gz')[0].split('/')[-1] + '.temp.fastq'

            # antibodies

            ab_r1 = R1_files[i]
            ab_r2 = R2_files[i]

            ab_r1_temp = temp_dir + ab_r1.split('.fastq.gz')[0].split('/')[-1] + '.temp.fastq'
            ab_r2_temp = temp_dir + ab_r2.split('.fastq.gz')[0].split('/')[-1] + '.temp.fastq'

            ab_reads = ab_dir + sample_name + '-' + str(tube_num) + '.ab_reads.tsv'
            umi_counts = ab_dir + sample_name + '-' + str(tube_num) + '.umi_counts.tsv'

            tubes.append(TapestriTube(tube_num,
                                      panel_r1,
                                      panel_r2,
                                      panel_r1_temp,
                                      panel_r2_temp,
                                      ab_r1,
                                      ab_r2,
                                      ab_r1_temp,
                                      ab_r2_temp,
                                      ab_reads,
                                      umi_counts))

    return tubes
