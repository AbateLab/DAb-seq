# -*- coding: utf-8 -*-
"""
mission bio single-cell pipeline code
written by ben 11.25.2018

"""

# modules
from __future__ import division
import os
import os.path
import csv
from itertools import product, combinations, groupby
from multiprocessing import Process, Queue
import json
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

class TapestriSample(object):
    # class for storing metadata for each tapestri sample (one tube, run, etc...)

    def __init__(self,
                 sample_num,
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
                 ab_reads):

        self.sample_num = sample_num            # number identifying sample (or tube)

        self.panel_r1 = panel_r1                # panel R1 fastq path
        self.panel_r2 = panel_r2                # panel R2 fastq path
        self.panel_r1_temp = panel_r1_temp      # temp panel R1 fastq path
        self.panel_r2_temp = panel_r2_temp      # temp panel R2 fastq path

        self.ab_r1 = ab_r1                      # ab R1 fastq path
        self.ab_r2 = ab_r2                      # ab R2 fastq path
        self.ab_r1_temp = ab_r1_temp            # temp ab R1 fastq path
        self.ab_r2_temp = ab_r2_temp            # temp ab R2 fastq path

        self.panel_barcodes = panel_barcodes    # file of cell barcodes for this sample (panel)
        self.ab_barcodes = ab_barcodes          # file of cell barcodes for this sample (abs)

        self.ab_reads = ab_reads                # file containing filtered ab reads

    def filter_valid_reads(self,
                           r1_start,
                           mb_barcodes,
                           bar_ind_1,
                           bar_ind_2,
                           sample_type):
        # filter r1 files to only keep reads with correct barcode structure in r1

        assert sample_type == 'ab' or sample_type == 'panel', 'Sample type must be panel or ab!'

        # set filenames according to sample type (panel or ab)
        if sample_type == 'panel':
            r1_in = self.panel_r1
            barcode_json = self.panel_barcodes

        else:
            r1_in = self.ab_r1
            barcode_json = self.ab_barcodes

        cmd = 'python3 /usr/local/bin/cutadapt' \
              ' -a r1_start=%s' \
              ' -j 16 -O 6 -e 0.2 %s --quiet' \
              % (r1_start,
                 r1_in)

        trim_process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

        read_id_dict = {}  # dict for storing read ids

        # iterate through all Read 1 records
        for line in trim_process.stdout:

            # R1
            header_1 = line.strip()  # header
            id_1 = header_1.split(' ')[0][1:]  # read id
            seq_1 = trim_process.stdout.next().strip()  # sequence string
            trim_process.stdout.next()  # +
            trim_process.stdout.next()  # qual

            # find barcodes and check that they are a valid MB barcode
            check = check_seq(seq_1, bar_ind_1, bar_ind_2, mb_barcodes)

            if check == 'fail':
                continue

            else:
                barcode = check[0] + check[1] + '-' + str(self.sample_num)
                read_id_dict[id_1] = barcode

        # export barcodes to json file
        json_export(read_id_dict, barcode_json)

    def barcode_reads(self,
                      r1_start,
                      r1_end,
                      r2_end,
                      r1_min_len,
                      r2_min_len,
                      sample_type):
        # for valid reads, add barcode header to fastq file and trim

        assert sample_type == 'ab' or sample_type == 'panel', 'Sample type must be panel or ab!'

        # set filenames according to sample type (panel or ab)
        if sample_type == 'panel':
            r1_in = self.panel_r1
            r2_in = self.panel_r2

            r1_out = open(self.panel_r1_temp, 'w')
            r2_out = open(self.panel_r2_temp, 'w')

            barcode_json = self.panel_barcodes

            # TODO make read 1 cutting parameter (-u 51) intelligent

            # cutadapt cmd for panel (trim read before finding adapter)
            cmd = 'python3 /usr/local/bin/cutadapt' \
                  ' -a %s' \
                  ' -A %s' \
                  ' --interleaved -j 16 -u 51 -U 5 -n 3 -O 8 -e 0.2 %s %s --quiet' \
                  % (r1_end,
                     r2_end,
                     r1_in,
                     r2_in)

        elif sample_type == 'ab':
            r1_in = self.ab_r1
            r2_in = self.ab_r2

            r1_out = open(self.ab_r1_temp, 'w')
            r2_out = open(self.ab_r2_temp, 'w')

            barcode_json = self.ab_barcodes

            # cutadapt cmd for abs
            cmd = 'python3 /usr/local/bin/cutadapt' \
                  ' -g %s' \
                  ' -a %s' \
                  ' -A %s' \
                  ' --interleaved -j 16 -n 3 -O 8 -e 0.2 %s %s --quiet' \
                  % (r1_start,
                     r1_end,
                     r2_end,
                     r1_in,
                     r2_in)

        read_id_dict = json_import(barcode_json)

        trim_process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

        bar_count = 0  # total count of all barcodes

        # iterate through all Read 1 records
        for line in trim_process.stdout:

            # R1
            header_1 = line.strip()
            id_1 = header_1.split(' ')[0][1:]
            seq_1 = trim_process.stdout.next().strip()
            trim_process.stdout.next()
            qual_1 = trim_process.stdout.next().strip()

            # R2
            header_2 = trim_process.stdout.next().strip()
            id_2 = header_2.split(' ')[0][1:]
            seq_2 = trim_process.stdout.next().strip()
            trim_process.stdout.next()
            qual_2 = trim_process.stdout.next().strip()

            assert id_1 == id_2, 'Read IDs do not match! Check input FASTQ files.'

            try:
                cell_barcode = read_id_dict[id_1]
            except KeyError:
                continue

            if len(seq_1) < r1_min_len or len(seq_2) < r2_min_len:
                continue

            else:
                # add barcode to header
                id = '@' + id_1 + '_' + cell_barcode
                header_1 = id
                header_2 = id

                # write to output fastq files
                r1_out.write('%s\n%s\n+\n%s\n' % (header_1, seq_1, qual_1))
                r2_out.write('%s\n%s\n+\n%s\n' % (header_2, seq_2, qual_2))

                bar_count += 1

        print '%d total valid trimmed pairs saved to file.' % bar_count

        r1_out.close()
        r2_out.close()

    def process_abs(self,
                    ab_barcodes,
                    barcode_descriptions,
                    ab_handles,
                    ab_bar_coord,
                    ab_umi_coord,
                    min_umi_qual):
        # extract ab barcodes and umis from raw ab reads

        # create dict for storing passed ab reads
        passed_ab_reads = {}
        passed_count = 0

        # use cutadapt to select reads with correct structure
        ab_cmd = 'python3 /usr/local/bin/cutadapt -j 24 %s -O 12 -e 0.2 -n 2 %s --quiet ' \
                 '--discard-untrimmed' % (ab_handles, self.ab_r2_temp)

        ab_process = subprocess.Popen(ab_cmd, stdout=subprocess.PIPE, shell=True)

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

            # if a read passes all filters, add it to a dictionary
            passed_count += 1
            passed_ab_reads[passed_count] = {}
            passed_ab_reads[passed_count]['cell barcode'] = cell_barcode
            passed_ab_reads[passed_count]['ab description'] = barcode_descriptions[bar]
            passed_ab_reads[passed_count]['raw umi'] = umi

        # write passed ab reads to tsv file
        with(open(self.ab_reads, 'w')) as f:
            for ab in passed_ab_reads:
                f.write(passed_ab_reads[ab]['cell barcode'] + '\t')
                f.write(passed_ab_reads[ab]['ab description'] + '\t')
                f.write(passed_ab_reads[ab]['raw umi'] + '\n')

class SingleCell(object):
    # class for storing metadata for each single cell file

    def __init__(self, cell_barcode, fastq_dir, bam_dir, vcf_dir):
        # initialize object by generating filenames

        self.cell_barcode = cell_barcode                    # cell barcode
        self.fastq = fastq_dir + cell_barcode + '.fastq'    # fastq file

        assert self.fastq, 'fastq file does not exist'

        self.bam = bam_dir + cell_barcode + '.bam'          # bam file
        self.bai = bam_dir + cell_barcode + '.bai'          # bam file index
        self.vcf = vcf_dir + cell_barcode + '.g.vcf'        # gvcf file

        self.valid = False      # marker for valid cells
        self.alignments = {}    # alignment counts for each interval

    def align_and_index(self, bt2_ref):

        # align the panel to the bowtie2 human index and generate sorted bam file
        # requirements: read mapped, mapq >= 30, primary alignment
        align_cmd = '/usr/local/bin/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -x %s --mm --interleaved %s' \
                    ' --rg-id %s --rg SM:%s --rg PL:ILLUMINA --rg CN:UCSF --quiet' \
                    ' | samtools view -b -q 30 -F 4 -F 0X0100' \
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

    def call_variants(self, fasta, interval_file, dbsnp_file):
        # call variants using gatk

        variants_cmd = 'gatk HaplotypeCaller -R %s -I %s -O %s -ERC BP_RESOLUTION -L %s -D %s --verbosity ERROR' \
                        ' --native-pair-hmm-threads 1 --max-reads-per-alignment-start 0' \
                       % (fasta,
                          self.bam,
                          self.vcf,
                          interval_file,
                          dbsnp_file)

        # process = subprocess.Popen(variants_cmd, shell=True)
        subprocess.call(variants_cmd, shell=True)

    @staticmethod
    def combine_gvcfs(cells, id, fasta, interval_file, dbsnp_file, merged_gvcf_dir, genotyping_dir, multi_sample=False):
        # combine single-cell gvcfs

        gvcf_list = genotyping_dir + 'gvcfs.txt'

        # if cells is a list of filenames
        if multi_sample:
            out_file = id
            with open(genotyping_dir + 'gvcfs.txt', 'w') as f:
                for c in cells:
                    f.write('--variant ' + c + '\n')

        # if cells is a SingleCell instance
        else:
            out_file = merged_gvcf_dir + str(id) + '_merged.g.vcf'
            with open(genotyping_dir + 'gvcfs.txt', 'w') as f:
                for c in cells:
                    f.write('--variant ' + c.vcf + '\n')

        combine_cmd = 'gatk CombineGVCFs -R %s --arguments_file %s -O %s -L %s -D %s' \
                       % (fasta,
                          gvcf_list,
                          out_file,
                          interval_file,
                          dbsnp_file)

        process = subprocess.Popen(combine_cmd, shell=True)

        return process

    @staticmethod
    def genotype_gvcfs(fasta, dbsnp_file, merged_gvcf_file, geno_gvcf_file):
        # genotype the multisample merged gvcf file

        genotype_cmd = 'gatk GenotypeGVCFs -R %s -V %s -O %s -D %s' \
                       ' --standard-min-confidence-threshold-for-calling 30' \
                       ' --max-alternate-alleles 2 --includeNonVariantSites' \
                       % (fasta,
                          merged_gvcf_file,
                          geno_gvcf_file,
                          dbsnp_file)

        process = subprocess.Popen(genotype_cmd, shell=True)

        return process

def count_alignments(r1_files, amplicon_file, fasta_file, tsv, dir):
    # align and count r1 reads for all barcodes, and save to tsv file

    # get fasta file from human genome for this interval
    insert_fasta = dir + 'amplicons.fasta'
    subprocess.call('bedtools getfasta -fi %s -bed %s -fo %s -name' % (fasta_file, amplicon_file, insert_fasta), shell=True)

    # build bt2 index for this fasta
    insert_bt2 = dir + 'inserts'
    subprocess.call('bowtie2-build %s %s' % (insert_fasta, insert_bt2), shell=True)

    # extract names of reference sequences
    refs = []
    get_refs = subprocess.Popen('bowtie2-inspect -n %s' % insert_bt2, stdout=subprocess.PIPE, shell=True)
    for line in get_refs.stdout:
        refs.append(line.strip())
    refs.sort()
    refs_dict = dict(zip(refs, [0] * len(refs)))

    # amplicon dict (key: cell barcode, value: list of amplicon counts)
    amplicons = {}

    # align reads with bowtie2
    bt2_input = ' -U '.join(r1_files)
    bt2_cmd = 'bowtie2 -p 24 -x %s -U %s' % (insert_bt2, bt2_input)
    bt2_align = subprocess.Popen(bt2_cmd, stdout=subprocess.PIPE, shell=True)

    # iterate through all reads
    for line in bt2_align.stdout:

        if line[0] == '@':  # ignore header lines
            continue

        # parse sam records
        record = line.split('\t')
        flag = int(record[1])
        unmapped = int(bin(flag)[-3])  # mapping bit (1: not mapped; 0: mapped)
        mapq = int(record[4])

        if flag >= 256:
            secondary = int(bin(flag)[-9])  # secondary bit (1: non-primary; 0: primary)
        else:
            secondary = 0

        # discard unmapped reads
        if unmapped == 1:
            continue

        # discard secondary alignments
        if secondary == 1:
            continue

        # discard reads with low mapping quality
        min_mapq = 30
        if mapq < min_mapq:
            continue

        # if read passes all filters, extract barcode
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

def json_import(filename):
    # imports json data into python. If json file does not exist, returns empty {}

    if not os.path.exists(filename):
        json_obj = {}
    else:
        with open(filename) as f:
            json_obj = json.load(f)

    return json_obj

def json_export(json_obj, filename, overwrite=True, update=False):
    # exports a json object to file. If overwrite is off, writing will fail. Can also update an existing json

    if os.path.exists(filename) and not overwrite:
        print 'File exists. Will not overwrite. Exiting...'
        raise SystemExit

    elif update:
        if not os.path.exists(filename):
            old_json = {}
        else:
            old_json = json_import(filename)

        json_obj.update(old_json)

    with open(filename, 'w') as out:
        json.dump(json_obj, out)

def load_barcodes(barcode_file, max_dist, check_barcodes=True):
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
            print 'Barcodes must all be same length! Exiting...'
            raise SystemExit

        # check pairwise hamming distances
        dist_req = 2 * max_dist + 1         # for this max_dist, need this distance between all barcodes
        pairs = list(combinations(barcodes.keys(), 2))
        for pair in pairs:
            if hd(pair[0], pair[1]) < dist_req:
                print 'Error: The edit distance between barcodes %s and %s is less than %d.\n' \
                      'An error correction of %d bases will not work.' % (pair[0], pair[1], dist_req, max_dist)

    return barcodes

def generate_hamming_dict(barcodes):
    # return dictionary of strings within 1 hamming distance of barcodes

    for barcode in barcodes:
        hd_1 = sorted(hamming_circle(barcode, 1))
        barcodes[barcode] = dict(zip(hd_1, [1] * len(hd_1)))

    return barcodes

# from https://codereview.stackexchange.com/a/88919
def hamming_circle(s, n):
    # generate strings over alphabet whose hamming distance from s is exactly n

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

def count_umis(ab_reads_file, umi_counts_file, n_procs=24):
    # count umis using selected clustering methods

    def extract_umis(ab_reads_file):
        # extracts lists of umis from ab reads tsv file and groups by cell and ab tag

        # store ab reads in list of lists
        # col 1: cell_barcode, col 2: ab_barcode_description, col 3: raw_umi
        ab_reads = []
        with(open(ab_reads_file, 'r')) as f:
            for line in f:
                ab_reads.append(line.strip().split('\t'))

        # concatenate cell barcode and ab barcode
        ab_reads = [['+'.join(r[:2]), r[2]] for r in ab_reads]
        ab_reads.sort()   # sort list by first element

        group_ids = []  # cell barcode + ab barcode defining a umi group
        umi_groups = []  # list of list of umis for each group

        # group by ab tag barcode (first element in list)
        keyfunc = lambda x: x[0]

        for k, g in groupby(ab_reads, keyfunc):
            group_ids += [k]  # add ab tag
            umi_groups += [[i[1] for i in list(g)]]  # add umi group

        return group_ids, umi_groups

    def cluster_umis(umi_dict, method='all'):
        # clusters the umis using the specified method (or all)
        # uses functions from umi-tools paper (Genome Research, 2017)

        # split umi dict into umis (keys) and counts
        umis = umi_dict.keys()
        counts = umi_dict

        # set up UMIClusterer functor with parameters specific to specified method
        # choose method = 'all' for all available methods
        # otherwise provide methods as a list of methods
        processor = network.UMIClusterer()  # initialize UMIclusterer

        # cluster the umis
        clusters = processor(
            umis,
            counts,
            threshold=1,
            cluster_method=method)

        return clusters

    def mp_umi_cluster(group_ids, umis, n_procs, method='all'):
        # multiprocess implementation of umi-tools clustering algorithms
        # assigns one umi group to one process, processing n_procs groups in parallel

        # make dicts of all umi groups and an index
        umi_groups = dict(zip(group_ids, umis))
        umi_groups_idx = dict(zip(range(len(group_ids)), group_ids))

        def worker(out_q, proc_id):
            # the function that sends umi groups to the umi-tools clustering functions

            for i in umi_groups_idx.keys():
                if i % n_procs == proc_id:  # delegate entries to each process
                    umi_counter = Counter(umi_groups[umi_groups_idx[i]])
                    clusters = cluster_umis(umi_counter, method=method)
                    clusters_by_group = {umi_groups_idx[i]: clusters}

                    out_q.put(clusters_by_group)

                else:
                    continue

            sys.stdout.flush()

        # each process will get random subset of umi groups
        out_q = Queue()
        procs = []

        for i in range(n_procs):
            p = Process(
                target=worker,
                args=(out_q, i))
            procs.append(p)
            p.start()

        # collect all results into a single output dict
        out_dict = {}

        while True:
            out_dict.update(out_q.get())
            if len(out_dict) == len(group_ids):
                break

        # wait for all worker processes to finish
        for p in procs:
            p.join()

        return out_dict

    def to_counts(groups_by_method):
        # counts umis for each method in umi dictionary

        counts_by_method = {}

        for barcode in groups_by_method:
            counts_by_method[barcode] = {}

            for method in groups_by_method[barcode]:
                counts_by_method[barcode][method] = len(groups_by_method[barcode][method])

        return counts_by_method

    # extract ab reads from tsv file
    group_ids, umi_groups = extract_umis(ab_reads_file)
    raw_counts_dict = dict(zip(group_ids, [len(l) for l in umi_groups]))

    # number of threads should be less than or equal to number of umi groups
    if n_procs > len(group_ids):
        n_procs = len(group_ids)

    # if no UMI groups found, quit program
    if len(group_ids) == 0:
        print 'No UMI groups found! Now exiting...\n'
        raise SystemExit

    print 'Found %s UMI groups. Now clustering using %d threads.' % (len(group_ids), n_procs)

    # send groups to multiprocessor handler and clustering functions
    # groups_by_method = {barcode1: {'adjacency': [[umi1],[umi1, umi2]], 'directional': ...}}

    groups_by_method = mp_umi_cluster(group_ids, umi_groups, n_procs, method='all')

    # count number of umis in each barcode for each method
    # counts_by_method = {barcode1: {'adjacency': 2, 'directional': ...}}
    counts_by_method = to_counts(groups_by_method)

    # save the counts to tsv file
    count_methods = ['unique', 'percentile', 'cluster', 'adjacency']
    with(open(umi_counts_file, 'w')) as f:
        f.write('cell_barcode\tab_description\traw\t' + '\t'.join(count_methods) + '\n')
        for group in counts_by_method:
            f.write('\t'.join(group.split('+')) + '\t' + str(raw_counts_dict[group]) + '\t')
            f.write('\t'.join([str(counts_by_method[group][m]) for m in count_methods]) + '\n')

    print 'All UMIs grouped and saved to %s.\n' % umi_counts_file

def bcftools_annotate(annotations_vcf, input_vcf, column_info, output_vcf):
    # uses bcftools to annotate a vcf file with annotation information from another
    # input and output vcf should both be uncompressed vcf

    # bgzip and index the input vcf
    subprocess.call('bgzip -@ 16 %s' % input_vcf, shell=True)
    subprocess.call('tabix %s' % input_vcf + '.gz', shell=True)

    # use bcftools to annotate the input
    bcf_cmd = 'bcftools annotate -a %s %s %s > %s' % (annotations_vcf,
                                                      column_info,
                                                      input_vcf + '.gz',
                                                      output_vcf)
    subprocess.call(bcf_cmd, shell=True)


def vcf_to_tables(vcf_file, genotype_file, variants_tsv):
    # parses a vcf file into a series of tables

    # load vcf file into numpy array
    # include annotation info from snpeff
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
    # GT: genotype (0: WT, 1: HET, 2: HOM, 3: no call)
    # DP: total read depth
    # GQ: genotype quality
    # AD: alt allele depth
    # RD: ref allele depth

    GT = np.sum(vcf['calldata/GT'], axis=2)
    GT[GT == -2] = 3
    DP = np.stack(vcf['calldata/DP'], axis=0)
    GQ = np.stack(vcf['calldata/GQ'], axis=0)
    AD = np.stack(vcf['calldata/AD'][:, :, 1], axis=0)
    RD = np.stack(vcf['calldata/AD'][:, :, 0], axis=0)

    # create variant names
    names = [vcf['variants/ANN_Gene_Name'][i] +
             ':' + vcf['variants/CHROM'][i] +
             ':' + str(vcf['variants/POS'][i]) +
             ':' + vcf['variants/REF'][i] +
             '/' + vcf['variants/ALT'][:, 0][i]
             for i in range((vcf['variants/REF'].shape[0]))]

    # assemble and save variant annotations to file
    variants_table = pd.DataFrame(data=names, columns=['Name'])

    # cosmic id
    variants_table['COSMIC_ID'] = vcf['variants/ID']

    # snpeff columns
    ANN_columns = [c for c in list(vcf) if '/ANN' in c]
    for ann in ANN_columns:
        variants_table['SnpEff_' + ann.split('/ANN_')[1]] = vcf[ann]

    # clinvar columns
    CLN_columns = [c for c in list(vcf) if '/CLN' in c]
    for cln in CLN_columns:
        variants_table['ClinVar_' + cln.split('/')[1]] = vcf[cln]

    variants_table.to_csv(path_or_buf=variants_tsv, sep='\t', index=False)

    # encode variant names and cell barcodes
    names = [n.encode('utf8') for n in names]
    barcodes = [b.encode('utf8') for b in vcf['samples']]

    # save genotyping information to compressed hdf5 file
    with h5py.File(genotype_file, 'w') as f:
        f.create_dataset('GT', data=GT, dtype='i1', compression='gzip')
        f.create_dataset('GQ', data=GQ, dtype='i1', compression='gzip')
        f.create_dataset('DP', data=DP, dtype='i2', compression='gzip')
        f.create_dataset('AD', data=AD, dtype='i2', compression='gzip')
        f.create_dataset('RD', data=RD, dtype='i2', compression='gzip')
        f.create_dataset('VARIANTS', data=names, compression='gzip')
        f.create_dataset('CELL_BARCODES', data=barcodes, compression='gzip')







