# -*- coding: utf-8 -*-
"""
mission bio single-cell pipeline code
written by ben 11.25.2018

"""

# modules
from __future__ import division
import os
import os.path
import collections
import csv
from itertools import product, combinations, izip, groupby
from multiprocessing import Process, Queue
import json
from collections import Counter
from itertools import izip
import subprocess
import sys

# plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.family'] = 'Helvetica'
plt.rcParams['axes.unicode_minus'] = False

# add the modified umi_tools directory to python path
# TODO copy files into repo
sys.path.insert(0, '/usr/local/lib/python2.7/dist-packages/umi_tools-0.5.3-py2.7-linux-x86_64.egg/umi_tools/')
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
                           cell_barcode_structure,
                           mb_barcodes,
                           bar_ind_1,
                           bar_ind_2,
                           sample_type):
        # filter r1 and r2 files to only keep reads with correct barcode structure in r1

        assert sample_type == 'ab' or sample_type == 'panel', 'Sample type must be panel or ab!'

        # set filenames according to sample type (panel or ab)
        if sample_type == 'panel':
            r1_in = self.panel_r1
            r2_in = self.panel_r2

            r1_out = self.panel_r1_temp
            r2_out = self.panel_r2_temp

            barcode_json = self.panel_barcodes

        else:
            r1_in = self.ab_r1
            r2_in = self.ab_r2

            r1_out = self.ab_r1_temp
            r2_out = self.ab_r2_temp

            barcode_json = self.ab_barcodes

        filter_cmd = 'python3 /usr/local/bin/cutadapt -g %s -m 100 -j 12 --interleaved' \
                     ' --discard-untrimmed --no-trim --pair-filter=both -e 0.1 %s %s --quiet' \
                     % (cell_barcode_structure,
                        r1_in,
                        r2_in)

        filter_process = subprocess.Popen(filter_cmd, stdout=subprocess.PIPE, shell=True)

        # create output fastq files for writing
        out_1 = open(r1_out, 'w')
        out_2 = open(r2_out, 'w')

        # create dict for storing barcode information (will be exported to json)
        barcodes = {}

        trim = len(cell_barcode_structure) - 1  # base indices to trim from read 1
        bar_count = 0

        # iterate through all reads
        for line in filter_process.stdout:

            # R1
            header_1 = line.strip()                           # header
            id_1 = header_1.split(' ')[0]                     # read id
            seq_1 = filter_process.stdout.next().strip()      # sequence string
            filter_process.stdout.next()                      # + sign
            qual_1 = filter_process.stdout.next().strip()     # quality string

            # R2
            header_2 = filter_process.stdout.next().strip()
            id_2 = header_2.split(' ')[0]
            seq_2 = filter_process.stdout.next().strip()
            filter_process.stdout.next()
            qual_2 = filter_process.stdout.next().strip()

            assert id_1 == id_2, 'read ids from paired-end files must match!'

            # extract barcodes and check that they are a valid MB barcode
            check = check_seq(seq_1, bar_ind_1, bar_ind_2, mb_barcodes)

            if check == 'fail':
                continue

            else:
                barcode = check[0] + check[1] + '-' + str(self.sample_num)
                bar_count += 1

                # add barcode to dict
                try:
                    barcodes[barcode].append(id_1)

                except KeyError:
                    barcodes[barcode] = [id_1]

                # add barcode to headers
                header_1 = header_1.split(' ')[0] + '_' + barcode + '_' + header_1.split(' ')[1]
                header_2 = header_2.split(' ')[0] + '_' + barcode + '_' + header_2.split(' ')[1]

                # save valid fastq records
                out_1.write('%s\n%s\n+\n%s\n' % (header_1, seq_1[trim:], qual_1[trim:]))
                out_2.write('%s\n%s\n+\n%s\n' % (header_2, seq_2, qual_2))

                # # print counter
                # if bar_count % 1e6 == 0:
                #     print '%d valid read-pairs extracted and written to new fastq.' % bar_count

        # close fastq and save barcode json
        out_1.close()
        out_2.close()
        json_export(barcodes, barcode_json)

        print '%d read-pairs written to %s and %s.\n' % (bar_count, r1_out, r2_out)

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

            # read_id = line.split(' ')[0]    # read id

            cell_barcode = line.split('_')[1]  # extract corrected barcode from header

            # # check that this ab data has a good cell barcode
            # if read_id not in read_id_dict:
            #     # print 'invalid cell barcode'
            #     break

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
            passed_ab_reads[passed_count]['cell barcode'] = cell_barcode#read_id_dict[read_id]
            passed_ab_reads[passed_count]['ab description'] = barcode_descriptions[bar]
            passed_ab_reads[passed_count]['raw umi'] = umi

        # write passed ab reads to tsv file
        with(open(self.ab_reads, 'w')) as f:
            for ab in passed_ab_reads:
                f.write(passed_ab_reads[ab]['cell barcode'] + '\t')
                f.write(passed_ab_reads[ab]['ab description'] + '\t')
                f.write(passed_ab_reads[ab]['raw umi'] + '\n')

    def fastq_split_by_cell(self, min_reads, by_cell_fastq_dir):
        # split fastq files into separate files for each valid cell

        # load panel barcode json from file
        barcodes = json_import(self.panel_barcodes)

        # add counts for barcodes
        for b in barcodes:
            barcodes[b] = len(barcodes[b])

        # dict of barcodes observed
        good_barcodes = {}

        r1 = open(self.panel_r1_temp, 'r')
        r2 = open(self.panel_r2_temp, 'r')

        reads_written = 0

        # iterate through lines of the FASTQ files together

        for line1, line2 in izip(r1, r2):
 
            barcode = line1.split('_')[1]   # extract corrected barcode from header

            if barcodes[barcode] >= min_reads:

                if barcode in good_barcodes:
                    i = len(good_barcodes[barcode]) + 1
                    good_barcodes[barcode][i] = ''.join([x for x in
                                                        [line1, r1.next(), r1.next(), r1.next(),
                                                         line2, r2.next(), r2.next(), r2.next()]])

                else:
                    good_barcodes[barcode] = {}
                    good_barcodes[barcode][0] = ''.join([x for x in
                                                        [line1, r1.next(), r1.next(), r1.next(),
                                                         line2, r2.next(), r2.next(), r2.next()]])

                reads_written += 1

                # if reads_written % 1e6 == 0:
                #     print '%d reads saved.' % reads_written

            else:
                r1.next(); r1.next(); r1.next();
                r2.next(); r2.next(); r2.next();

        r1.close()
        r2.close()

        for barcode in good_barcodes:
            f = open('%s%s.fastq' % (by_cell_fastq_dir, barcode), 'w')
            for i in good_barcodes[barcode]:
                f.write(good_barcodes[barcode][i])
            f.close()

        print 'Sample %d: %d reads saved to %d cell fastq files.' % (self.sample_num, reads_written, len(good_barcodes))

class SingleCell(object):
    # class for storing metadata for each single cell file

    def __init__(self, cell_barcode, fastq_dir, bam_dir, vcf_dir, interval_dir):
        # initialize object by generating filenames

        self.cell_barcode = cell_barcode                    # cell barcode
        self.fastq = fastq_dir + cell_barcode + '.fastq'    # fastq file

        assert self.fastq, 'fastq file does not exist'

        self.bam = bam_dir + cell_barcode + '.bam'          # bam file
        self.bai = bam_dir + cell_barcode + '.bai'          # bam file index
        self.vcf = vcf_dir + cell_barcode + '.g.vcf'        # gvcf file
        self.interval_aln = interval_dir + cell_barcode + '.tsv' # interval alignment file

        self.valid = False      # marker for valid cells
        self.alignments = {}    # alignment counts for each interval

    def align_sample(self, bt2_ref):
        # align the panel to the bowtie2 human index and generate sorted bam file

        align_cmd = '/usr/local/bin/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -x %s --mm --interleaved %s' \
                    ' --rg-id %s --rg SM:%s --rg PL:ILLUMINA --rg CN:UCSF --quiet' \
                    ' | samtools view -b' \
                    ' | samtools sort -o %s' \
                    % (bt2_ref,
                       self.fastq,
                       self.cell_barcode,
                       self.cell_barcode,
                       self.bam)

        process = subprocess.Popen(align_cmd, shell=True)

        return process

    def index_bam(self):
        # index all bam files using samtools

        index_cmd = 'samtools index %s %s' \
                    % (self.bam,
                       self.bai)

        process = subprocess.Popen(index_cmd, shell=True)

        return process

    def interval_alignments(self, interval_file):
        # get number of reads aligned in intervals

        int_aln_cmd = 'gatk CollectReadCounts -I %s -O %s -L %s --verbosity ERROR' \
                      ' --interval-merging-rule OVERLAPPING_ONLY --format TSV' \
                       % (self.bam,
                          self.interval_aln,
                          interval_file)

        process = subprocess.Popen(int_aln_cmd, shell=True)

        return process

    def call_variants(self, fasta, interval_file, dbsnp_file):
        # call variants using gatk

        variants_cmd = 'gatk HaplotypeCaller -R %s -I %s -O %s -ERC BP_RESOLUTION -L %s -D %s --verbosity ERROR' \
                        ' --native-pair-hmm-threads 1' \
                       % (fasta,
                          self.bam,
                          self.vcf,
                          interval_file,
                          dbsnp_file)

        process = subprocess.Popen(variants_cmd, shell=True)

        return process

    def cell_caller(self, interval_dict, min_coverage, min_fraction):
        # determine if a barcode is a valid cell using alignment counts

        aln_dict = self.alignments_from_tsv(self.interval_aln, interval_dict)
        self.alignments = aln_dict

        # mark valid cells based on coverage uniformity
        if len([a for a in aln_dict if aln_dict[a] >= min_coverage])/len(aln_dict) >= min_fraction:
            self.valid = True

    @staticmethod
    def alignments_from_tsv(tsv_file, interval_dict):
        # extracts alignment info from tsv file and imports into a dictionary

        aln_dict = {}

        with open(tsv_file) as f:
            for line in f:

                # skip header lines
                if line[0] == '@' or line[:6] == 'CONTIG':
                    continue

                else:
                    interval = line.strip().split('\t')[:3]
                    # need to convert between 0 and 1 based coordinates
                    interval = interval[0] + '_' + str(int(interval[1]) - 1) + '_' + str(int(interval[2]))
                    aln_dict[interval_dict[interval]] = int(line.strip().split('\t')[3])

        return aln_dict

    @staticmethod
    def generate_alignments_tsv(cells, interval_dict, out_file):
        # create tsv with number of alignments for each amplicon across cells

        amplicons = interval_dict.values()
        amplicons.sort()

        out = open(out_file, 'w')
        out.write('barcode\t' + '\t'.join(amplicons) + '\n')

        for c in cells:
            if c.valid:
                out.write(c.cell_barcode + '\t' + '\t'.join([str(c.alignments[a]) for a in amplicons]) + '\n')


        out.close()

    @staticmethod
    def combine_gvcfs(cells, id, fasta, interval_file, dbsnp_file, merged_gvcf_dir, genotyping_dir, multi_sample = False):
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
                       ' --max-alternate-alleles 2' \
                       % (fasta,
                          merged_gvcf_file,
                          geno_gvcf_file,
                          dbsnp_file)

        process = subprocess.Popen(genotype_cmd, shell=True)

        return process

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

def count_umis(ab_reads_file, umi_counts_file, n_procs=32):
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








