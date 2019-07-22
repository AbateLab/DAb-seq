'''
basespace fastq downloader
ben demaree 7.9.2019
'''

import os
import sys
import json
import math
import subprocess
import re
from urllib2 import Request, urlopen, URLError
import argparse
import shutil

base_url = "https://api.basespace.illumina.com/"

def restrequest(rawrequest):
    request = Request(rawrequest)

    try:
        response = urlopen(request)
        json_string = response.read()
        json_obj = json.loads(json_string)

    except URLError, e:
        print 'Got an error code:', e
        sys.exit()

    return json_obj

def download_files(href, access_token, output_dir):
    # downloads files for a given sample

    downloaded_files = []

    if not os.path.isdir(output_dir):
        raise OSError(output_dir + " is not a directory")

    request = "%s?access_token=%s&limit=1" % (href, access_token)
    json_obj = restrequest(request)

    total_file_count = json_obj["Paging"]["TotalCount"]

    num_offsets = int(math.ceil(float(total_file_count)/1000.0))

    for index in xrange(num_offsets):
        offset = 1000*index
        request = "%s?access_token=%s&limit=1000&offset=%s" % (href, access_token, offset)
        json_obj = restrequest(request)

        for f_json_obj in json_obj['Items']:
            target_path = os.path.abspath(os.path.join(output_dir, f_json_obj["Path"]))

            wget_url = "%s?access_token=%s" % (f_json_obj["HrefContent"], access_token)

            print "        Downloading %s to %s....(%i bytes)\n" % (f_json_obj["Path"], target_path, f_json_obj["Size"])
            downloaded_files.append(target_path)
            p = subprocess.Popen(["wget","-q","-O", target_path, wget_url])
            p.wait()

            print "\n"
            if os.path.getsize(target_path) != f_json_obj["Size"]:
                raise Exception("Downloading error. %s incorrect size (%i vs %i)" % (target_path, os.path.getsize(target_path), f_json_obj["Size"]))

    return downloaded_files

def main(project_id, access_token, merge, output_folder):

    # check no temp folder exists if lane merging is selected
    if merge:
        output_folder_temp = os.path.join(output_folder, 'temp')
        if os.path.exists(output_folder_temp):
            print 'Temp FASTQ input directory (%s) already exists! Exiting...\n' % output_folder_temp
            raise SystemExit
        else:
            os.mkdir(output_folder_temp)
            output_dir = output_folder_temp

    else:
        output_dir = output_folder

    request = base_url + ("v2/projects/%s/datasets?access_token=%s" % (project_id, access_token))
    json_obj = restrequest(request)

    print "Querying Project Samples: "

    total_count = json_obj["Paging"]["TotalCount"]

    if total_count == 0:
        print "    No samples found."
        sys.exit(0)

    print "    %s samples found." % total_count

    chunk_size = 50
    print "    Retrieving details in chunks of %s." % chunk_size

    sample_hrefs = {}
    for chunk in range(0, ((total_count - 1) / chunk_size) + 1):
        offset = chunk * chunk_size
        request = base_url + ("v2/projects/%s/datasets?access_token=%s&limit=%s&offset=%s" % (
            project_id,
            access_token,
            chunk_size,
            offset
        ))
        json_obj = restrequest(request)
        for sample_json in json_obj["Items"]:
            print "    Sample found - %s" % sample_json["Name"]
            sample_hrefs[sample_json["Name"]] = sample_json["HrefFiles"]

    print "Querying Individual Samples: "

    all_fastq = []
    sample_names = []

    for sample_href in sample_hrefs:

        print "    Getting Files For Sample: %s" % sample_href

        path_to_fastq = download_files(sample_hrefs[sample_href], access_token, output_dir)

        all_fastq.append(path_to_fastq)
        sample_names.append(sample_href[:-5])

    # option to combine samples from multiple lanes into a single file (lane 000)
    if merge:

        # for each sample, merge lanes

        samples = list(set(sample_names))

        for i in range(len(samples)):

            print 'Merging sample %s' % samples[i]

            sample_fastq = [all_fastq[j] for j in range(len(sample_names)) if sample_names[j] == samples[i]]
            sample_fastq = [item for sublist in sample_fastq for item in sublist]

            # merge R1
            r1_files = re.compile('.*_R1_.*')
            to_merge_r1 = filter(r1_files.match, sample_fastq)
            to_merge_r1.sort()
            r1_new = re.sub(r'_L\d\d\d_R1_', r'_L000_R1_', to_merge_r1[0])
            r1_new = os.path.join(output_folder, os.path.basename(r1_new))

            print "     Merging files %s into %s..." % (', '.join(to_merge_r1), r1_new)
            subprocess.call('cat %s > %s' % (' '.join(to_merge_r1), r1_new), shell=True)

            # merge R2 - if R2 files exist
            try:
                r2_files = re.compile('.*_R2_.*')
                to_merge_r2 = filter(r2_files.match, sample_fastq)
                to_merge_r2.sort()
                r2_new = re.sub(r'_L\d\d\d_R2_', r'_L000_R2_', to_merge_r2[0])
                r2_new = os.path.join(output_folder, os.path.basename(r2_new))

                print "     Merging files %s into %s..." % (', '.join(to_merge_r2), r2_new)
                subprocess.call('cat %s > %s' % (' '.join(to_merge_r2), r2_new), shell=True)

            except IndexError:
                pass

        # delete temp folder for unmerged samples
        shutil.rmtree(output_folder_temp)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''
    basespace fastq downloader
    ben demaree 2019
    
    Takes as input a BaseSpace Project (not Run!) ID and access token string and downloads the associated fastq.gz files
    to a chosen directory. Can also merge data across lanes, if selected. Output folder must already exist.
    
    A BaseSpace access token can be created by making an "app" at developer.basespace.illumina.com.

    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('project_id', type=str, help='BaseSpace PROJECT ID (numerical)')

    parser.add_argument('access_token', type=str, help='BaseSpace access token')

    parser.add_argument('output_folder', type=str, help='Output folder for FASTQ files')

    parser.add_argument('--merge', action='store_true', default=False, help='Option to merge lanes (default: no merging)')

    args = parser.parse_args()  # parse arguments

    project_id = args.project_id
    access_token = args.access_token
    output_folder = os.path.abspath(os.path.expanduser(args.output_folder))
    merge = args.merge

    # check that output folder exists already
    if not os.path.isdir(output_folder):
        raise OSError(output_folder + ' is not a directory. Exiting...')

    # download files
    main(project_id, access_token, merge, output_folder)





