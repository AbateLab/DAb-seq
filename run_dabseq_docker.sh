#!/bin/bash

# script for running the complete dabseq pipeline in a docker container
# Ben Demaree 2020

# username (must have read/write access to output directory)
USER=bdemaree

# base input/output directory
BASE_DIR=/drive3/dabseq/cohorts

# cohort name and directory
COHORT_NAME=hashing_pilot_1
COHORT_DIR=$BASE_DIR/$COHORT_NAME

# sample names
SAMPLE_NAME_1=all_tubes
#SAMPLE_NAME_2=timepoint2

# directory containing files for this run
# will be mounted into container at /dabseq/config
RUN_DIR=/home/bdemaree/code/dab-seq/runs/dabseq

# files in run directory needed for processing
CFG_FILE=dabseq.hg19.cfg
AB_CSV_FILE=ab_barcodes_hashpilot1.csv
DNA_PANEL_FILE=AML

# container to run
CONTAINER=bendemaree/dab-seq:human

# misc
SLACK_TOKEN="$(cat /home/bdemaree/.slack_token)"

### BARCODING ###

# note - to run using local (not GitHub) code, add the following argument to the docker calls:
#-v /home/bdemaree/code/dab-seq/:/dabseq/pipeline/DAb-seq/ \

docker run \
-it \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v /home/bdemaree/code/dab-seq/:/dabseq/pipeline/DAb-seq/ \
-v $BASE_DIR/$COHORT_NAME/$SAMPLE_NAME_1/:/input/$COHORT_NAME/$SAMPLE_NAME_1/ \
-v $RUN_DIR:/dabseq/config/ \
$CONTAINER $COHORT_NAME barcode /dabseq/config/$CFG_FILE \
--sample-name $SAMPLE_NAME_1 \
--dna-panel /dabseq/config/$DNA_PANEL_FILE \
--ab-csv /dabseq/config/$AB_CSV_FILE \
--chem V2 \
--slack-token $SLACK_TOKEN

#docker run \
#-it \
#-u $(id -u ${USER}):$(id -g ${USER}) \
#-v /home/bdemaree/code/dab-seq/:/dabseq/pipeline/DAb-seq/ \
#-v $BASE_DIR/$COHORT_NAME/$SAMPLE_NAME_2/:/input/$COHORT_NAME/$SAMPLE_NAME_2/ \
#-v $CFG_FILE_DIR/$CFG_FILE:/dabseq/config/$CFG_FILE \
#-v $AB_CSV_DIR/$AB_CSV_FILE:/dabseq/config/$AB_CSV_FILE \
#$CONTAINER $COHORT_NAME barcode /dabseq/config/$CFG_FILE \
#--sample-name $SAMPLE_NAME_2 \
#--dna-panel /dabseq/config/$DNA_PANEL_FILE \
#--ab-csv /dabseq/config/$AB_CSV_FILE \
#--chem V2 \
#--slack-token $SLACK_TOKEN

### GENOTYPING ###

docker run \
-it \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v /home/bdemaree/code/dab-seq/:/dabseq/pipeline/DAb-seq/ \
-v $COHORT_DIR/:/input/$COHORT_NAME/ \
-v $RUN_DIR:/dabseq/config/ \
$CONTAINER $COHORT_NAME genotype /dabseq/config/$CFG_FILE \
--dna-panel /dabseq/config/$DNA_PANEL_FILE \
--slack-token $SLACK_TOKEN