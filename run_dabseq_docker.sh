#!/bin/bash

# script for running the complete dabseq pipeline in a docker container
# Ben Demaree 2021

# username (must have read/write access to output directory)
USER=bdemaree

# base input/output directory
BASE_DIR=[]

# cohort name and directory
COHORT_NAME=[]
COHORT_DIR=$BASE_DIR/$COHORT_NAME

# sample names
SAMPLE_NAME_1=[]
# add more sample names here

# directory containing files for this run
# will be mounted into container at /dabseq/config
RUN_DIR=[]

# files in run directory needed for processing
CFG_FILE=dabseq.hg19.cfg
AB_CSV_FILE=[]
DNA_PANEL_FILE=[]

# container to run
CONTAINER=bendemaree/dab-seq:human

# misc
SLACK_TOKEN="$(cat /home/bdemaree/.slack_token)"

### BARCODING ###

docker run \
-it \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $BASE_DIR/$COHORT_NAME/$SAMPLE_NAME_1/:/input/$COHORT_NAME/$SAMPLE_NAME_1/ \
-v $RUN_DIR:/dabseq/config/ \
$CONTAINER $COHORT_NAME barcode /dabseq/config/$CFG_FILE \
--sample-name $SAMPLE_NAME_1 \
--dna-panel /dabseq/config/$DNA_PANEL_FILE \
--ab-csv /dabseq/config/$AB_CSV_FILE \
--chem V2 \
# add barcoding code blocks here for additional samples

### GENOTYPING ###

docker run \
-it \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $COHORT_DIR/:/input/$COHORT_NAME/ \
-v $RUN_DIR:/dabseq/config/ \
$CONTAINER $COHORT_NAME genotype /dabseq/config/$CFG_FILE \
--dna-panel /dabseq/config/$DNA_PANEL_FILE