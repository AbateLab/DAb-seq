#!/bin/bash

# script for running the complete dabseq pipeline in a docker container
# Ben Demaree 2020

# username (must have read/write access to output directory)
USER=bdemaree

# base input/output directory
BASE_DIR=/drive3/testing

# cohort name and directory
COHORT_NAME=cohort1
COHORT_DIR=$BASE_DIR/$COHORT_NAME

# sample names
SAMPLE_NAME_1=timepoint1
SAMPLE_NAME_2=timepoint2

# cfg file location
CFG_FILE_DIR=/home/bdemaree/code/dab-seq/cfg
CFG_FILE_NAME=dabseq.docker.cfg

# misc
SLACK_TOKEN="$(cat /home/bdemaree/.slack_token)"

### BARCODING ###

# note - to run using local (not GitHub) code, add the following argument to the docker call:
#-v /home/bdemaree/code/dab-seq/:/dabseq/pipeline/DAb-seq/ \

docker run \
-it \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $BASE_DIR/$COHORT_NAME/$SAMPLE_NAME_1/:/input/$COHORT_NAME/$SAMPLE_NAME_1/ \
-v $CFG_FILE_DIR/$CFG_FILE_NAME:/dabseq/config/$CFG_FILE_NAME \
dab-seq $COHORT_NAME barcode /dabseq/config/$CFG_FILE_NAME \
--sample-name $SAMPLE_NAME_1 \
--chem V1 \
--slack-token $SLACK_TOKEN

docker run \
-it \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $BASE_DIR/$COHORT_NAME/$SAMPLE_NAME_2/:/input/$COHORT_NAME/$SAMPLE_NAME_2/ \
-v $CFG_FILE_DIR/$CFG_FILE_NAME:/dabseq/config/$CFG_FILE_NAME \
dab-seq $COHORT_NAME barcode /dabseq/config/$CFG_FILE_NAME \
--sample-name $SAMPLE_NAME_2 \
--chem V1 \
--slack-token $SLACK_TOKEN

### GENOTYPING ###

docker run \
-it \
-u $(id -u ${USER}):$(id -g ${USER}) \
-v $COHORT_DIR/:/input/$COHORT_NAME/ \
-v $CFG_FILE_DIR/$CFG_FILE_NAME:/dabseq/config/$CFG_FILE_NAME \
dab-seq $COHORT_NAME genotype /dabseq/config/$CFG_FILE_NAME \
--slack-token $SLACK_TOKEN

