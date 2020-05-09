#!/bin/bash

# script for running the complete dabseq pipeline in a singularity container on the wynton cluster
# Ben Demaree 2020

# base input/output directory
BASE_DIR=/wynton/home/abate/bdemaree/dabseq/cohorts

# cohort name and directory
COHORT_NAME=cohort1
COHORT_DIR=$BASE_DIR/$COHORT_NAME

# sample names
SAMPLE_NAME_1=timepoint1
SAMPLE_NAME_2=timepoint2

# cfg file location
CFG_FILE_DIR=/wynton/home/abate/bdemaree/code/DAb-seq/cfg
CFG_FILE_NAME=dabseq.docker.cfg

# singularity container to run
CONTAINER=/wynton/home/abate/bdemaree/singularity/dabseq.human.sif

# misc
SLACK_TOKEN="$(cat /wynton/home/abate/bdemaree/.slack_token)"

### BARCODING ###

# note - to run using local (not GitHub) code, add the following argument to the docker call:
#-v /wynton/home/abate/bdemaree/code/DAb-seq/:/dabseq/pipeline/DAb-seq/ \

singularity run \
--bind $BASE_DIR/$COHORT_NAME/$SAMPLE_NAME_1/:/input/$COHORT_NAME/$SAMPLE_NAME_1/,\
$CFG_FILE_DIR/$CFG_FILE_NAME:/dabseq/config/$CFG_FILE_NAME,\
/wynton/home/abate/bdemaree/code/DAb-seq/:/dabseq/pipeline/DAb-seq/ \
$CONTAINER $COHORT_NAME barcode /dabseq/config/$CFG_FILE_NAME \
--sample-name $SAMPLE_NAME_1 \
--chem V1 \
--slack-token $SLACK_TOKEN

singularity run \
--bind $BASE_DIR/$COHORT_NAME/$SAMPLE_NAME_2/:/input/$COHORT_NAME/$SAMPLE_NAME_2/,\
$CFG_FILE_DIR/$CFG_FILE_NAME:/dabseq/config/$CFG_FILE_NAME,\
/wynton/home/abate/bdemaree/code/DAb-seq/:/dabseq/pipeline/DAb-seq/ \
$CONTAINER $COHORT_NAME barcode /dabseq/config/$CFG_FILE_NAME \
--sample-name $SAMPLE_NAME_2 \
--chem V1 \
--slack-token $SLACK_TOKEN

### GENOTYPING ###

singularity run \
--bind $COHORT_DIR/:/input/$COHORT_NAME/,\
$CFG_FILE_DIR/$CFG_FILE_NAME:/dabseq/config/$CFG_FILE_NAME,\
/wynton/home/abate/bdemaree/code/DAb-seq/:/dabseq/pipeline/DAb-seq/ \
$CONTAINER $COHORT_NAME genotype /dabseq/config/$CFG_FILE_NAME \
--slack-token $SLACK_TOKEN

