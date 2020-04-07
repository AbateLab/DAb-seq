# DAb-seq
DAb-seq: combined single-cell DNA and Antibody sequencing.

## Introduction

DAb-seq is a multiomic platform combining targeted genotyping and immunophenotyping at the single-cell level. Through the use of DNA-antibody conjugates, phenotypic signal is encoded into next-generation sequencing data, providing a readout analagous to that of flow cytometry. The result is a dataset of linked proteogenomic information from thousands of single cells.

![DAb-seq Workflow](https://i.imgur.com/2Z2GTey.png)
<p align="center"><i>The experimental DAb-seq workflow, from sample collection to bioinformatics.</i><br></p>

The DAb-seq data analysis pipeline incoporates elements of targeted DNA genotyping and digital fragment counting, as in RNA-seq. The primary output of the pipeline is a genotyping matrix of calls by cell, and a counts matrix of antibody UMI counts by cell.

## Running the Pipeline

The pipeline is run through the main Python script, `dabseq_pipeline.py`, which must be run sequentially in `barcode` and `genotype` modes.

### Pipeline Input Files



### Software Dependencies

The following software packages should be installed and located on the user's PATH. The version numbers shown are those used in the DAb-seq publication.

* GATK (4.1.3.0)
* bowtie2 (2.3.4.1)
* ITDseek (1.2)
* samtools (1.8)
* bedtools (2.27.1)
* bcftools (1.9)
* cutadapt (2.4)
* BBMap (38.57)
* snpEff (4.3t)

To simplify installation and enhance data reproducibility, the pipeline can also be run in a Docker container. A Dockerfile is included in the repo and instructions for building and running the DAb-seq image are listed in the section [DAb-seq in Docker](###dab-seq-in-docker).

### `barcode` Mode

In `barcode` mode, the pipeline processes raw FASTQ files according to settings in a configuration file (see the provided example, `dabseq.cfg`). The script demultiplexes DNA panel amplicons and antibody tags into single cells, aligns panel reads to the human genome, and generates a GVCF file for each cell. For specific usage instructions and a list of arguments, use the `-h` option.

### `genotype` Mode

In `genotype` mode, the pipeline calls variants for samples that have been processed in `barcode` mode. The script imports GVCFs into a GenomicsDB database (GATK GenomicsDBImport) and calls genotypes (GATK GenotypeGVCFs) separately for each genomic interval in parallel.

## DAb-seq in Docker

Running the DAb-seq pipeline in a Docker container is recommended and ensures that the pipeline dependencies are installed and configured properly.
