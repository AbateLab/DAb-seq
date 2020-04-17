# DAb-seq
<i>DAb-seq: combined single-cell DNA and Antibody sequencing.</i>

View the current version of the [DAb-seq manuscript](https://www.biorxiv.org/content/10.1101/2020.02.26.967133v1) on BioRxiv. 

## Introduction

DAb-seq is a multiomic platform combining targeted genotyping and immunophenotyping at the single-cell level. Through the use of DNA-antibody conjugates, phenotypic signal is encoded into next-generation sequencing data, providing a readout analagous to that of flow cytometry. The result is a dataset of linked proteogenomic information from thousands of single cells.

![DAb-seq Workflow](https://i.imgur.com/2Z2GTey.png)
<p align="center"><i>The experimental DAb-seq workflow, from sample collection to bioinformatics.</i><br></p>

The DAb-seq data analysis pipeline incoporates elements of targeted DNA genotyping and digital fragment counting, as in RNA-seq. The primary output of the pipeline is a genotyping matrix of calls by cell, and a counts matrix of antibody UMI counts by cell.

## Input File Requirements

The pipeline requires as input, at a minimum, a configuration file and paired-end, compressed FASTQ files (ending in ".fastq.gz").

The configuration file is subdivided into sections for each analysis module (see the provided example, `dabseq.cfg`). Users should be sure to indicate the correct panel files (.bed and .amplicons). Mission Bio's standard AML and tumor hotspot (THP) panels are built-in and need not be provided separately.

The folder structure of the input FASTQ files needs to be set up in a specific way to allow the program to find everything. For a given cohort (CohortA) of samples to joint-genotype (Timepoint1 and Timepoint2), the folder should look like:
```
Timepoint1 DNA Panel:   .../CohortA/Timepoint1/fastq/panel/<filename.fastq.gz>
Timepoint1 Abs:         .../CohortA/Timepoint1/fastq/abs/<filename.fastq.gz>

Timepoint2 DNA Panel:   .../CohortA/Timepoint2/fastq/panel/<filename.fastq.gz>
Timepoint2 Abs:         .../CohortA/Timepoint2/fastq/abs/<filename.fastq.gz>
```
For each `<filename.fastq.gz>`, both R1 and R2 files need to be present. When sequencing multiple tubes of Tapestri output (e.g. grouping tubes 1-4 and 5-8 into two libraries), FASTQ files from multiple tubes should be placed in the same folder. Users should verify that the panel and antibody filenames remain in the same order when sorted lexicographically (a simple filenaming scheme like panel-A/abs-A, panel-B/abs-B, etc.. works well).

When running the pipeline in `dna-only` or `ab-only` modes, the user is not required to create the folders for the missing file types.

## Software Dependencies

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

To simplify installation and enhance data reproducibility, the pipeline can also be run in a Docker container. Instructions for building and running the DAb-seq image are listed in the section [DAb-seq in Docker](##dab-seq-in-docker).

## Running the Pipeline

The pipeline is run through the main Python script, `dabseq_pipeline.py`, which must be run sequentially in `barcode` and `genotype` modes.

### `barcode` Mode

In `barcode` mode, the pipeline processes raw FASTQ files according to settings in a configuration file . The script demultiplexes DNA panel amplicons and antibody tags into single cells, aligns panel reads to the human genome, and generates a GVCF file for each cell. All samples belonging to the same cohort must be individually barcoded before joint genotyping.

### `genotype` Mode

In `genotype` mode, the pipeline calls variants for a single cohort, containing samples that have been individually processed in `barcode` mode. The script imports GVCFs into a GenomicsDB database (GATK GenomicsDBImport) and calls genotypes (GATK GenotypeGVCFs) separately for each genomic interval in parallel.

## DAb-seq in Docker

Running the DAb-seq pipeline in a Docker container is recommended and ensures that the pipeline dependencies are installed and configured properly.

1. Build the DAb-seq image using Docker (it will take some time to build the Bowtie2 index):
```
docker build <path_to_dabseq_repo> -t dab-seq:latest
```

2. Organize the input files on the host machine in the same file structure as described in [Input File Requirements](##input-file-requirements).

3. Edit the included bash script `run_dabseq_docker.sh` with the appropriate cohort and sample information. Samples to barcode can be added as needed. The user running the pipeline should also be specified. This user should also have read/write access to the input FASTQ files.

4. Run `run_dabseq_docker.sh`. You may need to change the file permissions to allow execution before running (e.g. `chmod +x run_dabseq_docker.sh`).

That's it! The pipeline will run in the Docker container and produce output files at the mounted locations on the host machine.

## Memory and CPU Considerations

DNA genotyping is CPU and memory-intensive for a single bulk sample, let alone thousands of single cells. The DAb-seq pipeline is implemented with tunable parallelization to scale to the resources available on different systems. In its default configuration, the pipeline requires at least 64 gb of memory and 16 physical threads to run. Changing the amount of cells aligned simultaneously or the number of genomic intervals genotyped in parallel can reduce this hardware requirement significantly, with a corresponding increase in total processing time.

## Output Files


![DAb-seq Output](https://i.imgur.com/rUMK27M.png)
<p align="center"><i>The DAb-seq output can be represented as linked genotyping and antibody UMI count tables.</i><br></p>

Output genotyping data is saved in a new directory labeled `GENOTYPING` in the root of each cohort. A compresed HDF5 file contains matrices of discrete genotyping calls and antibody UMI counts for all cells in the cohort. Genotyping calls for each variant assume a genome ploidy of 2. Therefore, there are four possible genotypes per matrix entry:

* Wildtype (0)
* Heterozygous Mutant (1)
* Homozygous Mutant (2)
* No Call (3)

Further details on reading and manipulating data from this file can be found in the included Python notebooks.

## Non-Human Organism Support

The DAb-seq pipeline includes support for non-human organisms using th flag `--non-human`. When combined with additional settings such as `--ploidy 1`, the pipeline can be used to perform single-cell genotyping on haploid bacteria and yeast cells.
