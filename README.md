# DAb-seq
DAb-seq: combined single-cell DNA genotyping and protein quantification

### Dependencies

The following software packages should be installed and located on the user's PATH. The version numbers represent those used in the DAb-seq publication.

* GATK (4.1.3.0)
* bowtie2 (2.3.4.1)
* ITDseek (1.2)
* samtools (1.8)
* bedtools (2.27.1)
* bcftools (1.9)
* cutadapt (2.4)
* BBMap (38.57)
* snpEff (4.3t)

### Usage

The pipeline consists of two primary Python scripts, which must be run sequentially: `mb_pipeline.py` and `genotype_cohort.py`.

#### `mb_pipeline.py`

The first script, `mb_pipeline.py` processes raw FASTQ files according to settings in a configuration file (see the provided example, `dabseq.cfg`). The script demultiplexes DNA panel amplicons and antibody tags into single cells, aligns panel reads to the human genome, and generates a GVCF file for each cell. For specific usage instructions and a list of arguments, use the `-h` option.

#### `genotype_cohort.py`

The second script operates on the single-cell GVCF files produced by `mb_pipeline.py` for one or multiple experiments. The script imports GVCFs into a GenomicsDB database (GATK GenomicsDBImport) and calls genotypes (GATK GenotypeGVCFs) separately for each genomic interval in parallel.
