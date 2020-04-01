# DAb-seq Dockerfile
# Ben Demaree 3.30.2020

# start with ubuntu base
FROM ubuntu:18.04

### install linux dependencies
RUN apt-get -y update
RUN apt-get install -y wget
RUN apt-get -y install python-pip
RUN apt-get -y install python3-pip
RUN apt-get install -y pigz
RUN apt-get install -y libncurses5-dev
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y libbz2-dev
RUN apt-get install -y liblzma-dev
RUN apt-get install -y unzip
RUN apt-get install -y git
RUN apt-get install -y default-jdk

### install bioinformatics programs

WORKDIR /dabseq/programs

# htslib
RUN wget -q --show-progress --progress=bar:force:noscroll https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
RUN tar -xjf htslib-1.9.tar.bz2
RUN mkdir htslib
WORKDIR /dabseq/programs/htslib-1.9
RUN ./configure --prefix=/dabseq/programs/htslib
RUN make
RUN make install
ENV PATH "$PATH:/dabseq/programs/htslib/bin"

# samtools
WORKDIR /dabseq/programs
RUN wget -q --show-progress --progress=bar:force:noscroll https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
RUN tar -xjf samtools-1.8.tar.bz2
RUN mkdir samtools
WORKDIR /dabseq/programs/samtools-1.8
RUN ./configure --prefix=/dabseq/programs/samtools
RUN make
RUN make install
ENV PATH "$PATH:/dabseq/programs/samtools/bin"

# bcftools
WORKDIR /dabseq/programs
RUN wget -q --show-progress --progress=bar:force:noscroll https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
RUN tar -xjf bcftools-1.9.tar.bz2
RUN mkdir bcftools
WORKDIR /dabseq/programs/bcftools-1.9
RUN ./configure --prefix=/dabseq/programs/bcftools
RUN make
RUN make install
ENV PATH "$PATH:/dabseq/programs/bcftools/bin"

# gatk
WORKDIR /dabseq/programs
RUN wget -q --show-progress --progress=bar:force:noscroll https://github.com/broadinstitute/gatk/releases/download/4.1.3.0/gatk-4.1.3.0.zip
RUN unzip gatk-4.1.3.0.zip
ENV PATH "$PATH:/dabseq/programs/gatk-4.1.3.0"

# bowtie2
RUN wget -q --show-progress --progress=bar:force:noscroll https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip
RUN unzip bowtie2-2.3.4.1-linux-x86_64.zip
ENV PATH "$PATH:/dabseq/programs/bowtie2-2.3.4.1-linux-x86_64"

# bedtools
RUN mkdir bedtools
WORKDIR /dabseq/programs/bedtools
RUN wget -q --show-progress --progress=bar:force:noscroll https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools
ENV PATH "$PATH:/dabseq/programs/bedtools"

# bbmap
WORKDIR /dabseq/programs
RUN wget -q --show-progress --progress=bar:force:noscroll https://sourceforge.net/projects/bbmap/files/BBMap_38.57.tar.gz
RUN gunzip BBMap_38.57.tar.gz
RUN tar -xf BBMap_38.57.tar
ENV PATH "$PATH:/dabseq/programs/bbmap"

# idtseek
RUN git clone https://github.com/tommyau/itdseek
ENV PATH "$PATH:/dabseq/programs/itdseek"

# snpeff
RUN wget -q --show-progress --progress=bar:force:noscroll http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
RUN unzip snpEff_latest_core.zip
ENV PATH "$PATH:/dabseq/programs/snpEff"
ENV PATH "$PATH:/dabseq/programs/snpEff/scripts"

### download genome reference files

WORKDIR /dabseq/references

# get clinvar db
RUN wget -q --show-progress --progress=bar:force:noscroll https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20200329.vcf.gz
RUN gunzip clinvar_20200329.vcf.gz
RUN awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' clinvar_20200329.vcf > clinvar_20200329.chr.vcf
RUN bgzip -f -@ 16 clinvar_20200329.chr.vcf
RUN tabix clinvar_20200329.chr.vcf.gz

# get hg19 fasta and build indices
RUN wget -q --show-progress --progress=bar:force:noscroll ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
RUN gunzip hg19.fa.gz
RUN mv hg19.fa hg19.fasta
RUN gatk CreateSequenceDictionary -R hg19.fasta
RUN samtools faidx hg19.fasta
RUN bowtie2-build hg19.fasta hg19 --threads 24

### install missing python modules
RUN pip3 install cutadapt
RUN pip install numpy
RUN pip install pandas
RUN pip install regex
RUN pip install scipy
RUN pip install scikit-allel
RUN pip install h5py
RUN pip install future
RUN pip install matplotlib

### get dab-seq private repo files

# run before building: export SSH_PRIVATE_KEY="$(cat /home/bdemaree/.ssh/dab_seq_key)"

ARG SSH_PRIVATE_KEY
RUN mkdir /root/.ssh/
RUN echo "${SSH_PRIVATE_KEY}" > /root/.ssh/id_rsa
RUN chmod 400 /root/.ssh/id_rsa

# add domain to known hosts
RUN touch /root/.ssh/known_hosts
RUN ssh-keyscan github.com >> /root/.ssh/known_hosts
WORKDIR /dabseq/pipeline
RUN git clone git@github.com:AbateLab/DAb-seq
#ENV PATH "$PATH:/dabseq/pipeline/DAb-seq"
#RUN echo $PATH
#ENV PYTHONPATH "$PYTHONPATH:/dabseq/pipeline/DAb-seq"
RUN python --version
RUN python /dabseq/pipeline/DAb-seq/mb_pipeline.py -h

