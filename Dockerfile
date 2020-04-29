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
RUN apt-get install -y libcurl4-openssl-dev
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
RUN chmod 777 /dabseq/programs/bedtools/bedtools

# bbmap
WORKDIR /dabseq/programs
RUN wget -q --show-progress --progress=bar:force:noscroll https://sourceforge.net/projects/bbmap/files/BBMap_38.57.tar.gz
RUN gunzip BBMap_38.57.tar.gz
RUN tar -xf BBMap_38.57.tar
ENV PATH "$PATH:/dabseq/programs/bbmap"

# the following section is only needed for builds using a human reference

# idtseek
RUN git clone https://github.com/tommyau/itdseek
ENV PATH "$PATH:/dabseq/programs/itdseek"

# snpeff
RUN wget -q --show-progress --progress=bar:force:noscroll http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
RUN unzip snpEff_latest_core.zip
ENV PATH "$PATH:/dabseq/programs/snpEff"
ENV PATH "$PATH:/dabseq/programs/snpEff/scripts"
# install snpEff database for hg19
RUN snpEff download hg19

### download genome reference files
WORKDIR /dabseq/references

# get clinvar db
RUN wget -q --show-progress --progress=bar:force:noscroll https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2020/clinvar_20200329.vcf.gz
RUN gunzip clinvar_20200329.vcf.gz
RUN awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' clinvar_20200329.vcf > clinvar_20200329.chr.vcf
RUN bgzip -f -@ 16 clinvar_20200329.chr.vcf
RUN tabix clinvar_20200329.chr.vcf.gz

# get hg19 fasta and pre-built indices
# need to download pre-built files to avoid out-of-memory error on Docker Hub
RUN wget -q --show-progress --progress=bar:force:noscroll ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
RUN gunzip hg19.fa.gz
RUN mv hg19.fa hg19.fasta
#RUN gatk CreateSequenceDictionary -R hg19.fasta
#RUN samtools faidx hg19.fasta
RUN wget -q --show-progress --progress=bar:force:noscroll https://www.dropbox.com/s/6yzq1n5fwtp06hf/hg19.dict?dl=0 -O hg19.dict
RUN wget -q --show-progress --progress=bar:force:noscroll https://www.dropbox.com/s/jbfyt5uz7jnfrm1/hg19.fasta.fai?dl=0 -O hg19.fasta
RUN wget -q --show-progress --progress=bar:force:noscroll https://www.dropbox.com/s/asqiuaiyzvqnkj0/hg19_bt2.tar.gz?dl=0 -O hg19_bt2.tar.gz
RUN gunzip hg19_bt2.tar.gz
RUN tar -xf hg19_bt2.tar
RUN rm hg19_bt2.tar
RUN mv hg19_bt2/*.bt2 .

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
RUN pip install slackclient

# set input directories to fully open - so any user can write
RUN mkdir /input
RUN chmod -R 777 /input

### get dab-seq private repo files
# before building: export SSH_PRIVATE_KEY="$(cat /home/bdemaree/.ssh/dab_seq_key)"
ARG SSH_PRIVATE_KEY
RUN mkdir /root/.ssh/
RUN echo "${SSH_PRIVATE_KEY}" > /root/.ssh/id_rsa
RUN chmod 400 /root/.ssh/id_rsa
RUN touch /root/.ssh/known_hosts
RUN ssh-keyscan github.com >> /root/.ssh/known_hosts
WORKDIR /dabseq/pipeline
RUN git clone git@github.com:AbateLab/DAb-seq
ENTRYPOINT ["python", "/dabseq/pipeline/DAb-seq/dabseq_pipeline.py"]