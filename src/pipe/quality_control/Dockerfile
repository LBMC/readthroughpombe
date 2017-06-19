FROM ubuntu:16.04
LABEL maintainer "laurent.modolo@ens-lyon.fr"

RUN apt-get update && apt-get install -y \
    wget=1.17.1-1ubuntu1.2 \
    pigz=2.3.1-2 \
    git=1:2.7.4-0ubuntu1.1 \
    python-pip=8.1.1-2ubuntu0.4 \
    python3-pip=8.1.1-2ubuntu0.4 \
    build-essential=12.1ubuntu2 \
    unzip=6.0-20ubuntu1 \
    trimmomatic=0.35+dfsg-1 \
    zlib1g-dev=1:1.2.8.dfsg-2ubuntu4.1 \
    libfreetype6=2.6.1-0.1ubuntu2.3 \
    libbz2-dev=1.0.6-8 \
    liblzma-dev=5.1.1alpha+20120614-2ubuntu2 \
   && rm -rf /var/lib/apt/lists/*
RUN pip3 install --system argparse==1.4.0
RUN pip install --system numpy==1.12.1 Cython==0.25.2 && \
  pip install --system multiqc==0.9 htseq==0.7.2 cutadapt==1.14

RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.3.21-Linux-x86_64.sh && \
  bash Miniconda3-4.3.21-Linux-x86_64.sh -b && \
  export PATH=/root/miniconda3/bin:$PATH && \
  conda install --yes -c bioconda cutadapt=1.13

#install UrQt
RUN git clone https://github.com/l-modolo/UrQt.git && \
  cd UrQt && \
  make && \
  cp UrQt /usr/local/bin/

#install fastqc
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
  unzip fastqc_v0.11.5.zip && \
  mv FastQC /opt/FastQC && \
  ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
  chmod +x /usr/local/bin/fastqc

COPY file_handle.py /usr/bin/local/file_handle.py
