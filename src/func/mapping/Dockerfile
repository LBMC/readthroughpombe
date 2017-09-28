FROM ubuntu:16.04
MAINTAINER Laurent Modolo

ENV SALMON_VERSION 0.8.2
ENV KALLISTO_VERSION 0.43.1
ENV BOWTIE2_VERSION 2.3.2
ENV SAMTOOLS_VERSION 1.5

ENV PACKAGES curl=7.47.0-1ubuntu2.2 \
    git=1:2.7.4-0ubuntu1.1 \
    unzip=6.0-20ubuntu1 \
    pigz=2.3.1-2 \
    ca-certificates=20160104ubuntu1 \
    python-pip=8.1.1-2ubuntu0.4 \
    python3-pip=8.1.1-2ubuntu0.4 \
    g++=4:5.3.1-1ubuntu1 \
    gcc=4:5.3.1-1ubuntu1 \
    autoconf=2.69-9 \
    pkg-config=0.29.1-0ubuntu1 \
    make=4.1-6 \
    cmake=3.5.1-1ubuntu3 \
    libboost-all-dev=1.58.0.1ubuntu1 \
    liblzma-dev=5.1.1alpha+20120614-2ubuntu2 \
    libbz2-dev=1.0.6-8 \
    zlib1g-dev=1:1.2.8.dfsg-2ubuntu4.1 \
    libjemalloc-dev=3.6.0-9ubuntu1 \
    libdivsufsort-dev=2.0.1-2 \
    libtbb2=4.4~20151115-0ubuntu3 \
    libncurses5=6.0+20160213-1ubuntu1 \
    libncurses5-dev=6.0+20160213-1ubuntu1 \
    bedtools=2.25.0-1

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN pip3 install --system argparse==1.4.0
RUN pip install --system numpy==1.12.1 setuptools==36.0.1 && \
    pip install --system multiqc==0.9 htseq==0.8.0

RUN curl -k -L  https://github.com/COMBINE-lab/salmon/archive/v${SALMON_VERSION}.tar.gz -o salmon-v${SALMON_VERSION}.tar.gz && \
    tar xzf salmon-v${SALMON_VERSION}.tar.gz && \
    cd salmon-${SALMON_VERSION} && \
    mkdir build && \
    cd build && \
    sh ../scripts/fetchRapMap.sh && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install

RUN  curl -k -L  https://github.com/pachterlab/kallisto/releases/download/v${KALLISTO_VERSION}/kallisto_linux-v${KALLISTO_VERSION}.tar.gz -o kallisto_linux-v${KALLISTO_VERSION}.tar.gz && \
  tar xzf kallisto_linux-v${KALLISTO_VERSION}.tar.gz && \
  cp kallisto_linux-v${KALLISTO_VERSION}/kallisto /usr/local/bin/ && \
  cd .. && rm -Rf kallisto_linux-v${KALLISTO_VERSION}*

RUN curl -k -L  https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip/download -o bowtie2-${BOWTIE2_VERSION}.zip && \
  unzip bowtie2-${BOWTIE2_VERSION}.zip && \
  cd bowtie2-${BOWTIE2_VERSION} && \
  cp bowtie2* /usr/local/bin/ && \
  cd .. && rm -Rf bowtie2*

RUN  git clone https://github.com/deweylab/RSEM.git && \
  cd RSEM && \
  git checkout tags/v1.3.0 && \
  make && \
  make install && \
  cd .. && rm -Rf RSEM

RUN curl -k -L https://sourceforge.net/projects/samtools/files/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2/download -o samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
  tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
  cd samtools-${SAMTOOLS_VERSION} && \
  ./configure && \
  make && \
  make install && \
  cd .. && rm -Rf samtools-${SAMTOOLS_VERSION}*

COPY file_handle.py /usr/bin/local/file_handle.py
