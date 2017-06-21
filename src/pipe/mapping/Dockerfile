FROM ubuntu:16.04
MAINTAINER Laurent Modolo

RUN apt-get update && apt-get install -y \
    wget=1.17.1-1ubuntu1.2 \
    git=1:2.7.4-0ubuntu1.1 \
    python-pip=8.1.1-2ubuntu0.4 \
    python3-pip=8.1.1-2ubuntu0.4 \
    build-essential=12.1ubuntu2 \

   && rm -rf /var/lib/apt/lists/*
RUN pip3 install argparse==1.4.0
RUN pip install numpy==1.12.1 && \
  pip install multiqc==0.9 htseq==0.7.2

RUN git clone https://github.com/COMBINE-lab/salmon.git && \
  cd salmon && \
  git checkout tags/v0.8.2 -b master && \
  mkdir build && \
  cmake -DFETCH_BOOST=TRUE -DCMAKE_INSTALL_PREFIX=/usr/  && \
  make  && \
  make install && \
  cd .. && rm -Rf salmon

RUN  git clone https://github.com/pachterlab/kallisto.git && \
  cd kallisto && \
  git checkout tags/v0.43.1 -b master && \
  mkdir build && \
  cmake && \
  make && \
  make install && \
  cd .. && rm -Rf kallisto

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip/download && \
  unzip download && \
  cd bowtie2-2.3.2 && \
  cp bowtie2* /usr/bin/ && \
  cd .. && rm -Rf bowtie2* download

RUN  git clone https://github.com/deweylab/RSEM.git
  cd RSEM && \
  git checkout tags/v1.3.0 -b master && \
  make && \
  sudo make install && \
  cd .. && rm -Rf RSEM
