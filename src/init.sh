#!/bin/sh

################################### all ########################################
# install nextflow
cd bin
wget -qO- get.nextflow.io | bash

################################## Docker ######################################
# build docker image

docker build src/func/docker -t 'pipeline:0.0.1'

################################### PSMN #######################################
# echo $SHELL must return /bin/bash for this to work
ln -s /Xnfs/site/lbmcdb/common/modules/modulefiles ~/privatemodules

# then add the two following line to the file ~/.bashrc
source /usr/share/modules/init/bash
module load use.own

################################## local #######################################

# install python3 modules
pip install argparse
pip install docker
pip install numpy
sudo pip install cutadapt

# install python2 modules
pip2 install multiqc
pip2 install click
pip2 install jinja2
sudo pip2 install htseq

# install pigz
sudo apt-get install pigz

# install salmon
cd /tmp
git clone https://github.com/COMBINE-lab/salmon.git
cd /tmp/salmon
mkdir build
cmake -DFETCH_BOOST=TRUE -DCMAKE_INSTALL_PREFIX=/usr/local/
make
sudo make install
rm -Rf /tmp/salmon

# install kallisto
cd /tmp
git clone https://github.com/pachterlab/kallisto.git
cd /tmp/kallisto
mkdir build
cd build
cmake ..
make
sudo make install
rm -Rf /tmp/kallisto

# install bowtie2
cd /tmp
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip/download
unzip download
cd bowtie2-2.3.2
sudo cp bowtie2* /usr/local/bin/
rm -Rf /tmp/bowtie2* /tmp/download

# install RSEM
mkdir /tmp/rsem
cd /tmp/rsem
git clone https://github.com/deweylab/RSEM.git
cd RSEM
make -j 8
sudo make install
rm -Rf /tmp/rsem

#install fastqc
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
  unzip fastqc_v0.11.5.zip && \
  mv FastQC /opt/FastQC && \
  ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
  chmod +x /usr/local/bin/fastqc
