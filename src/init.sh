#!/bin/sh

############################ we rename the data ################################
mv /home/laurent/data/readthroughpombe/data/fastq /home/laurent/data/readthroughpombe/data/fastq_original
ln -s /home/laurent/data/readthroughpombe/data/fastq_original data/fastq_original
mkdir /home/laurent/data/readthroughpombe/data/fastq
cp data/fastq_original/PLBD1.fastq.gz	data/fastq/2017_04_14_PLBD1_wt_R1.fastq.gz
cp data/fastq_original/PLBD2.fastq.gz	data/fastq/2017_04_14_PLBD2_cut14_208_R1.fastq.gz
cp data/fastq_original/PLBD3.fastq.gz	data/fastq/2017_04_14_PLBD3_cut14_208 cdc15_118_R1.fastq.gz
cp data/fastq_original/PLBD4.fastq.gz	data/fastq/2017_04_14_PLBD4_cdc15_118_R1.fastq.gz
cp data/fastq_original/PLBD5.fastq.gz	data/fastq/2017_04_14_PLBD5_rrp6D_R1.fastq.gz
cp data/fastq_original/PLBD6.fastq.gz	data/fastq/2017_04_14_PLBD6_wt_R2.fastq.gz
cp data/fastq_original/PLBD7.fastq.gz	data/fastq/2017_04_14_PLBD7_cut14_208_R2.fastq.gz
cp data/fastq_original/PLBD8.fastq.gz	data/fastq/2017_04_14_PLBD8_cut14_208_cdc15_118_R2.fastq.gz
cp data/fastq_original/PLBD9.fastq.gz	data/fastq/2017_04_14_PLBD9_cdc15_118_R2.fastq.gz
cp data/fastq_original/PLBD10.fastq.gz	data/fastq/2017_04_14_PLBD10_rrp6D_R2.fastq.gz
cp data/fastq_original/PLBD11.fastq.gz	data/fastq/2017_04_14_PLBD11_wt_R3.fastq.gz
cp data/fastq_original/PLBD12.fastq.gz	data/fastq/2017_04_14_PLBD12_cut14_208_R3.fastq.gz
cp data/fastq_original/PLBD13.fastq.gz	data/fastq/2017_04_14_PLBD13_cut14_208_cdc15_118_R3.fastq.gz
cp data/fastq_original/PLBD14.fastq.gz	data/fastq/2017_04_14_PLBD14_cdc15_118_R3.fastq.gz
cp data/fastq_original/PLBD15.fastq.gz	data/fastq/2017_04_14_PLBD15_rrp6D_R3.fastq.gz

################################### all ########################################

# install nextflow
mkdir bin
cd bin
wget -qO- get.nextflow.io | bash
cd ..

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


################################################################################
# readthrough detection
docker build src/func/docker_readthrough -t 'readthrough:0.0.1'
