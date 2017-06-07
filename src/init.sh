#!/bin/sh

# install nextflow
cd bin
wget -qO- get.nextflow.io | bash
pip install argparse
pip install docker
pip2 install multiqc
pip2 install click
pip2 install jinja2
sudo pip install cutadapt

cd /tmp
git clone https://github.com/COMBINE-lab/salmon.git
cd /tmp/salmon
mkdir build
cmake -DFETCH_BOOST=TRUE -DCMAKE_INSTALL_PREFIX=/usr/
make
sudo make install
rm -Rf /tmp/salmon

cd /tmp
git clone https://github.com/pachterlab/kallisto.git
cd /tmp/kallisto
mkdir build
cd build
cmake ..
make
sudo make install
rm -Rf /tmp/kallisto

cd /tmp
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip/download
unzip download
cd bowtie2-2.3.2
sudo cp bowtie2* /usr/bin/
rm -Rf /tmp/bowtie2* /tmp/download


apt-get install build-essential python2.7-dev python2-numpy python2-matplotlib
mkdir /tmp/HTseq
cd /tmp/HTseq
wget  https://pypi.python.org/packages/46/f7/6105848893b1d280692eac4f4f3c08ed7f424cec636aeda66b50bbcf217e/HTSeq-0.7.2.tar.gz#md5=8ddaaf53e83547fbca3bba7b84c9dde8
tar xvf HTSeq*
cd HTSeq-0.7.2
python setup.py build
sudo python setup.py install
rm -Rf /tmp/HTseq

mkdir /tmp/rsem
cd /tmp/rsem
git clone https://github.com/deweylab/RSEM.git
cd RSEM
make -j 8
sudo make install
rm -Rf /tmp/rsem
