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
cd build
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
