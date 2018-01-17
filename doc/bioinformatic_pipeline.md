# Reads Quality control

The first part of the pipeline is the quality control of the data.
For this we need 4 steps:

- check the quality of the fastq files
- remove (if any) the adaptor of the reads
- trim the reads to remove nucleotid of poor quality
- check the quality of the results

# Reads Mapping

The second part of the pipeline is to map the RNASeq data on the reference
genomes. For this we need 2 steps:

- indexing the reference genomes
- mapping the reads on the references

all the steps of quality control and mapping are computed with the
`src/pipeline.nf` script. To use this pipeline, one must first run the
appropriate sections of the `src/init.sh` script, then run the commands from the
`/src/pipe/quantif.sh` script.

# Readthrough Detection

The detection of readthrough is done with the `src/readthrough.nf` script. To
use this pipeline, one must first un the appropriate sections of the
`src/init.sh` script, then run the commands from the `src/pipe/readthrough.sh`
script.
