# project initialisation

Start by executing the `src/init.sh` script to install nextflow and the
necessary dependencies. This script also perform the following renaming of the
original fastq files:

|old name        | new names                                           |
|----------------|-----------------------------------------------------|
|PLBD1.fastq.gz	| 2017_04_14_PLBD1_wt_R1.fastq.gz                      |
|PLBD2.fastq.gz	| 2017_04_14_PLBD2_cut14_208_R1.fastq.gz               |
|PLBD3.fastq.gz	| 2017_04_14_PLBD3_cut14_208 cdc15_118_R1.fastq.gz     |
|PLBD4.fastq.gz	| 2017_04_14_PLBD4_cdc15_118_R1.fastq.gz               |
|PLBD5.fastq.gz	| 2017_04_14_PLBD5_rrp6D_R1.fastq.gz                   |
|PLBD6.fastq.gz	| 2017_04_14_PLBD6_wt_R2.fastq.gz                      |
|PLBD7.fastq.gz	| 2017_04_14_PLBD7_cut14_208_R2.fastq.gz               |
|PLBD8.fastq.gz	| 2017_04_14_PLBD8_cut14_208_cdc15_118_R2.fastq.gz     |
|PLBD9.fastq.gz	| 2017_04_14_PLBD9_cdc15_118_R2.fastq.gz               |
|PLBD10.fastq.gz	| 2017_04_14_PLBD10_rrp6D_R2.fastq.gz                |
|PLBD11.fastq.gz	| 2017_04_14_PLBD11_wt_R3.fastq.gz                   |
|PLBD12.fastq.gz	| 2017_04_14_PLBD12_cut14_208_R3.fastq.gz            |
|PLBD13.fastq.gz	| 2017_04_14_PLBD13_cut14_208_cdc15_118_R3.fastq.gz  |
|PLBD14.fastq.gz	| 2017_04_14_PLBD14_cdc15_118_R3.fastq.gz            |
|PLBD15.fastq.gz	| 2017_04_14_PLBD15_rrp6D_R3.fastq.gz                |

# Quality control analysis

To run the quality control analysis you can launch the `src/pipe/QC.sh` script

# Mapping and quantification

To run the mapping and quantification analysis you can launch the
`src/pipe/quantif.sh` script

# RNASeq pipeline

The RNASeq pipeline is based on [nextflow](https://www.nextflow.io/). This pipeline aim is to automatize the first bioinformatic steps of an NGS analysis with a focus on RNASeq.
With this pipeline you will be able to perform:

- Quality control
- Indexing
- Mapping
- Quantification

Using different tools.

If you want another tool included or if you find a bug. Please open a new issue for it [at this address](http://gitlab.biologie.ens-lyon.fr/pipelines/RNASeq/issues).


## project initialisation

Start by executing selectively the content of the `src/init.sh` script.

### Docker installation

With the Docker installation every step of the pipeline is launched within a Docker container. This is the recommended installation for your personal computer as you don't have to manage the different version and install of the pipeline tools.

To install nextflow with docker you need to run the following commands, which require Java 7 or higher to be installed on your computer:

```sh
mkdir bin
cd bin
wget -qO- get.nextflow.io | bash
cd ..
```

Then you can build the pipeline docker image with the following command:

```sh
docker build src/func/docker -t 'pipeline:0.0.1'
```

### local installation

With the local installation, every step of the pipeline is launched locally. For this, you need to install and manage all the tools used by the pipeline.

Follow the instruction of the `src/init.sh` script to install all the required tools.

### PSMN installation

With the PSMN installation, every step of the pipeline is launched as a `qsub` command. The nextflow process needs to stay alive for the duration of the pipeline. Thus, we advise you to launch the pipeline within a `tmux` instance.

To use the pipeline on the PSMN you need to be able to load the LBMC modules.
To do so follow the instructions on the [PSMN/modules repository](http://gitlab.biologie.ens-lyon.fr/PSMN/modules/tree/master).

Then you can use the following command to load the nextflow module.

```sh
module load nextflow/0.25.1
```

## Testing

You can test the pipeline with the scripts from the folder `tests`.

First get the test data set with the following commands:

```sh
git submodule init
git submodule update
```

Then you can launch

- `bash tests/pipeline_docker_test.sh` to test your docker installation
- `bash tests/pipeline_local_test.sh` to test your local installation
- `bash tests/pipeline_docker_test.sh` to test your PSMN installation

## Project personnalisation

You can modify the configuration file `src/pipeline.config`. If you are working with single-end sequencing data, don't forget to specify the `mean` and `sd` parameter for the size distribution of your library fragments length (not read size).

## Pipeline commands

To launch the pipeline you need to use the following command:

```sh
bin/nextflow src/pipeline.nf -c src/pipeline.config -profile <profile> --fastq "path_to_fastq_files" --todo "fastqc+cutadapt+urqt+fastqc+multiqc"
```

with :

- `<profile>`: `docker`, `local` or `sge` (on the PSMN).
- `path_to_fastq_files`: the path to your fastq files. You can use `*.fastq` for single-end or `*_R{1,2}.fastq` for paired-end data.
- The `--todo` argument specify the chain of steps you wan't the pipeline to take. Here we run `fastqc`, then `cutadapt` to remove the adaptor, then `urqt` to trim by quality then `fastq` on the trimmed data and finally `multiqc`.
- On the PSMN one the nextflow modules is loaded it is directly available so you don't need to specify `bin/nextflow`, only `nextflow`.

For the mapping and quantification step, you need to specify a reference (in the `fasta` format) and an annotation file (in `gff(3)` or `gtf` format).

```sh
bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "path_to_fastq_files" --fasta "path_to_fasta_files" --annot "path_to_annotation_files" --todo "fastqc+cutadapt+cutadapt+fastqc+multiqc+split_ref+kallisto"
```

Here are some example of the `--todo` parameter:

- `"fastqc+cutadapt+urqt+fastqc+multiqc"`
- `"fastqc+cutadapt+cutadapt+fastqc+multiqc"`
- `"fastqc+cutadapt+cutadapt+fastqc+multiqc+split_ref+kallisto"`
- `"fastqc+cutadapt+cutadapt+fastqc+multiqc+bowtie2+htseq"`
- `"fastqc+cutadapt+cutadapt+fastqc+multiqc+bowtie2+rsem"`

You can see more examples in the `tests` folder scripts.
