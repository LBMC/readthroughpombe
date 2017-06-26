# Quality control

The first part of the pipeline is about quality control of the data.
For this we need 4 steps:
- check the quality of the fastq files
- remove (if any) the adaptor of the reads
- trim the reads to remove nucleotid of poor quality
- check the quality of the results

Those steps are executed with the `pipe/quality_control.nf` script

There are two testing script regarding quality control in the `tests` directory.

You can either run:

```
sh tests/quality_control_test.sh
```
To run the analysis locally. You will need to install every required tools
and configure the pipeline to find them in the respective path.

Otherise, you can run:

```
sh tests/quality_control_docker_test.sh
```
To run the analysis in a Docker container, the tools will be installed
automatically within the container and the default paths will work.
