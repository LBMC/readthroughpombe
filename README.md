
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
