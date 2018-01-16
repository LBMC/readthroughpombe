#$ -S /bin/bash
### nom du job:
#$ -N preprocess_MUSIC
### file d'attente:
#$ -q E5-2670deb128*
### parallel environnement & nslots
#$ -pe openmp16 16
### charger l'environnement utilisateur pour SGE
#$ -cwd
### exporte les variables d'environnement sur les noeuds d'exécution
#$ -V
### change logs folder
#$ -o /home/cburny/logs
#$ -e /home/cburny/logs


source /applis/PSMN/modeles/Bashrc
source /usr/share/modules/init/bash
module use /applis/PSMN/Modules
module load Base/psmn
module load use.own
module load MUSIC/6613c53
module load SAMtools/1.5

##### Preprocess

# cut14
samtools view /scratch/cburny/Input_Music/cut14_forward/2017_10_26_PLBD2_cut14-208_R1_trim_sort_forward_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_cut14_forward/ &&\

samtools view /scratch/cburny/Input_Music/cut14_forward/2017_10_26_PLBD7_cut14-208_R2_trim_sort_forward_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_cut14_forward/ &&\

samtools view /scratch/cburny/Input_Music/cut14_forward/2017_10_26_PLBD12_cut14-208_R3_trim_sort_forward_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_cut14_forward/ &&\


samtools view /scratch/cburny/Input_Music/cut14_reverse/2017_10_26_PLBD2_cut14-208_R1_trim_sort_reverse_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_cut14_reverse/ &&\

samtools view /scratch/cburny/Input_Music/cut14_reverse/2017_10_26_PLBD7_cut14-208_R2_trim_sort_reverse_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_cut14_reverse/ &&\

samtools view /scratch/cburny/Input_Music/cut14_reverse/2017_10_26_PLBD12_cut14-208_R3_trim_sort_reverse_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_cut14_reverse/ &&\


# wt
samtools view /scratch/cburny/Input_Music/wt_forward/2017_10_26_PLBD1_wt_R1_trim_sort_forward_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_wt_forward/ &&\

samtools view /scratch/cburny/Input_Music/wt_forward/2017_10_26_PLBD11_wt_R3_trim_sort_forward_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_wt_forward/ &&\

samtools view /scratch/cburny/Input_Music/wt_forward/2017_10_26_PLBD6_wt_R2_trim_sort_forward_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_wt_forward/ &&\


samtools view /scratch/cburny/Input_Music/wt_reverse/2017_10_26_PLBD1_wt_R1_trim_sort_reverse_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_wt_reverse/ &&\

samtools view /scratch/cburny/Input_Music/wt_reverse/2017_10_26_PLBD11_wt_R3_trim_sort_reverse_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_wt_reverse/ &&\

samtools view /scratch/cburny/Input_Music/wt_reverse/2017_10_26_PLBD6_wt_R2_trim_sort_reverse_without.bam | MUSIC -preprocess SAM stdin /scratch/cburny/Output_Music/output_wt_reverse/ 


##### Sort

MUSIC -sort_reads /scratch/cburny/Output_Music/output_wt_reverse/ /scratch/cburny/Output_Music/output_wt_reverse/sorted &&\

MUSIC -sort_reads /scratch/cburny/Output_Music/output_wt_forward/ /scratch/cburny/Output_Music/output_wt_forward/sorted &&\


MUSIC -sort_reads /scratch/cburny/Output_Music/output_cut14_reverse/ /scratch/cburny/Output_Music/output_cut14_reverse/sorted &&\

MUSIC -sort_reads /scratch/cburny/Output_Music/output_cut14_forward/ /scratch/cburny/Output_Music/output_cut14_forward/sorted &&\


##### Remove duplicates

MUSIC -remove_duplicates /scratch/cburny/Output_Music/output_wt_reverse/sorted 2 /scratch/cburny/Output_Music/output_wt_reverse/dedup &&\

MUSIC -remove_duplicates /scratch/cburny/Output_Music/output_wt_forward/sorted 2 /scratch/cburny/Output_Music/output_wt_forward/dedup &&\


MUSIC -remove_duplicates /scratch/cburny/Output_Music/output_cut14_reverse/sorted 2 /scratch/cburny/Output_Music/output_cut14_reverse/dedup &&\

MUSIC -remove_duplicates /scratch/cburny/Output_Music/output_cut14_forward/sorted 2 /scratch/cburny/Output_Music/output_cut14_forward/dedup &&\
