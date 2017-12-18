#$ -S /bin/bash
### nom du job:
#$ -N test_MUSIC_cut14_wt
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

# You need to tun  this script in the directory where you want to store your peak calling outputs per analysis 

# Obtain window size per analysis
MUSIC -get_per_win_p_vals_vs_FC -chip /scratch/cburny/Output_Music/output_cut14_forward/dedup -control /scratch/cburny/Output_Music/output_wt_forward/dedup -l_win_step 50 -l_win_min 200 -l_win_max 16000 &&\

MUSIC -get_TF_peaks -chip /scratch/cburny/Output_Music/output_cut14_forward/dedup -control /scratch/cburny/Output_Music/output_wt_forward/dedup -mapp /scratch/cburny/Input_Music/temp -begin_l 100 -end_l 2000 -step 1.2 -l_mapp 50 -l_frag 363 -q_val 1 -l_p 0
