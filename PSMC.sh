#!/bin/sh
#
#SBATCH --job-name=psmc
#SBATCH --output=psmc.log
#SBATCH --error=psmc.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

module load Python/Anaconda_v10.2019
source activate project_secondary

NAME=LP6005441-DNA_F05.srt.aln.4X.clean

srun  angsd -i BAMs/hg38/${NAME}.bam -dopsmc 1  -out psmc_input/${NAME}.old -GL 2 -minQ 20 -minMapQ 30
srun ./ngsPSMC_original/ngsPSMC psmc_input/${NAME}.psmc.idx -p "1*4+25*2+1*4+1*6" -dospline 0 -nthreads 16 -nIter 20 -init 1 -theta 0.000233095 -rho 0.005357 > ${NAME}.psmc
