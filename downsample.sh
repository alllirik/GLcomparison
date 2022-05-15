#!/bin/sh
#
#SBATCH --job-name=downsample
#SBATCH --output=ds3.log
#SBATCH --error=ds3.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

module load Python/Anaconda_v10.2019
source activate project

srun samtools view -bs 42.3 17X.ZVEJ25.chr20.bam > 5X.ZVEJ25.chr20.bam
srun samtools index 5X.ZVEJ25.chr20.bam

srun samtools view -bs 42.24 17X.ZVEJ25.chr20.bam > 4X.ZVEJ25.chr20.bam
srun samtools index 4X.ZVEJ25.chr20.bam

srun samtools view -bs 42.18 17X.ZVEJ25.chr20.bam > 3X.ZVEJ25.chr20.bam
srun samtools index 3X.ZVEJ25.chr20.bam

srun samtools view -bs 42.12 17X.ZVEJ25.chr20.bam > 2X.ZVEJ25.chr20.bam
srun samtools index 2X.ZVEJ25.chr20.bam

srun samtools view -bs 42.06 17X.ZVEJ25.chr20.bam > 1X.ZVEJ25.chr20.bam
srun samtools index 1X.ZVEJ25.chr20.bam



