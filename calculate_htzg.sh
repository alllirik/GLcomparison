#!/bin/sh
#
#SBATCH --job-name=htzg_recal
#SBATCH --output=run_htzg_recal.log
#SBATCH --error=error_htzg_recal.err
#SBATCH --ntasks=1

module load Python/Anaconda_v10.2019
source activate project
module load gnu10
module load armadillo

#ATLAS: PMD and GLF for each bam file
for BAM in bams/low_ZVEJ31*.bam
do
        #samtools index $BAM
        bamname=${BAM##*/}
        newbam=${bamname%%.*}_recal
        #IF ANCIENT [fasta should be in script location file]
        #srun ~/atlas/./atlas task=PMD bam=$BAM fasta=human_g1k_v37.fasta chr=20 out=PMD/${bamname%%.*} length=20
        srun ~/atlas/./atlas task=recal bam=$BAM chr=chrX haploid=chrX minDepth=2 out=recal/${newbam}
        srun ~/atlas/./atlas task=GLF bam=$BAM chr=chr20 recal=recal/${newbam}_recalibrationEM.txt out=GLF2_chr20/${newbam}
        srun ~/atlas/./atlas task=printGLF glf=GLF2_chr20/${newbam}.glf.gz | tail -n +17 | head -n -3 | awk '{$3=$4=""; print}'> GLF2_chr20/${newbam}.glf
done

#ATLAS GLF format to ANGSD GLF format

runname=low_ZVEJ31_recal

#ONLY FOR FILES WITH 'CHR' INSTEAD OF NUMBER
awk '{$1 = 20; print}' GLF2_chr20/low_ZVEJ31*softclip_recal.glf > ${runname}_chr20_atlas.glf
#Log-scaled likelihood and add one to every position
awk '{printf $1"\t"$2+1"\t"; for ( i=3; i<=NF; i++ ) {printf $i*-0.0023025850"\t" } }{ print""}' ${runname}_chr20_atlas.glf > ${runname}_chr20_atlas_log.glf

source activate project_secondary

srun angsd -glf10_text ${runname}_chr20_atlas_log.glf -nInd 1 -doSaf 1 -fai human_g1k_v37.fasta.fai -anc ~/chimpHg19.fa.gz -out  SAF/${runname}_chr20_atlas -P 4

srun realSFS  SAF/${runname}_chr20_atlas.saf.idx -maxIter 100 -P 4 >  SFS/${runname}_chr20_atlas.sfs
