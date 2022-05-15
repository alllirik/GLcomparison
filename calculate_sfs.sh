#!/bin/sh
#
#SBATCH --job-name=atlas
#SBATCH --output=atlas.log
#SBATCH --error=atlas.err
#SBATCH --ntasks=1

module load Python/Anaconda_v10.2019
source activate project
module load gnu10
module load armadillo

#COV=2X

#Simulate BAM files
#~/atlas/./atlas task=simulate out=${COV} chrLength=1000000 depth=2 writeTrueGenotypes writeTrueAlleleFreq theta=0.005 type=SFS fixedSeed=10 sampleSize=10

#samtools faidx ${COV}.fasta
#mv ${COV}_trueSFS_chr1.txt ${COV}.truth.sfs

#Script for GL calculation using ATLAS recal
#Using script for merging GLF files into ANGSD format
#Prepare index files for bams
#MODERN SAMPLES ONLY (NO PMD)	
#Launch from parent directory

#Calculate recal files? (0/1)
recal_ans=0

#Clear/Recal? (0/1)
launch_mode=0

#File name (COVERAGE)
runname=HG0

#Number of samples?
nInd=10

#	SPECIFY BAMS WITH CORRECT COVERAGE IN THE NEXT LINE
for BAM in ~/BAMs/${runname}*.bam
do
        bamname=${BAM##*/}
        newbam=${bamname%%.*}_clear
	if [ $recal_ans == 1 ]; then
		srun ~/atlas/./atlas task=recal bam=$BAM haploid=X minDepth=2 out=${newbam}
	fi
	if [ $launch_mode == 0 ]; then
        	srun ~/atlas/./atlas task=GLF chr=20 bam=$BAM out=${newbam}
        	srun ~/atlas/./atlas task=printGLF glf=${newbam}.glf.gz | tail -n +17 | head -n -3 | awk '{$3=$4=""; print}'> ${newbam}.glf
	else
		srun ~/atlas/./atlas task=GLF bam=$BAM recal=${newbam}_recalibrationEM.txt out=${newbam}
		srun ~/atlas/./atlas task=printGLF glf=${newbam}.glf.gz | tail -n +17 | head -n -3 | awk '{$3=$4=""; print}'> ${newbam}.glf
	fi
done


#ATLAS GLF format to ANGSD GLF format [CHANGE NAME OF FILES IF NECESSARY]
ls ${runname}*_clear.glf > ${runname}.filelist


file1=$(sed '1q;d' ${runname}.filelist)
file2=$(sed '2q;d' ${runname}.filelist)

N=$(awk '{print NF}' $file1 | head -n 1)
flds=$(eval echo 1.{3..$N})
join -e 0  -a1 -a2 -1 2 -2 2 -o 1.1 0 $flds 2.{3..12} $file1 $file2 | awk '{$1 = 20; print}' > ${runname}_chr20_atlas_merged.glf

if [ $(wc -l < ${runname}.filelist) -ge 3 ]; then
        for file2 in $(tail -n +3 ${runname}.filelist)
        do
                file1=${runname}_chr20_atlas_merged.glf

                N=$(awk '{print NF}' $file1 | head -n 1)
                flds=$(eval echo 1.{3..$N})
                join -e 0  -a1 -a2 -1 2 -2 2 -o 1.1 0 $flds 2.{3..12} $file1 $file2 | awk '{$1 = 20; print}' > atlas_merged_tmp.glf
                cp atlas_merged_tmp.glf ${runname}_chr20_atlas_merged.glf
                rm atlas_merged_tmp.glf
        done
fi

#At this point all the files are merged into "{runname}_chr20_atlas_merged.glf"

#Log-scaled likelihood and add one to every position
awk '{printf $1"\t"$2+1"\t"; for ( i=3; i<=NF; i++ ) {printf $i*-0.0023025850"\t" } }{ print""}' ${runname}_chr20_atlas_merged.glf > ${runname}_chr20_atlas_log.glf

#samtools faidx ${COV}.fasta

#Now calculating SFS using ANGSD
source activate project_secondary

srun angsd -glf10_text ${runname}_chr20_atlas_log.glf -nInd ${nInd} -doSaf 1 -fai human_g1k_v37.fasta.fai -anc ~/chimpHg19.fa.gz -out ${runname}_atlas -P 4
srun realSFS ${runname}_atlas.saf.idx -maxIter 100 -P 4 > ${runname}.atlas.sfs

ls ~/BAMs/${runname}*.bam > bam.filelist


#GL1
srun angsd -bam bam.filelist -GL 1 -nInd 10 -doSaf 1 -anc ~/chimpHg19.fa.gz -out ${runname}.samtools -P 4
srun realSFS ${runname}.samtools.saf.idx -maxIter 100 -P 4 > ${runname}.samtools.sfs
#GL2
srun angsd -bam bam.filelist -GL 2 -nInd 10 -doSaf 1 -anc ~/chimpHg19.fa.gz -out  ${runname}.gatk -P 4
srun realSFS ${runname}.gatk.saf.idx -maxIter 100 -P 4 > ${runname}.gatk.sfs

#Delete useless files
rm *.saf*
rm *.arg
#rm *.glf*

#source activate project
#python plot_sfs.py -c ${COV}
