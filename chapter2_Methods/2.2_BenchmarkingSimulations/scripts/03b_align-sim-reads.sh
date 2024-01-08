#!/bin/bash
#SBATCH --job-name="star_align"
#SBATCH --err=logs/star.%j.err
#SBATCH --output=logs/star.%j.out
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=12
#SBATCH --time=00:45:00
#SBATCH --array=0-49

# Modules required

module load nixpkgs/16.09
module load gcc/7.3.0
module load star
module load samtools/1.17

GENDIR=../references/gencode_GRCh38.p13/STAR_genomeIndex

pargrid=()
i=0
while read u
do
	pargrid[i]=$u
	i=$((i+1))
done < parameters.txt

pars=${pargrid[${SLURM_ARRAY_TASK_ID}]}
sample=$(echo ${pars} | cut -d " " -f1)
run=$(echo ${pars} | cut -d " " -f2)

FSQDIR=../output/rsem_readsimulator/run_${run}
BAMDIR=../output/star/run_${run}

mkdir -p ${BAMDIR} 

bam=${BAMDIR}/${sample}
fq1=${FSQDIR}/${sample}_1.fq
fq2=${FSQDIR}/${sample}_2.fq

STAR --runMode alignReads \
     --genomeLoad  NoSharedMemory \
     --outSAMtype BAM Unsorted \
     --genomeDir ${GENDIR} \
     --readFilesIn ${fq1} ${fq2} \
     --runThreadN 11 \
     --outFileNamePrefix ${bam}. 

echo "Finished this alignment!"

echo "Sorting bamfile by position and converting into samfile"

bam=${BAMDIR}/${sample}.Aligned.out.bam
sam=${BAMDIR}/${sample}.Aligned.out.pos.sorted.sam

samtools sort -u -O sam -o ${sam} -@ 12 -m 2G ${bam}

echo "Finished sorting!"
