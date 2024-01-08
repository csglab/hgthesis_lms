#!/bin/bash
#SBATCH --job-name="rsem"
#SBATCH --err=logs/rsem.%j.err
#SBATCH --output=logs/rsem.%j.out
#SBATCH --mem=36gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=3:00:00
#SBATCH --array=0-3

# Input
# Raw fastq files downloaded from https://www.encodeproject.org/experiments/ENCSR972AZD/
# to re run this steps download raw fastqs and place them  under input/fastqs

SAMPLEID=ENCSR972AZD_ENCSR341TTW
SAMPLES=samples.${SAMPLEID}.txt
REFDIR=../references/gencode_GRCh38.p13
REFNAME=${REFDIR}/RSEM_genomeIndex/ref

# Modules

module load StdEnv/2020
module load rsem/1.3.3
module load star

ulimit -v 31751044396

# Build job array 

samples=()
i=0
while read u
do
    samples[i]=$u
    i=$((i+1))
done < ${SAMPLES}

sample=${samples[${SLURM_ARRAY_TASK_ID}]}
OUTDIR=../input/rsem
FSQDIR=../input/fastqs 

mkdir -p ${OUTDIR} ${OUTDIR}/${sample}

# Run RSEM to measure expression

fq1=${FSQDIR}/${sample}_R1.fastq.gz
fq2=${FSQDIR}/${sample}_R2.fastq.gz

rsem-calculate-expression --paired-end \
                          --num-threads 11 \
                          --star \
                          --seed 7 \
                          --no-bam-output \
                          --star-gzipped-read-file \
                          ${fq1} \
                          ${fq2} \
                          ${REFNAME} \
                          ${OUTDIR}
