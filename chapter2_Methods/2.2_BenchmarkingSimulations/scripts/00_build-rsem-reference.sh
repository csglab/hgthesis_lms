#!/bin/bash
#SBATCH --job-name="buildref"
#SBATCH --err=ref.%j.err
#SBATCH --output=ref.%j.out
#SBATCH --mem=70gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1:30:00

module load StdEnv/2020
module load star
module load rsem/1.3.3

REFDIR=../references/gencode_GRCh38.p13
GTF=${REFDIR}/gencode.v37.primary_assembly.annotation.gtf
GENOME=${REFDIR}/GRCh38.primary_assembly.genome.fa
REFNAME=${REFDIR}/RSEM_genomeIndex/ref

rsem-prepare-reference --gtf ${GTF} \
                       --star \
                       -p 11 \
                       ${GENOME} \
                       ${REFNAME}