#!/usr/bin/bash

#SBATCH --job-name="txi"
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --err=logs/txi.%j.err
#SBATCH --output=outs/txi.%j.out
#SBATCH --mem-per-cpu=24G
#SBATCH --account=rrg-hsn
#SBATCH --array=0-30

# Input 

CANCERS=cancers_valid.tsv

# Modules

module load StdEnv/2020 r/4.2.1

# Build jobarray

cancers=()
i=0
while read p
do
  cancers[i]=$p
  i=$((i+1))
done < ${CANCERS}
CANCER=${cancers[${SLURM_ARRAY_TASK_ID}]}
OUTDIR=../results/${CANCER}/tximport

mkdir -p ${OUTDIR}

# Run tximport

echo "Processing sample ${CANCER}"
INDIR=../output/${CANCER}/salmon
MTDAT=../input/${CANCER}/metadata.csv

# Transcript counts
Rscript build-count-matrix.r ${INDIR} ${OUTDIR} ${MTDAT} ${CANCER}

# Gene counts
Rscript build-gene-count-matrix.r ${INDIR} ${OUTDIR} ${MTDAT} ${CANCER}