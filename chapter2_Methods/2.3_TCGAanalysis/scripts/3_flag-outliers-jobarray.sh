#!/usr/bin/bash

#SBATCH --job-name="outflag"
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --err=tmp/outflag.%j.err
#SBATCH --output=tmp/outflag.%j.out
#SBATCH --mem-per-cpu=24G
#SBATCH --array=0-30

# Load required modules

module load StdEnv/2020 r/4.1.2

# Generate tximport objects
###########################

cancers=()
i=0
while read p
do
  cancers[i]=$p
  i=$((i+1))
done < cancers_valid.tsv

CANCER=${cancers[${SLURM_ARRAY_TASK_ID}]}
TXFILE=../../2.1_TRexMethod/tcga/output/${CANCER}/tximport/${CANCER}.txiobject.RData
MTFILE=../output/processed_metadata/${CANCER}.csv
OUTDIR=../output/outlier_tests

mkdir -p ${OUTDIR}

Rscript flag_outliers.R  ${CANCER} ${TXFILE} ${MTFILE} ${OUTDIR}/${CANCER}.pdf 