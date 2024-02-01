#!/usr/bin/bash
#SBATCH --job-name="flag"
#SBATCH --err=logs/flag.%j.err
#SBATCH --output=outs/flag.%j.out
#SBATCH --account=rrg-hsn
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#SBATCH --time=01:00:00
#SBATCH --array=0-139

# Load modules

module load StdEnv/2020 r/4.1.2

params=()
i=0
while read p
do
  params[i]="$p"
  i=$((i+1))      
done < flag_params.tsv

CANCER=$(echo ${params[${SLURM_ARRAY_TASK_ID}]} | cut -d "-" -f1 )
EVENT=$(echo ${params[${SLURM_ARRAY_TASK_ID}]} | cut -d "-" -f2 )

CANCER=STAD
EVENT=AF
TREXDIR=../../2.1_TRexMethod/tcga/output/${CANCER}/trex
OUTDIR=../output/PSI_flags/${CANCER}

mkdir -p ${OUTDIR}

echo "Input settings:"
echo "- Cancer = $CANCER"
echo "- Events = $EVENT"

Rscript get_psi_flags.r ${CANCER} ${EVENT} ${TREXDIR} ${OUTDIR}
