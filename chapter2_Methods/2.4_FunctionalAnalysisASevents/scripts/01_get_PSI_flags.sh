#!/usr/bin/bash
#SBATCH --job-name="flag"
#SBATCH --err=logs/flag.%j.err
#SBATCH --output=outs/flag.%j.out
#SBATCH --account=rrg-hsn
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --time=01:00:00
#SBATCH --array=0-85

# Load modules

module load StdEnv/2020 r/4.1.2

# Define variables

CTS_DIR=~/lmprojects/hgthesis_lms/chapter2_Methods/2.1_TRexMethod/tcga/output
OUT_DIR=../output/psi_flags

params=()
i=0
while read p
do
  params[i]="$p"
  i=$((i+1))      
done < flag_params.mis

cancer=$(echo ${params[${SLURM_ARRAY_TASK_ID}]} | cut -d "-" -f1 )
event=$(echo ${params[${SLURM_ARRAY_TASK_ID}]} | cut -d "-" -f2 )

mkdir -p ${OUT_DIR}

echo "Input settings:"
echo "- Cancer = $cancer"
echo "- Events = $event"

Rscript get_psi_flags.r ${cancer} ${event} ${CTS_DIR} ${OUT_DIR}
