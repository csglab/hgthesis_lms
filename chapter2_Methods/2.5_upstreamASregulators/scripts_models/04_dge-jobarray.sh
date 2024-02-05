#!/usr/bin/bash
#SBATCH --job-name="dge"
#SBATCH --err=logs/dge.%j.err
#SBATCH --output=outs/dge.%j.out
#SBATCH --account=rrg-hsn
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --time=4:00:00
#SBATCH --array=0-19

# Load modules
module load StdEnv/2020
module load r/4.1.2 

# Define variables
args=()
i=0
while read p
do
  args[i]="$p"
  i=$((i+1))      
done < cancers_valid_condition.tsv

cancer=${args[${SLURM_ARRAY_TASK_ID}]}

# Run script

echo "Input args:"
echo " - Cancer = $cancer"

INOBJ=../input/gene_expression/tcga.gene.expression.inputs.RData
OUTDIR=../output/rbp_dge/${cancer}

mkdir -p ${OUTDIR}

Rscript dge_analysis.r -c ${cancer} -i ${INOBJ} -o ${OUTDIR}
