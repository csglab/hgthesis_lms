#!/bin/bash
#SBATCH --job-name="sample-control"
#SBATCH --account=rrg-hsn
#SBATCH --err=logs/sample.%j.err
#SBATCH --output=logs/sample%j.out
#SBATCH --mem=12gb
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=larimst@gmail.com
#SBATCH --mail-type=ALL

# Inputs 

SAMPLE=K562_control_rep1 # Sample file to be used as a template
RSEM_RESDIR=../input/rsem/${SAMPLE}
INDIR=../input
OUTDIR=../output/countsimulator

# Modules

module load r/4.2.1

# Main program

mkdir -p ${OUTDIR}

Rscript simulate-counts.r ${INDIR} ${OUTDIR} ${RSEM_RESDIR}/rsem.isoforms.results 

