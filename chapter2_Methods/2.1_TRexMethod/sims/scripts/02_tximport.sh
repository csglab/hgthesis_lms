#!/usr/bin/bash
#SBATCH --job-name="tximport"
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --err=err.counts.%j.log
#SBATCH --output=out.counts.%j.log
#SBATCH --mem-per-cpu=8G
#SBATCH --account=rrg-hsn
#SBATCH --array=0-4

# Inputs

PARS=run_parameters.txt
INDIR=../output/run_${run}/salmon
MTDAT=../input/metadata.tsv

# Modules
module load StdEnv/2020 r/4.2.1

# Build jobarray
pargrid=()
i=0
while read u
do
	pargrid[i]=$u
	i=$((i+1))
done < ${PARS}

run=${pargrid[${SLURM_ARRAY_TASK_ID}]}
OUTDIR=../output/run_${run}/tximport

mkdir -p ${OUTDIR}

# Run tximport
Rscript build-count-matrix-sims.r ${INDIR} ${OUTDIR} ${MTDAT}




