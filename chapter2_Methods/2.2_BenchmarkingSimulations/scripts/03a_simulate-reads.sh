#!/bin/bash
#SBATCH --job-name="sim-rsem"
#SBATCH --account=rrg-hsn
#SBATCH --err=logs/sim.%j.err
#SBATCH --output=logs/sim.%j.out
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH --array=0-249

# Inputs

PARAMS=parameters.txt
REFDIR=../references/gencode_GRCh38.p13
REFNAME=${REFDIR}/RSEM_genomeIndex/ref

# Modules

module load perl/5.30.2 java/13.0.2 r/4.0.0 python/2.7.18

# Build jobarray
pargrid=()
i=0
while read u
do
	pargrid[i]=$u
	i=$((i+1))
done < ${PARAMS}

pars=${pargrid[${SLURM_ARRAY_TASK_ID}]}
sample=$(echo ${pars} | cut -d " " -f1)
run=$(echo ${pars} | cut -d " " -f2)

RSEM_RES=../input/rsem/${sample}
OUTDIR=../output/rsem_readsimulator/run_${run}
SIM_CTS=../output/countsimulator/${sample}/${run}_rsem.isoforms.results

mkdir -p ${OUTDIR}

rsem-simulate-reads ${REFNAME} \
                    ${RSEM_RES}/rsem.stat/rsem.model \
                    ${SIM_CTS} \
                    0 \
                    30000000 \
                    ${OUTDIR}/${sample} \
                    --seed 0

