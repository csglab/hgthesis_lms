#!/usr/bin/bash
#SBATCH --job-name="trex"
#SBATCH --err=trex.%j.err
#SBATCH --output=trex.%j.out
#SBATCH --account=rrg-hsn
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --array=0-24

# Inputs
PARS=run_parameters.txt
SPLICEMAP=../references/gencode_GRCh38.p13/SUPPA2_splicemap/gencode.v37.primary_assembly.annotation_SE_strict.ioe
METADATA=../input/metadata.tsv

# Modules
module load r/4.2.1 

# Build jobarray
params=()
i=0
while read p
do
  params[i]=$p
  i=$((i+1))
done < ${PARS}
run=${params[${SLURM_ARRAY_TASK_ID}]}


# Run TRex
event=SE
TXIFORMAT=tximeta
COUNTS=../output/run_${run}/tximport/txiobject.raw.RData
OUTDIR=../output/run_${run}/trex

mkdir -p ${OUTDIR}

# Get TRex event counts
echo "run $run"
Rscript trex_general-event_counts.r -s ${run} \
                                     -e ${event} \
                                     -t ${TXIFORMAT} \
                                     -c ${COUNTS} \
                                     -m ${METADATA} \
                                     -r ${SPLICEMAP} \
                                     -o ${OUTDIR}

# Fit models
Rscript trex_general_sims.r -s ${run} \
                            -e ${event} \
                            -d ${OUTDIR} \
                            -r ${SPLICEMAP}