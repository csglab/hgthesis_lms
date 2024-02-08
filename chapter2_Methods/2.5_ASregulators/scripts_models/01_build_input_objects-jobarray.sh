#!/usr/bin/bash
#SBATCH --job-name="getobjs"
#SBATCH --err=logs/lmfit.%j.err
#SBATCH --output=outs/lmfit.%j.out
#SBATCH --cpus-per-task=10
#SBATCH --mem=40gb
#SBATCH --time=2:00:00
#SBATCH --array=0-1
#SBATCH --account=rrg-hsn


# Load modules
module load StdEnv/2020
module load r/4.1.2 

# Define variables

#AFF_DIR=../output/affimx
AFF_DIR=~/lmprojects/splicing-pancancer/results/upstream_rbp_AffiMx_forward/by_event
WS_ARR=("100")
#WS_ARR=("50" "100" "200")
#EVENT_ARR=("SE" "AF" "MX" "RI" "A5" "A3" "AL")
EVENT_ARR=("SE" "A5")
COEF="res.lfcShrink"
args=()
i=0

for w in ${WS_ARR[@]}
do
    for e in ${EVENT_ARR[@]}
    do
        str=${w}-${e}
        args[i]=${str}
        (( i++ ))
    done
done

arg=${args[${SLURM_ARRAY_TASK_ID}]}
WS=$(echo $arg | cut -d "-" -f1)
EVENT=$(echo $arg | cut -d "-" -f2)

echo "Input args:"
echo " - Window size = $WS"
echo " - Event type = $EVENT"

# Run scripts

OUTDIR=../output/rbp_models

mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/objects

Rscript build_motif_data_objects.r --affDir ${AFF_DIR} --windowSize ${WS} --eventType ${EVENT} --outDir ${OUTDIR}
