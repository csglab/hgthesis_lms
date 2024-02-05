#!/usr/bin/bash
#SBATCH --job-name="lmfit"
#SBATCH --err=logs/lmfit.%j.err
#SBATCH --output=outs/lmfit.%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=60gb
#SBATCH --time=5:00:00
#SBATCH --array=0-41
#SBATCH --account=rrg-hsn

# Load modules
module load StdEnv/2020
module load r/4.1.2 

# Define variables
WS_ARR=("50" "100" "200")
COEF_ARR=("res.lfcShrink")
EVENT_ARR=("SE" "AF" "MX" "RI" "A5" "A3" "AL")
CANCER_ARR=("BLCA" "BRCA" "CESC" "CHOL" "COAD" "ESCA" "HNSC" "KICH" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "PAAD" "PCPG" "PRAD" "READ" "STAD" "THCA" "UCEC")
#CANCER_ARR=("COAD" "HNSC" "LIHC" "PRAD" "READ" "STAD")

args=()
i=0
for cancer in ${CANCER_ARR[@]}
do
    for w in ${WS_ARR[@]}
    do
        for c in ${COEF_ARR[@]}
        do
            for e in ${EVENT_ARR[@]}
            do
                str=${cancer}-${w}-${c}-${e}
                args[i]=${str}
                (( i++ ))
            done
        done
    done
done

arg=${args[${SLURM_ARRAY_TASK_ID}]}
CANCER=$(echo $arg | cut -d "-" -f1)
WS=$(echo $arg | cut -d "-" -f2)
COEF=$(echo $arg | cut -d "-" -f3)
EVENT=$(echo $arg | cut -d "-" -f4)

echo "Input args:"
echo " - Cancer = $CANCER"
echo " - Window size = $WS"
echo " - Coefficient = $COEF"
echo " - Event type = $EVENT"

# Run scripts

TREXOBJ=../input/trex_objects/${CANCER}.condition.${COEF}.conditiontumor.RDS
OUTDIR=../output/rbp_models

mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/coefficients ${OUTDIR}/training_coefficients ${OUTDIR}/performance_reports 

Rscript fit_RBP_models.r --windowSize ${WS} --cancer ${CANCER} --trexObject ${TREXOBJ} --eventType ${EVENT} --outDir ${OUTDIR}
