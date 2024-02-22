#!/usr/bin/bash
#SBATCH --job-name="suppa2"
#SBATCH --err=logs/suppa.%j.err
#SBATCH --output=logs/suppa.%j.out
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --array=0-24

module load python/3.6
module load r

SUPPA=~/.local/lib/python3.6/site-packages/SUPPA-2.3
IOE=../references/gencode_GRCh38.p13/SUPPA2_splicemap/modID.GRCh38.p13.gencode.v37.primary_assembly.annotation_SE_strict.ioe

pargrid=()
i=0
while read u
do
    pargrid[i]=$u
    i=$((i+1))
done < run_parameters.txt
run=${pargrid[${SLURM_ARRAY_TASK_ID}]}

TXDIR=../output/tximport/run_${run}
OUTDIR=../output/suppa2/run_${run}

mkdir -p ${OUTDIR}

# Modify count file format

echo "Reformatting count data..."
INFILE=${TXDIR}/counts.raw.csv
CTSFILE=${OUTDIR}/counts.raw.tsv

sed 's/\"\",//' ${INFILE} | sed 's/\"//g' | sed 's/,/\t/g' > ${CTSFILE}

# Compute PSI 

echo "Computing PSI..."

OUT_PREFIX=${OUTDIR}/SE.events
python ${SUPPA}/suppa.py psiPerEvent -e ${CTSFILE} -i ${IOE} -o ${OUT_PREFIX}

# Separate files 

echo "Splitting files..."

SAM_CTS=${OUTDIR}/sample_counts/
SAM_PSI=${OUTDIR}/sample_psi/

mkdir -p ${SAM_CTS} ${SAM_PSI}

Rscript --vanilla sepconditions.R ${CTSFILE} ${OUT_PREFIX}.psi ${SAM_CTS} ${SAM_PSI}

# Differential testing

echo "Testing differential splicing..."

python ${SUPPA}/suppa.py diffSplice \
        --input ${IOE} \
        --method empirical \
        --psi $(ls ${SAM_PSI}/*psi | sort ) \
        --tpm $(ls ${SAM_CTS}/*tsv | sort ) \
        -o ${OUT_PREFIX}

echo "Finished successfully!"