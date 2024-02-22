#!/usr/bin/bash
#SBATCH --job-name="trex"
#SBATCH --cpus-per-task=20
#SBATCH --time=4:00:00
#SBATCH --err=logs/trex.%j.err
#SBATCH --output=logs/trex.%j.out
#SBATCH --mem=40G
#SBATCH --array=0-24
#SBATCH --account=rrg-hsn

#----------------------------------------------------------#
echo "The job "${SLURM_JOB_ID}" is running on "${SLURM_JOB_NODELIST}
#----------------------------------------------------------#

# Modules required

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3
module load salmon/1.7.0

# Variables

FAIDX=../references/gencode_GRCh38.p13/salmon_genomeIndex
SPLICEMAP=../references/gencode_GRCh38.p13/SUPPA2_splicemap/modID.GRCh38.p13.gencode.v37.primary_assembly.annotation_ASALL_strict.ioe
METADATA=../data/metadata.tsv
EVENT=SE
TXIFORMAT=tximeta
declare -a reps=("rep1" "rep2" "rep3" "rep4" "rep5")

# Parameter array

pargrid=()
i=0
while read u
do
	pargrid[i]=$u
	i=$((i+1))
done < run_parameters.txt
run=${pargrid[${SLURM_ARRAY_TASK_ID}]}

# Run salmon
# =================================

FSQDIR=../output/rsem_readsimulator/run_${run}
OUTDIR=../output/salmon/run_${run}
S1=${FSQDIR}/K562_control
S2=${FSQDIR}/K562_SRSF9.KD

mkdir -p ${OUTDIR}

echo "Running salmon..."
for rep in "${reps[@]}"
do
    OUTSAM=${OUTDIR}/K562_control_${rep}
    mkdir -p ${OUTSAM}
    salmon quant -p 19 \
                 -i ${FAIDX} \
                 -l A --gcBias \
                 -1 ${S1}_${rep}_1.fq.gz \
                 -2 ${S1}_${rep}_2.fq.gz \
                 --validateMappings \
                 -o ${OUTSAM} 

    OUTSAM=${OUTDIR}/K562_SRSF9.KD_${rep}
    mkdir -p ${OUTSAM}
    salmon quant -p 19 \
                 -i ${FAIDX} \
                 -l A --gcBias \
                 -1 ${S2}_${rep}_1.fq.gz  \
                 -2 ${S2}_${rep}_2.fq.gz  \
                 --validateMappings \
                 -o ${OUTSAM} 
done

# Load with tximport
# =================================

INDIR=../output/salmon/run_${run}
OUTDIR=../output/tximport/run_${run}

mkdir -p ${OUTDIR}

echo "Loading salmon outputs with tximport..."
Rscript build-count-matrix-sims.R ${INDIR} ${OUTDIR} ${METADATA}


# Run TRex 
# =================================

COUNTS=../output/tximport/run_${run}/txiobject.raw.RData
OUTDIR=../output/trex/run_${run}

mkdir -p ${OUTDIR}

# Get event counts
echo "Estimating event counts..."
Rscript trex_general-event_counts.R -s ${run} \
                                     -e ${EVENT} \
                                     -t ${TXIFORMAT} \
                                     -c ${COUNTS} \
                                     -m ${METADATA} \
                                     -r ${SPLICEMAP} \
                                     -o ${OUTDIR}

# Run TRex model fit
echo "Running TRex model fits..."
Rscript trex_general_batch.r -s ${run} \
                             -e ${EVENT} \
                             -d ${OUTDIR} \
                             -r ${SPLICEMAP}
 
echo "Finished successfully!"