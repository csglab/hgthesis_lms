#!/usr/bin/bash
#SBATCH --account=rrg-hsn
#SBATCH --job-name="salmon_batch"
#SBATCH --cpus-per-task=20
#SBATCH --time=4:00:00
#SBATCH --err=salmon.%j.err
#SBATCH --output=salmon.%j.out
#SBATCH --mem=40G
#SBATCH --array=0-24

# Inputs
FAIDX=../references/gencode_GRCh38.p13/salmon_genomeIndex # salmon index 
FSQDIR=../input/fastqs # directory with paired fastq files for all samples
G1=K562_control # name of group 1
G2=K562_SRSF9.KD # name of group 2
reps=("rep1" "rep2" "rep3" "rep4" "rep5") # list of replicates per group
S1=${FSQDIR}/${G1} # Directory with replicates of group 1
S2=${FSQDIR}/${G2} # Directory with replicates of group 2
PARS=run_parameters.txt

# Modules 
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3
module load salmon/1.7.0

# Build jobarray 
pargrid=()
i=0
while read u
do
    pargrid[i]=$u
    i=$((i+1))
done < ${PARS}
run=${pargrid[${SLURM_ARRAY_TASK_ID}]}
OUTDIR=../output/run_${run}/salmon

mkdir -p ${OUTDIR}

for rep in "${reps[@]}"
do
	echo "Quantifying replicate ${rep}"
	echo "==========================================="

	OUTSAM=${OUTDIR}/${G1}_${rep}
	mkdir -p ${OUTSAM}
	salmon quant -p 19 \
				 -i ${FAIDX} \
				 -l A --gcBias \
				 -1 ${S1}_${rep}_1.fq \
				 -2 ${S1}_${rep}_2.fq \
				 --validateMappings \
				 -o ${OUTSAM} 

	OUTSAM=${OUTDIR}/${G2}_${rep}
	mkdir -p ${OUTSAM}
	salmon quant -p 19 \
				 -i ${FAIDX} \
				 -l A --gcBias \
				 -1 ${S2}_${rep}_1.fq  \
				 -2 ${S2}_${rep}_2.fq  \
				 --validateMappings \
				 -o ${OUTSAM} 
done




 
