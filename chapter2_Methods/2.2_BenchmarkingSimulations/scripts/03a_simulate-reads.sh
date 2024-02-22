#!/bin/bash
#SBATCH --job-name=reads
#SBATCH --account=rrg-hsn
#SBATCH --err=logs/reads.%j.err
#SBATCH --output=logs/reads.%j.out
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=12
#SBATCH --time=03:00:00
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

# Generate fastq files 

RSEM_RES=../input/rsem/${sample%_rep*}_rep1 # use the RSEM model parameters of replicate 1 in all samples 
#OUTDIR=../output/rsem_readsimulator/run_${run}
OUTDIR=/scratch/lmoral7/hgthesis_lms/rsem_readsimulator/run_${run}
SIM_CTS=../output/countsimulator/${sample}/${run}_rsem.isoforms.results

mkdir -p ${OUTDIR}
echo "Generating sim fastqs from $SIM_CTS"
echo "Will write output to $OUTDIR"

rsem-simulate-reads ${REFNAME} \
                    ${RSEM_RES}/rsem.stat/rsem.model \
                    ${SIM_CTS} \
                    0 \
                    30000000 \
                    ${OUTDIR}/${sample} \
                    --seed 0
                    
echo "Compressing fastq files..."                    
fq1=${OUTDIR}/${sample}_1.fq
fq2=${OUTDIR}/${sample}_2.fq
gzip $fq1
gzip $fq2

# Align reads

#FSQDIR=../output/rsem_readsimulator/run_${run}
#BAMDIR=../output/star/run_${run}

#mkdir -p ${BAMDIR} 

#bam=${BAMDIR}/${sample}
#STAR --runMode alignReads \
#     --genomeLoad  NoSharedMemory \
#     --outSAMtype BAM Unsorted \
#     --genomeDir ${GENDIR} \
#     --readFilesCommand zcat \
#     --readFilesIn ${fq1}.gz ${fq2}.gz \
#     --runThreadN 11 \
#     --outFileNamePrefix ${bam}. 

# Sort bamfile by position 

#echo "Sorting bamfile by position and converting into samfile"

#module load StdEnv/2020
#module load samtools/1.17

#bam=${BAMDIR}/${sample}.Aligned.out.bam
#sam=${BAMDIR}/${sample}.Aligned.out.pos.sorted.sam

#samtools sort -u -O sam -o ${sam} -@ 12 -m 2G ${bam}

#echo "Finished successfully!"
