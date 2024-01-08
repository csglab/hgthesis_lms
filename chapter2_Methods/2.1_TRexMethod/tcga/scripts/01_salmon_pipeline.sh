#!/usr/bin/bash
#SBATCH --job-name="tcga_salmon"
#SBATCH --cpus-per-task=10
#SBATCH --time=3:00:00
#SBATCH --err=pipe_salmon.%j.err
#SBATCH --output=pipe_salmon.%j.out
#SBATCH --mem=40G
#SBATCH --array=0-10918

# Input files 

REF=../references/GRCh38/GRCh38.d1.vd1.fa # Genome fasta
REFBED=../references/GRCh38/hg38_RefSeq.bed # Genome RefSeq annotations BED file 
FAIDX=../references/gencode_GRCh38.p13/salmon_genomeIndex # Directory with salmon index
CRAMINFO=full_sample_info.txt # File with cancer, file_id and path to cram file

# Modules

module --force purge
module load StdEnv/2020 gcc/9.3.0 samtools/1.15.1

# Define global variables

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export WKDIR=${PWD}
export global_start_time=$SECONDS

# Build jobarray

pargrid=()
i=0
while read u
do
	pargrid[i]=$u
	i=$((i+1))
done < ${CRAMINFO} 
pars=${pargrid[${SLURM_ARRAY_TASK_ID}]}

# Retrieve inputs

CANCER=$(echo ${pars} | cut -d " " -f1)
FILEID=$(echo ${pars} | cut -d " " -f2)
CRAM=$(echo ${pars} | cut -d " " -f3)

export TMPDIR=~/scratch/tmp/${CANCER}

mkdir -p ${TMPDIR}

### CRAM to BAM sorted by name
###################################################

echo "------------------------------------"
echo "> Convert CRAM to namesorted BAM"

# Inputs

BAMDIR=../output/${CANCER}/namesorted_bamfiles
BAMFILE=${BAMDIR}/${FILEID}.bam

mkdir -p ${BAMDIR}

# Main

start_time=$SECONDS
if [ -s ${BAMFILE} ]
then 
    echo "Found name sorted bam, skipping to next stage."
else

    if [ -s ${CRAM} ]
    then
      echo "Converting CRAM file ..."
      samtools view -h -T ${REF} -b ${CRAM} -@ ${SLURM_CPUS_PER_TASK} | \
      samtools sort -o $BAMDIR/${FILEID}.bam -O bam -n -T $TMPDIR/${FILEID}.bam -@ ${SLURM_CPUS_PER_TASK} -

      elapsed=$(( SECONDS - start_time ))
      echo "Finished conversion and sorting in : $(date -ud "@$elapsed" +'%H hr %M min %S sec')"
    else
        if [ -s ${CRAM%.cram} ]
        then
              echo  "Found existing bamfile, sorting..."
              samtools sort -o $BAMDIR/${FILEID}.bam -O bam -n -T $TMPDIR/${FILEID}.bam -@ ${SLURM_CPUS_PER_TASK} ${CRAM%.cram}
              elapsed=$(( SECONDS - start_time ))
              echo "Finished sorting in : $(date -ud "@$elapsed" +'%H hr %M min %S sec')"
        else
              echo "File doesn't exist, need to re download. Exiting pipeline"
              echo $FILEID >> missing_cramfiles_to_redownload.txt
              exit
        fi
    fi
fi 

### QC bam files 
###################################################

echo "----------------------------------"
echo "> QC bam files"

# Inputs

STATSDIR=../output/${CANCER}/libstats
BAMSTAT=${STATSDIR}/summary_type_strand.txt
BAMLEN=${STATSDIR}/summary_read_length.txt

mkdir -p ${STATSDIR}

# Main

start_time=$SECONDS
bamcheck=$(samtools quickcheck -v "$BAMFILE")

if [ -s ${BAMLEN} ]
then
    echo "Found existing QC stats file, skipping to next stage."
else
    if [ -z "$bamcheck" ]
    then
        ./library-stats.sh ${REFBED} ${BAMFILE} ${TMPDIR}/${FILEID}.stats ${BAMSTAT} ${FILEID}
        ./read-length.sh ${BAMFILE} ${BAMLEN} ${FILEID}
        echo "Valid bamfile, proceding ..."
    else
      echo $FILEID >> corrupted_bamfiles.txt
      echo "Ivalid bamfile, writing to 'corrupted_bamfiles.txt' and exiting ..."
      exit
    fi
    elapsed=$(( SECONDS - start_time ))
    echo "---> Finished QC in : $(date -ud "@$elapsed" +'%H hr %M min %S sec')"
fi

### BAM to FASTQ
###################################################

echo "----------------------------------"
echo " > Convert BAM to FASTQ"

# Inputs

FQDIR=../output/${CANCER}/fastqs
FQ1=${FQDIR}/${FILEID}_R1.fastq
FQ2=${FQDIR}/${FILEID}_R2.fastq

mkdir -p ${FQDIR}

# Modules

module load bedtools
module load perl

# Main

if [ -s ${FQ2} ]
then
    echo "Found existing fastqs, skipping to next stage."
else
    if [ -z "$bamcheck" ]
    then
        echo "Converting bam file ..."
        start_time=$SECONDS
        
        FQS_TMP=${TMPDIR}/${FILEID}
        
        samtools collate -u -O ${BAMFILE} | \
        samtools fastq -n -F 2304 -1 ${FQS_TMP}_R1.fastq -2 ${FQS_TMP}_R2.fastq -0 ${FQS_TMP}_R0.fastq -s ${FQS_TMP}_Rs.fastq
        elapsed=$(( SECONDS - start_time ))
        echo "---> Finished bam to fastq conversion in : $(date -ud "@$elapsed" +'%H hr %M min %S sec')"
        
        echo "Shuffling FASTQ files ..."
        start_time=$SECONDS
        ./fastq-shuffle.pl -1 ${FQS_TMP}_R1.fastq -2 ${FQS_TMP}_R2.fastq
        mv ${FQS_TMP}_R1.fastq.shuffled ${FQ1}
        mv ${FQS_TMP}_R2.fastq.shuffled ${FQ2}
        elapsed=$(( SECONDS - start_time ))
        echo "---> Finished shuffling  in : $(date -ud "@$elapsed" +'%H hr %M min %S sec')"
        
    fi
fi

### Quantify using Salmon
###################################################

echo "-----------------------------"
echo "> Quantify with salmon"

# Inputs

OUTDIR=../output/${CANCER}/salmon
SALMONOUT=${OUTDIR}/${FILEID}

mkdir -p ${OUTDIR}

# Modules

module --force purge 
module load nixpkgs/16.09
module load gcc/7.3.0
module load openmpi/3.1.4
module load salmon/1.3.0 

# Main

start_time=$SECONDS
salmon quant -p ${SLURM_CPUS_PER_TASK} \
             -l A \
             -1 ${FQ1} \
             -2 ${FQ2} \
             -o ${SALMONOUT} \
             -i ${FAIDX} \
             --validateMappings \
             --gcBias
             
elapsed=$(( SECONDS - start_time ))
echo "---> Finished salmon  in : $(date -ud "@$elapsed" +'%H hr %M min %S sec')"


echo "Inspect salmon output"
FILENAME=${SALMONOUT}/quant.sf
if [ -f "${FILENAME}" ];then
    if [ -s "${FILENAME}" ];then
        echo "Corresponding quant.sf file exists and is not empty. Removing intermediate files"
        rm ${FQ1} ${FQ2} ${BAMFILE}
        rm ${TMPDIR}/${FILEID}*
    else
        echo "Corresponding quant.sf file exists but is empty"
        echo "Reanalyze ${FILEID}"
        exit
    fi
else
    echo "File ${FILENAME} not exists"
    echo "Reanalyze ${FILEID}"
    exit
fi


global_elapsed=$(( SECONDS - global_start_time ))
echo "################################################################### "
echo "Completed pipeline in: $(date -ud "@$global_elapsed" +'%H hr %M min %S sec')"
echo "################################################################### "
