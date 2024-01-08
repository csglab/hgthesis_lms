#!/usr/bin/bash
#SBATCH --job-name="refs"
#SBATCH --mem=64gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=6:00:00
#SBATCH --err=refs.%j.err
#SBATCH --output=refs.%j.out
#SBATCH --account=rrg-hsn

###########################################################################
# DESCRIPTION
# This job script generates two reference annotations:
# 1) Exon skipping events file with SUPPA2
# 2) Salmon genome index
# The scripts require as input three large files that were downloaded from:
# 1. GTF file of comprehensive gene annotation of primary assembly of 
# GRCh38.p13 https://www.gencodegenes.org/human/release_38.html
# 2. The gentrome.fa and decoys.txt files for salmon were generated as in
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
###########################################################################

# Inputs

REFDIR=../references/gencode_GRCh38.p13 # directory with all reference files
SUPPA=~/.local/lib/python3.9/site-packages/SUPPA-2.3 # Local installation of SUPPA2 

GTF=${REFDIR}/gencode.v37.primary_assembly.annotation.gtf # GTF file
GENTR=${REFDIR}/gentrome.fa  # gentrome fasta
DECOY=${REFDIR}/decoys.txt   # decoys for salmon

# Build genome index for Salmon
###################################################

# Modules
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 salmon/1.3.0

# Run salmon index
OUT=${REFDIR}/salmon_genomeIndex
salmon index -t ${GENTR} \
             -d ${DECOY} \
             -p 12 \
             -i ${OUT} \
             --tmpdir _indextmp \
             --gencode

# Build splice map reference annotation with SUPPA2
###################################################

# Modules
module load python/3.6

# Run generateEvents
OUT=${REFDIR}/gencode.v37.primary_assembly.annotation
python ${SUPPA}/suppa.py generateEvents \
    -i ${GTF} \
    -o ${OUT} \
    -f ioe \ 
    -e SE \
    --pool-genes 