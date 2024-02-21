#!/usr/bin/bash

#SBATCH --job-name="sash"
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --err=logs/sash.%j.err
#SBATCH --output=outs/sash.%j.out
#SBATCH --mem-per-cpu=20G
#SBATCH --account=rrg-hsn

module load python/3.11.5 r/4.1.2

GGSASHIMI=~/lib/ggsashimi/ggsashimi.py
GTF=../input/gencode.v37.primary_assembly.annotation.gtf
BAMS=figures_v3_data/KIRC_impurity_bins/CD46_input_bams.tsv # File with paths to bamfiles


# Impurity event

echo "Plotting impurity event"
${GGSASHIMI}  -b ${BAMS} \
              -c chr1:207790253-207790345 \
              -M 10 \
              -C 3 \
              -O 3 \
              --strand=SENSE \
              --out-strand=+\
              --alpha 0.25 \
              --fix-y-scale \
              -o impurity 

# Cancer event
echo "Plotting cancer event..."
${GGSASHIMI}  -b ${BAMS} \
              -c chr1:207767607-207767651 \
              -M 10 \
              -C 3 \
              -O 3 \
              --alpha 0.25 \
              --fix-y-scale \
              -o cancer
          