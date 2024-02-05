#!/usr/bin/bash
#SBATCH --account=rrg-hsn
#SBATCH --job-name="fetch_seqs"
#SBATCH --cpus-per-task=5
#SBATCH --time=6:00:00
#SBATCH --err=logs/err.fetch.%j.log
#SBATCH --output=logs/out.fetch.%j.log
#SBATCH --mem=60G
#SBATCH --account=rrg-hsn

module load r/4.1.2

Rscript fetch_seqs.r
