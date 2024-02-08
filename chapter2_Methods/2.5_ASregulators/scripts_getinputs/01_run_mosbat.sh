#!/usr/bin/bash
#SBATCH --account=rrg-hsn
#SBATCH --job-name="mosbat"
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --err=logs/err.mbat.%j.log
#SBATCH --output=logs/out.mbat.%j.log
#SBATCH --mem=12G
#SBATCH --account=rrg-hsn

cd ~/lib/MoSBAT # Needs to run inside the directory where MoSBAT is installed

CISBP_FILE=${PWD%scripts_getinputs}/input/cisbp_Homo_sapiens_2023_06/all_cisbp_mosbat.txt
OUT_DIR=${PWD%scripts_getinputs}/output/mosbat
ID=Homo_sapiens_2023_06

mkdir -p ${OUT_DIR}/${ID}
bash MoSBAT.sh ${ID} ${CISBP_FILE} ${CISBP_FILE} rna 50 50000

cp -R out/${ID} ${OUT_DIR}/