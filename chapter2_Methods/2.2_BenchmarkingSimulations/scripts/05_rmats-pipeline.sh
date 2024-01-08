#SBATCH --err=logs/rmats.%j.err
#SBATCH --output=logs/rmats.%j.out
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --array=0-24

# Modules required

module load StdEnv/2020 gcc/9.3.0 
module load star
module load rmats/4.1.1

# Variables

G1=K562_control
G2=K562_SRSF9.KD
STAR=../references/gencode_GRCh38.p13/STAR_genomeIndex
GTF=../references/gencode_GRCh38.p13/gencode.v37.primary_assembly.annotation.gtf
TMP=tmp_output

params=()
i=0
while read p
do
	params[i]=$p
	i=$((i+1))
done < run_parameters.txt

run=${params[${SLURM_ARRAY_TASK_ID}]}

OUT_DIR=../output/rmats/run_${run}
PATH1=${OUT_DIR}/${G1}_bams.txt # files with paths to inpute bamfiles of each group bamfiles 
PATH2=${OUT_DIR}/${G2}_bams.txt

# Main

echo "Running rMATS..."
python $EBROOTRMATS/rmats.py \
        --b1 ${PATH1} \
        --b2 ${PATH2} \
        --gtf ${GTF} \
        --readLength 100 \
        --anchorLength 2 \
        --nthread 7 \
        --cstat 0 \
        --od ${OUT_DIR} --tmp ${TMP}

echo "Converting identifiers to SUPPA2 IDs..."
RMATSIN=../output/run_${run}/rmats/SE.MATS.JC.txt
RMATSOUT=../output/run_${run}/rmats/SE.MATS.JC.suppa2.tsv

./rmats_to_suppa_ids_SE_events.pl ${RMATSIN} > ${RMATSOUT}.tmp # Custom script from SUPPA2s github
head -1 ${RMATSIN} | sed 's/^ID/SUPPA2_ID\tID/' > headers.txt 
cat headers.txt ${RMATSOUT}.tmp > ${RMATSOUT} && rm ${RMATSOUT}.tmp  

echo "Finished successfully!"
