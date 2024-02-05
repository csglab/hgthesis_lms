events=("SE" "MX" "RI" "A5" "A3" "AL" "AF")
AFFDIR=~/projects/rrg-hsn/ahcorcha/ahcorcha/Collaborations/C15_rbp_binding_splicing_larisa/data/05_processed_AffiMx/
OUTDIR=/home/lmoral7/lmprojects/splicing-pancancer/results/upstream_rbp_AffiMx_forward/by_event

readarray -t motifs < 25112023_representative_motif_ids.txt

#cd ../results/upstream_rbp_AffiMx_forward/by_motif
cd ${AFFDIR}

for event in ${events[@]}
do
    echo "Filtering events $event"
    for motif in ${motifs[@]}
    do
        for file in $(ls ${motif}*affinity.txt)
        do
            outfile=${OUTDIR}/${file%.txt}_${event}.txt
            if ! [ -f ${outfile} ]
            then
              echo "... $file"
              echo "File does not exist, copying..."
              grep "${event}" $file > ${outfile}
            fi
        done
    done
done
