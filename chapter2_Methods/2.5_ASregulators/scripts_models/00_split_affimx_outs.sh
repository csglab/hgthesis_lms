#events=("SE" "MX" "RI" "A5" "A3" "AL" "AF")
events=("SE" "A5")
AFFDIR=../output/affimx_by_motif
OUTDIR=../output/affimx

readarray -t motifs < 25112023_representative_motif_ids.txt

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
