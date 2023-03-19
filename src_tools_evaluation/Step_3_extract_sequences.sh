SAMPLE=$1
# TARGET="NEG_noncan_junctions"

# for OUTPUTFILE in "./dataset/pos" "./dataset/neg_1" "./dataset/neg_1_random"
for OUTPUTFILE in "./dataset/pos_refseq_protein_all" "./dataset/pos_refseq_protein_alts"
# for OUTPUTFILE in "./dataset/outlier_test"
do  
    for TARGET in "splam" "spliceai"
    do  
        if [[ $TARGET == "splam" ]]
        then
            echo "bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.donor.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.donor.seq.fa"

            bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.donor.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.donor.seq.fa

            echo "bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.acceptor.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.acceptor.seq.fa"

            bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.acceptor.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.acceptor.seq.fa
        fi

        if [[ $TARGET == "spliceai" ]]
        then
            echo "bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.seq.fa"

            bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/$TARGET/$TARGET.juncs.bed -fo $OUTPUTFILE/$TARGET/$TARGET.juncs.seq.fa
        fi
    done
done