THRESHOLD=250
file='../Dataset/geneid_2_name.txt'
cat "$file" | jq -r '. | keys[]' | 
while IFS= read -r GENE_ID; do
    echo "GeneID: $GENE_ID"
    lines=$(samtools view ../Dataset/tranx_bam/$GENE_ID.bam | wc -l)
    echo $lines
    if [  $lines -gt $THRESHOLD ] ; then
        python parse_feature_parse_BAM.py $GENE_ID
        echo $GENE_ID >> ../Dataset/tranx_feature/processed_gene.txt
    fi
    echo ""
done
