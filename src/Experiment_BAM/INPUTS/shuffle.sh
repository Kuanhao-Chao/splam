awk '{OFS="\t"; getline seq; print $0,seq}' input_$1.fa | shuf | awk '{OFS="\n"; print $1,$2}' > input_$1.shuffle.fa
