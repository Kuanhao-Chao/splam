# mkdir ../NEG_junctions/pre_filter/
SEQ_LENGTH=600
for DONOR_F in "../NEG_junctions/${SEQ_LENGTH}bp/donor"/*
do
    echo "$DONOR_F"
    sort -k2,2n -k3,3n "$DONOR_F" | uniq -u >> ../NEG_junctions/${SEQ_LENGTH}bp/donor.bed
done

for ACCEPTOR_F in "../NEG_junctions/${SEQ_LENGTH}bp/acceptor"/*
do
    echo "$ACCEPTOR_F"
    sort -k2,2n -k3,3n "$ACCEPTOR_F" | uniq -u  >> ../NEG_junctions/${SEQ_LENGTH}bp/acceptor.bed
done

for D_A in "../NEG_junctions/${SEQ_LENGTH}bp/d_a"/*
do
    echo "$D_A"
    sort -k2,2n -k3,3n "$D_A" | uniq -u  >> ../NEG_junctions/${SEQ_LENGTH}bp/d_a.bed
done