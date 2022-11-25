# mkdir ../NEG_noncan_junctions/pre_filter/
for DONOR_F in "../NEG_noncan_junctions/donor"/*
do
    echo "$DONOR_F"
    sort -k2,2n -k3,3n "$DONOR_F" >> ../NEG_noncan_junctions/donor.bed
done

for ACCEPTOR_F in "../NEG_noncan_junctions/acceptor"/*
do
    echo "$ACCEPTOR_F"
    sort -k2,2n -k3,3n "$ACCEPTOR_F" >> ../NEG_noncan_junctions/acceptor.bed
done

for D_A in "../NEG_noncan_junctions/d_a"/*
do
    echo "$D_A"
    sort -k2,2n -k3,3n "$D_A"  >> ../NEG_noncan_junctions/d_a.bed
done