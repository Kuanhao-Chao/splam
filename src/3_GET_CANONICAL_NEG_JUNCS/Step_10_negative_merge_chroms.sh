mkdir ./NEG_junctions/pre_filter/
for DONOR_F in "NEG_junctions/donor"/*
do
    echo "$DONOR_F"
    sort -k2,2n -k3,3n "$DONOR_F" | uniq -u >> ./NEG_junctions/pre_filter/donor.bed
done

for ACCEPTOR_F in "NEG_junctions/acceptor"/*
do
    echo "$ACCEPTOR_F"
    sort -k2,2n -k3,3n "$ACCEPTOR_F" | uniq -u  >> ./NEG_junctions/pre_filter/acceptor.bed
done

for D_A in "NEG_junctions/d_a"/*
do
    echo "$D_A"
    sort -k2,2n -k3,3n "$D_A" | uniq -u  >> ./NEG_junctions/pre_filter/d_a.bed
done