# mkdir ../TEST/NEG_junctions/pre_filter/
# TARGET="NEG_junctions"
TARGET="NEG_noncan_junctions"

for DONOR_F in "../TEST/${TARGET}/donor"/*
do
    echo "$DONOR_F"
    sort -k2,2n -k3,3n "$DONOR_F" | uniq -u >> ../TEST/$TARGET/donor.bed
done

for ACCEPTOR_F in "../TEST/${TARGET}/acceptor"/*
do
    echo "$ACCEPTOR_F"
    sort -k2,2n -k3,3n "$ACCEPTOR_F" | uniq -u  >> ../TEST/$TARGET/acceptor.bed
done

for D_A in "../TEST/${TARGET}/d_a"/*
do
    echo "$D_A"
    sort -k2,2n -k3,3n "$D_A" | uniq -u  >> ../TEST/$TARGET/d_a.bed
done