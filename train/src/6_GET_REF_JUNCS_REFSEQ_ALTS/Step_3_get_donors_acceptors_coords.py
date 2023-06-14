import pandas as pd
import os

def get_hg38_chrom_size():
    f_chrs = open("../hg38_chrom_size.tsv", "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = int(eles[1])
    return chrs

def main():
    THRESHOLD = "100"
    SEQ_LEN="800"
    QUATER_SEQ_LEN = int(SEQ_LEN) // 4

    ref_juncs = pd.read_csv('./REF_junctions/ref_d_a.sort.bed', sep="\t", header=None)
    
    print(ref_juncs)
    
    # Calling merge() function
    os.makedirs('./BAM_REF_Intersection/'+SEQ_LEN+"bp/"+THRESHOLD+'_juncs/', exist_ok=True)
    d_a_out = './BAM_REF_Intersection/'+SEQ_LEN+"bp/"+THRESHOLD+'_juncs/d_a.bed'
    d_out = './BAM_REF_Intersection/'+SEQ_LEN+"bp/"+THRESHOLD+'_juncs/donor.bed'
    a_out = './BAM_REF_Intersection/'+SEQ_LEN+"bp/"+THRESHOLD+'_juncs/acceptor.bed'

    # intersect_df = pd.merge(ref_juncs, bam_juncs, how ='inner', on =[0, 1, 2, 5])
    # print("intersect_df: ", intersect_df)
    # out_df = intersect_df[[0, 1, 2, "3_x", "4_x", 5]]
    ref_juncs = ref_juncs.rename(columns={0:"chr",1:"start", 2:"end", 3:"junc", 4:"score", 5:"strand"})

    ref_juncs.to_csv(d_a_out, sep="\t", index=False, header=None)
    chrs = get_hg38_chrom_size()
    fw_donor = open(d_out, "w")
    fw_acceptor = open(a_out, "w")

    JUNCS = set()

    skip_counter = 0
    with open(d_a_out, "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            eles = line.split("\t")
            chr = eles[0]
            junc_name = eles[3]
            score = eles[4]
            strand = eles[5]
            donor = 0
            acceptor = 0
            if (strand == "+"):
                donor = int(eles[1])
                acceptor = int(eles[2])-1
                splice_junc_len = acceptor - donor
            elif (strand == "-"):
                donor = int(eles[2])-1
                acceptor = int(eles[1])
                splice_junc_len = donor - acceptor
    
            flanking_size = QUATER_SEQ_LEN
            if splice_junc_len < QUATER_SEQ_LEN:
                flanking_size = splice_junc_len

            if (strand == "+"):
                donor_s = donor - QUATER_SEQ_LEN
                donor_e = donor + flanking_size
                acceptor_s = acceptor - flanking_size
                acceptor_e = acceptor + QUATER_SEQ_LEN

            elif (strand == "-"):
                donor_s = donor - flanking_size
                donor_e = donor + QUATER_SEQ_LEN
                acceptor_s = acceptor - QUATER_SEQ_LEN
                acceptor_e = acceptor + flanking_size

            if donor_e >= chrs[chr] or acceptor_e >= chrs[chr]:
                # skip_counter += 1
                continue
            if donor_s < 0 or acceptor_s < 0:
                # skip_counter += 1
                continue
            new_junc = (chr, str(donor_s), str(donor_e), str(acceptor_s), str(acceptor_e), strand)
            if new_junc in JUNCS:
                skip_counter += 1
                continue
            else:
                JUNCS.add(new_junc)
                fw_donor.write(chr + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + junc_name+"_donor" + "\t" + score + "\t" + strand + "\n")

                fw_acceptor.write(chr + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + score + "\t" + strand + "\n")

    fw_donor.close()
    fw_acceptor.close()

if __name__ == "__main__":
    main()