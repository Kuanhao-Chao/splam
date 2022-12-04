import pandas as pd
import os 
import sys

def get_hg38_chrom_size():
    f_chrs = open("../hg38_chrom_size.tsv", "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = int(eles[1])
    return chrs


def main(argv):
    chrs = get_hg38_chrom_size()

    threshold = "100"

    #################################
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    #################################
    os.makedirs("../../results/"+argv[0]+"/juncs/", exist_ok=True)
    fw_donor = open("../../results/"+argv[0]+"/juncs/donor.bed", "w")
    fw_acceptor = open("../../results/"+argv[0]+"/juncs/acceptor.bed", "w")
    
    d_a_bed = "../../results/"+argv[0]+"/juncs/d_a.bed"
    fw_da = open(d_a_bed, "w")
    # fw_d = open("BAM_junctions/d.bed", "w")
    # fw_a = open("BAM_junctions/a.bed", "w")
    JUNCS = set()

    with open("../../Dataset/"+argv[0]+"/"+argv[0]+".bed", "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            eles = line.split("\t")

            chr = eles[0]
            junc_name = eles[3]
            score = eles[4]
            strand = eles[5]

            lengths = eles[10].split(',')
            len_1 = int(lengths[0])
            len_2 = int(lengths[1])
            if (strand == "+"):
                donor = int(eles[1]) + len_1
                acceptor = int(eles[2]) - len_2
                splice_junc_len = acceptor - donor
            elif (strand == "-"):
                acceptor = int(eles[1]) + len_1
                donor = int(eles[2]) - len_2
                splice_junc_len = donor - acceptor

            flanking_size = 200
            if splice_junc_len < 400:
                flanking_size = splice_junc_len // 2

            if (strand == "+"):
                donor_s = donor - 200
                donor_e = donor + flanking_size
                acceptor_s = acceptor - flanking_size
                acceptor_e = acceptor + 200

            elif (strand == "-"):
                donor_s = donor - flanking_size
                donor_e = donor + 200
                acceptor_s = acceptor - 200
                acceptor_e = acceptor + flanking_size
                
            # if chr == "chr22_KI270733v1_random" or chr == "chr22_KI270734v1_random":
            #     continue

            if donor_e >= chrs[chr] or acceptor_e >= chrs[chr]:
                continue
            if donor_s < 0 or acceptor_s < 0:
                continue
            new_junc = (chr, str(donor_s), str(donor_e), str(acceptor_s), str(acceptor_e), strand)
            if new_junc in JUNCS:
                continue
            else:
                JUNCS.add(new_junc)
                fw_donor.write(chr + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + junc_name+"_donor" + "\t" + score + "\t" + strand + "\n")
                fw_acceptor.write(chr + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + score + "\t" + strand + "\n")

                if (strand == "+"):
                    fw_da.write(chr + "\t" + str(donor) + "\t" + str(acceptor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")
                elif (strand == "-"):
                    fw_da.write(chr + "\t" + str(acceptor) + "\t" + str(donor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")

    fw_donor.close()
    fw_acceptor.close()
    fw_da.close()

if __name__ == "__main__":
    main(sys.argv[1:])