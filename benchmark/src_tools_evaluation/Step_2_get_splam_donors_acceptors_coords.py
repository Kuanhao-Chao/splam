import pandas as pd
import os 

def get_hg38_chrom_size():
    f_chrs = open("../src/hg38_chrom_size.tsv", "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = int(eles[1])
    return chrs


def main():
    chrs = get_hg38_chrom_size()
    SEQ_LEN="800"
    QUATER_SEQ_LEN = int(SEQ_LEN) // 4

    #################################
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    #################################
    output_dir = "./dataset/"
    output_files = [output_dir+"pos/", output_dir+"pos_MANE/", output_dir+"pos_ALTS/", output_dir+"neg_1/", output_dir+"neg_random/"]

    splam_dir = "splam/"
    for output_file in output_files:
        fw_donor = open(output_file + splam_dir + "splam.juncs.donor.bed", "w")
        fw_acceptor = open(output_file + splam_dir + "splam.juncs.acceptor.bed", "w")
        JUNCS = set()

        with open(output_file + splam_dir + "splam.juncs.bed", "r") as f:
            lines = f.read().splitlines()
            for line in lines:
                eles = line.split("\t")
                chr = eles[0]
                junc_name = eles[3]
                score = eles[4]
                strand = eles[5]
                label = eles[6]
                donor = 0
                acceptor = 0
                if (strand == "+"):
                    donor = int(eles[1])
                    acceptor = int(eles[2])
                    # -1
                    splice_junc_len = acceptor - donor
                elif (strand == "-"):
                    donor = int(eles[2])
                    # -1
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
                    continue
                if donor_s < 0 or acceptor_s < 0:
                    continue
                new_junc = (chr, str(donor_s), str(donor_e), str(acceptor_s), str(acceptor_e), strand)
                if new_junc in JUNCS:
                    continue
                else:
                    JUNCS.add(new_junc)
                    fw_donor.write(chr + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + junc_name+"_donor" + "\t" + score + "\t" + strand + "\t" + label + "\n")
                    fw_acceptor.write(chr + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + score + "\t" + strand + "\t" + label + "\n")
        fw_donor.close()
        fw_acceptor.close()

if __name__ == "__main__":
    main()