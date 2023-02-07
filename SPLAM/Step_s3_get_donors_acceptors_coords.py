import os
import pandas as pd

def get_hg38_chrom_size():
    f_chrs = open("../hg38_chrom_size.tsv", "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = int(eles[1])
    return chrs


def main():
    chrs = get_hg38_chrom_size()

    #################################
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    #################################
    os.makedirs("./juncs/", exist_ok=True)
    fw_donor = open("./juncs/donor.bed", "w")
    fw_acceptor = open("./juncs/acceptor.bed", "w")

    d_a_bed = "./juncs/d_a.bed"
    fw_da = open(d_a_bed, "w")
    JUNCS = set()
    with open("introns.sort.bed", "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            eles = line.split(" ")

            chr = eles[0]
            junc_name = eles[3]
            score = eles[4]
            strand = eles[5]

            if (strand == "+"):
                donor = int(eles[1])
                acceptor = int(eles[2])
                splice_junc_len = acceptor - donor
            elif (strand == "-"):
                acceptor = int(eles[1])
                donor = int(eles[2])
                splice_junc_len = donor - acceptor
            # print("donor   : ", donor)
            # print("acceptor: ", acceptor)
            # print("strand  : ", strand)
            # print("\n\n")

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
    main()