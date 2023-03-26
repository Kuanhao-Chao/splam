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
    chrs = get_hg38_chrom_size()
    JUNC_POSITIONS = set()

    THRESHOLD = "1"
    SEQ_LEN = "800"
    HALF_SEQ_LEN = int(SEQ_LEN) // 2
    QUATER_SEQ_LEN = int(SEQ_LEN) // 4

    #################################
    # Adding all 100-alignment-supported splice sites into the 
    # 'D_A_POSITIONS' dictionary.
    #   I need to avoid these splice sites.
    #################################
    with open("../1_GET_BAM_JUNCS/BAM_junctions/"+SEQ_LEN+"bp/100_juncs/d_a.bed", "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            eles = line.split("\t")
            JUNC_POSITIONS.add((eles[0], eles[1], eles[2], eles[5]))


    #################################
    # Adding all reference-supported splice sites into the 
    # 'D_A_POSITIONS' dictionary.
    #   I need to avoid these splice sites.
    #################################
    with open("../2_GET_REF_JUNCS_REFSEQ/REF_junctions/ref_d_a.sort.bed", "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            eles = line.split("\t")
            JUNC_POSITIONS.add((eles[0], eles[1], eles[2], eles[5]))


    print("JUNC_POSITIONS: ", len(JUNC_POSITIONS))
    
    #################################
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    #################################
    os.makedirs("./BAM_junctions/"+SEQ_LEN+"bp/"+THRESHOLD+"_juncs/", exist_ok=True)
    fw_donor = open("./BAM_junctions/"+SEQ_LEN+"bp/"+THRESHOLD+"_juncs/donor.bed", "w")
    fw_acceptor = open("./BAM_junctions/"+SEQ_LEN+"bp/"+THRESHOLD+"_juncs/acceptor.bed", "w")
    
    d_a_bed = "./BAM_junctions/"+SEQ_LEN+"bp/"+THRESHOLD+"_juncs/d_a.bed"
    fw_da = open(d_a_bed, "w")
    JUNCS = set()

    
    print("./BAM_junctions/junctions_"+THRESHOLD+"_cleaned.bed")
    with open("./BAM_junctions/junctions_"+THRESHOLD+"_cleaned.bed", "r") as f:
        lines = f.read().splitlines()
        counter = 0
        for line in lines:
            eles = line.split("\t")

            chr = eles[0]
            junc_name = eles[3]
            score = eles[4]
            strand = eles[5]
            donor = 0
            acceptor = 0

            if (eles[0], eles[1], eles[2], eles[5]) in JUNC_POSITIONS:
                continue


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
            else:
                continue

            if (splice_junc_len < 0 and strand == "+"):
                print("splice_junc_len: ", splice_junc_len)

                print("lengths: ", lengths)
                print("eles[1] : ", eles[1])
                print("eles[2] : ", eles[2])
                print("donor   : ", donor)
                print("acceptor: ", acceptor)

            flanking_size = QUATER_SEQ_LEN
            if splice_junc_len < QUATER_SEQ_LEN:
                flanking_size = splice_junc_len
            if (flanking_size < 0):
                print("flanking_size: ", flanking_size)

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

                counter += 1
                fw_donor.write(chr + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + junc_name+"_donor" + "\t" + score + "\t" + strand + "\n")
                fw_acceptor.write(chr + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + score + "\t" + strand + "\n")

                if (strand == "+"):
                    fw_da.write(chr + "\t" + str(donor) + "\t" + str(acceptor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")
                elif (strand == "-"):
                    fw_da.write(chr + "\t" + str(acceptor) + "\t" + str(donor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")

    print("D_A_POSITIONS After: ", len(JUNC_POSITIONS))

    fw_donor.close()
    fw_acceptor.close()
    fw_da.close()

if __name__ == "__main__":
    main()