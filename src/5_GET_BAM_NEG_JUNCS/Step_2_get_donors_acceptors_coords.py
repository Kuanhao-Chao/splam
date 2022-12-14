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
    D_A_POSITIONS=set()

    THRESHOLD = "1"
    SEQ_LEN = "800"
    HALF_SEQ_LEN = int(SEQ_LEN) // 2
    QUATER_SEQ_LEN = int(SEQ_LEN) // 4

    #################################
    # Adding all splice sites into the 
    # 'D_A_POSITIONS' dictionary.
    #################################
    with open("../BAM_junctions/"+SEQ_LEN+"bp/"+str(100)+"_juncs/d_a.bed", "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            # print(line)
            eles = line.split("\t")
            # print("eles[0], eles[1]: ", eles[0], eles[1], eles[2])
            D_A_POSITIONS.add((eles[0], eles[1]))
            D_A_POSITIONS.add((eles[0], eles[2]))

    with open("../REF_junctions/ref_d_a.sort.bed", "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            # print(line)
            eles = line.split("\t")
            # print("eles[0], eles[1]: ", eles[0], eles[1], eles[2])
            D_A_POSITIONS.add((eles[0], eles[1]))
            D_A_POSITIONS.add((eles[0], eles[2]))

    print("D_A_POSITIONS: ", len(D_A_POSITIONS))
    
    #################################
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    #################################
    os.makedirs("../BAM_junctions/"+SEQ_LEN+"bp/"+THRESHOLD+"_juncs/", exist_ok=True)
    fw_donor = open("../BAM_junctions/"+SEQ_LEN+"bp/"+THRESHOLD+"_juncs/donor.bed", "w")
    fw_acceptor = open("../BAM_junctions/"+SEQ_LEN+"bp/"+THRESHOLD+"_juncs/acceptor.bed", "w")
    
    d_a_bed = "../BAM_junctions/"+SEQ_LEN+"bp/"+THRESHOLD+"_juncs/d_a.bed"
    fw_da = open(d_a_bed, "w")
    JUNCS = set()

    with open("../BAM_junctions/junctions_"+str(THRESHOLD)+".bed", "r") as f:
        # lines = f.read().splitlines()
        # for line in lines:
        #     eles = line.split("\t")

        #     chr = eles[0]
        #     junc_name = eles[3]
        #     score = eles[4]
        #     strand = eles[5]

        #     lengths = eles[10].split(',')
        #     len_1 = int(lengths[0])
        #     len_2 = int(lengths[1])
        #     if (strand == "+"):
        #         donor = int(eles[1]) + len_1
        #         acceptor = int(eles[2]) - len_2
        #         splice_junc_len = acceptor - donor
        #     elif (strand == "-"):
        #         acceptor = int(eles[1]) + len_1
        #         donor = int(eles[2]) - len_2
        #         splice_junc_len = donor - acceptor

        #     flanking_size = 250
        #     if splice_junc_len < 250:
        #         flanking_size = splice_junc_len
        #         # flanking_size = splice_junc_len // 2

        #     if (strand == "+"):
        #         donor_s = donor - 250
        #         donor_e = donor + flanking_size
        #         acceptor_s = acceptor - flanking_size
        #         acceptor_e = acceptor + 250

        #     elif (strand == "-"):
        #         donor_s = donor - flanking_size
        #         donor_e = donor + 250
        #         acceptor_s = acceptor - 250
        #         acceptor_e = acceptor + flanking_size
                

        #     if donor_e >= chrs[chr] or acceptor_e >= chrs[chr]:
        #         continue
        #     if donor_s < 0 or acceptor_s < 0:
        #         continue
        #     new_junc = (chr, str(donor_s), str(donor_e), str(acceptor_s), str(acceptor_e), strand)
        #     if new_junc in JUNCS:
        #         continue
        #     else:
        #         JUNCS.add(new_junc)
        #         fw_donor.write(chr + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + junc_name+"_donor" + "\t" + score + "\t" + strand + "\n")
        #         fw_acceptor.write(chr + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + score + "\t" + strand + "\n")

        #         if (strand == "+"):
        #             fw_da.write(chr + "\t" + str(donor) + "\t" + str(acceptor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")
        #         elif (strand == "-"):
        #             fw_da.write(chr + "\t" + str(acceptor) + "\t" + str(donor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")










        lines = f.read().splitlines()
        counter = 0
        for line in lines:
            if counter >= 100000: break
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

            flanking_size = QUATER_SEQ_LEN
            if splice_junc_len < QUATER_SEQ_LEN:
                flanking_size = splice_junc_len
                # flanking_size = splice_junc_len // 2

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
                #################################
                # Check if the junction is in the 'D_A_POSITIONS' 
                #################################
                in_da_set = False
                for d in range(donor_s, donor_e):
                    if in_da_set: break
                    if (chr, d) in D_A_POSITIONS:
                        in_da_set = True
                for a in range(acceptor_s, acceptor_e):
                    if in_da_set: break
                    if (chr, a) in D_A_POSITIONS:
                        in_da_set = True
                # print("in_da_set: ", in_da_set)
                D_A_POSITIONS.add((chr, donor_s))
                D_A_POSITIONS.add((chr, donor_e))
                D_A_POSITIONS.add((chr, acceptor_s))
                D_A_POSITIONS.add((chr, acceptor_e))
                if not in_da_set:
                    counter += 1
                    fw_donor.write(chr + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + junc_name+"_donor" + "\t" + score + "\t" + strand + "\n")
                    fw_acceptor.write(chr + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + score + "\t" + strand + "\n")

                    if (strand == "+"):
                        fw_da.write(chr + "\t" + str(donor) + "\t" + str(acceptor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")
                    elif (strand == "-"):
                        fw_da.write(chr + "\t" + str(acceptor) + "\t" + str(donor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")

    print("D_A_POSITIONS After: ", len(D_A_POSITIONS))

    fw_donor.close()
    fw_acceptor.close()
    fw_da.close()

if __name__ == "__main__":
    main()