import pandas as pd
import random
import os
from concurrent.futures import ThreadPoolExecutor

chrs = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415
}
threshold = 5
EACH_CHR = 10000
D_A_POSITIONS = set()
MIN_GAP = 400
MAX_GAP = 2000

def task(chromosome):
    fw_donor = open("../NEG_noncan_junctions/donor/"+chromosome+"_donor.bed", "w")
    fw_acceptor = open("../NEG_noncan_junctions/acceptor/"+chromosome+"_acceptor.bed", "w")
    fw_da = open("../NEG_noncan_junctions/d_a/"+chromosome+"_d_a.bed", "w")
    print("chromosome: ", chromosome)
    for idx in range(EACH_CHR):
        
        in_da_set = False
        select_donor = 0
        select_acceptor = 0
        donor_s = 0
        donor_e = 0
        acceptor_s = 0
        acceptor_e = 0
        while True:
            in_da_set = False
            select_donor = random.randint(0, chrs[chromosome]-10000)
            select_acceptor = select_donor+random.randint(MIN_GAP, MAX_GAP)

            for d in range(select_donor-200, select_donor+200):
                # print("d: ", d)
                if (chromosome, d) in D_A_POSITIONS:
                    in_da_set = True
            for a in range(select_acceptor-200, select_acceptor+200):
                # print("a: ", a)
                if (chromosome, a) in D_A_POSITIONS:
                    in_da_set = True
            if not in_da_set:
                break
        donor_s = select_donor-200
        donor_e = select_donor+200
        acceptor_s = select_acceptor-200
        acceptor_e = select_acceptor+200

        fw_da.write(chromosome + "\t" + str(select_donor) + "\t" + str(select_acceptor+1) + "\t" + "JUNC\t1\t+\n")

        if donor_e > chrs[chromosome] or acceptor_e > chrs[chromosome]:
            i -= 1
            continue
        fw_donor.write(chromosome + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + "JUNC_donor\t1\t+\n")
        fw_acceptor.write(chromosome + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + "JUNC_acceptor\t1\t+\n")

def main():
    os.mkdir("../NEG_noncan_junctions/donor/")
    os.mkdir("../NEG_noncan_junctions/acceptor/")
    os.mkdir("../NEG_noncan_junctions/d_a/")
    with open("../BAM_junctions/junctions_"+str(threshold)+".bed", "r") as f:
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
            
            D_A_POSITIONS.add((chr, donor))
            D_A_POSITIONS.add((chr, acceptor))

            # if (strand == "+"):
            #     donor_s = donor - 200
            #     donor_e = donor + flanking_size
            #     acceptor_s = acceptor - flanking_size
            #     acceptor_e = acceptor + 200

            # elif (strand == "-"):
            #     donor_s = donor - flanking_size
            #     donor_e = donor + 200
            #     acceptor_s = acceptor - 200
            #     acceptor_e = acceptor + flanking_size
            # print("eles: ", eles)
    #         if chr == "chr22_KI270733v1_random" or chr == "chr22_KI270734v1_random":
    #             continue
    #         if donor_s > 0:
    #             new_junc = (chr, str(donor_s), str(donor_e), str(acceptor_s), str(acceptor_e), strand)
    #             if new_junc in JUNCS:
    #                 continue
    #             else:
    #                 JUNCS.add(new_junc)
    #                 fw_donor.write(chr + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + junc_name+"_donor" + "\t" + score + "\t" + strand + "\n")
    #                 fw_acceptor.write(chr + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + score + "\t" + strand + "\n")

    #                 if (strand == "+"):
    #                     fw_da.write(chr + "\t" + str(donor) + "\t" + str(acceptor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")
    #                 elif (strand == "-"):
    #                     fw_da.write(chr + "\t" + str(acceptor) + "\t" + str(donor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")

    # fw_donor.close()
    # fw_acceptor.close()
    # fw_da.close()

    for chromosome in chrs.keys():
        task(chromosome)


            # print(select_num)
        

if __name__ == "__main__":
    main()