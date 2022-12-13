import pandas as pd
import random
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

targets = {
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

def get_hg38_chrom_size():
    f_chrs = open("../hg38_chrom_size.tsv", "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = int(eles[1])
    return chrs


THRESHOLD = "100"
EACH_CHR = 10000
D_A_POSITIONS = set()
MIN_GAP = 400
MAX_GAP = 2000
SEQ_LEN="600"
QUATER_SEQ_LEN = int(SEQ_LEN) // 4
chrs = get_hg38_chrom_size()

#################################
# 1. Confirming that the junction does not overlap with any BAM junctions.
# 2. Confirming that the junction doesn not overlap with any REF junctions.
#################################

#################################
# For 'd_a.bed': 0-based, 1-based
# For 'donor.bed': 0-based, 0-based
# For 'acceptor.bed': 0-based, 0-based
#################################
def task(chromosome):
    fw_donor = open("../NEG_noncan_junctions/"+SEQ_LEN+"bp/donor/"+chromosome+"_donor.bed", "w")
    fw_acceptor = open("../NEG_noncan_junctions/"+SEQ_LEN+"bp/acceptor/"+chromosome+"_acceptor.bed", "w")
    fw_da = open("../NEG_noncan_junctions/"+SEQ_LEN+"bp/d_a/"+chromosome+"_d_a.bed", "w")
    print("chromosome: ", chromosome)
    # for idx in range(EACH_CHR):
    idx = 0    
    while True:
        if idx >= EACH_CHR:
            break
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

            for d in range(select_donor-QUATER_SEQ_LEN, select_donor+QUATER_SEQ_LEN):
                if (chromosome, d) in D_A_POSITIONS:
                    in_da_set = True
            for a in range(select_acceptor-QUATER_SEQ_LEN, select_acceptor+QUATER_SEQ_LEN):
                if (chromosome, a) in D_A_POSITIONS:
                    in_da_set = True
            if not in_da_set:
                break
        donor_s = select_donor-QUATER_SEQ_LEN
        donor_e = select_donor+QUATER_SEQ_LEN
        acceptor_s = select_acceptor-QUATER_SEQ_LEN
        acceptor_e = select_acceptor+QUATER_SEQ_LEN
        if donor_e > chrs[chromosome] or acceptor_e > chrs[chromosome] or donor_s <= 0 or acceptor_s <= 0:
            continue

        D_A_POSITIONS.add((chromosome, select_donor))
        D_A_POSITIONS.add((chromosome, select_acceptor))

        fw_da.write(chromosome + "\t" + str(select_donor) + "\t" + str(select_acceptor+1) + "\t" + "JUNC\t1\t+\n")
        fw_donor.write(chromosome + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + "JUNC_donor\t1\t+\n")
        fw_acceptor.write(chromosome + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + "JUNC_acceptor\t1\t+\n")
        idx += 1

def main():
    Path("../NEG_noncan_junctions/"+SEQ_LEN+"bp/donor/").mkdir(parents=True, exist_ok=True)
    Path("../NEG_noncan_junctions/"+SEQ_LEN+"bp/acceptor/").mkdir(parents=True, exist_ok=True)
    Path("../NEG_noncan_junctions/"+SEQ_LEN+"bp/d_a/").mkdir(parents=True, exist_ok=True)

    with open("../BAM_junctions/"+SEQ_LEN+"bp/"+str(THRESHOLD)+"_juncs/d_a.bed", "r") as f:
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

    for chromosome in targets.keys():
        task(chromosome)

if __name__ == "__main__":
    main()