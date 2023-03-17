from Bio import SeqIO
import random
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

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

chrs = get_hg38_chrom_size()

SEQ_LENGTH="800"
QUATER_SEQ_LEN = int(SEQ_LENGTH) // 4
EACH_JUNC_PER_LOCUS = 1000
MIN_JUNC = 300
THRESHOLD = "100"
hg38_ref = "../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa"
os.makedirs("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/", exist_ok=True)
output_bed = "./NEG_rev_junctions/"+SEQ_LENGTH+"bp/neg_junctions.bed"
output_file = "../INPUTS/"+SEQ_LENGTH+"bp/input_can_neg.fa"
D_A_POSITIONS = set()
os.makedirs("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/", exist_ok=True)


fw_donor = open("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/donor/donor.bed", "w")
fw_acceptor = open("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/acceptor/acceptor.bed", "w")
fw_da = open("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/d_a/d_a.bed", "w")

def task(chromosome, sequence, start, end ,strand):
    idx = 0
    while True:
        idx += 1
        # print("Current idx: ", idx)
        if idx > EACH_JUNC_PER_LOCUS:
            break
        # if idx % 1000 == 0:
        #     print("\t", idx)
        base_num = random.randint(start, end)
        acceptor_skip_idx = random.randint(0, 3)
        # Finding the 'GT'
        donor_idx = 0
        acceptor_idx = MIN_JUNC
        
        no_donor = False
        no_acceptor = False

        donor_1 = ""
        donor_2 = ""
        acceptor_1 = ""
        acceptor_1 = ""

        strand_target = "."
        if strand == "+":
            strand_select = "-"
        elif strand == "-":
            strand_select = "+"

        if strand_select == "-":
            # I need to find junctions on the reverse strand (CT-AC pairs).
            donor_1 = "C"
            donor_2 = "T"
            acceptor_1 = "A"
            acceptor_2 = "C"

        elif strand_select == "+":
            # I need to find junctions on the forward strand (GT-AG pairs).
            donor_1 = "G"
            donor_2 = "T"
            acceptor_1 = "A"
            acceptor_2 = "G"

        
        while not((base_num+donor_idx+1) < end and sequence[base_num+donor_idx] == donor_1 and sequence[base_num+donor_idx+1] == donor_2):
            # if chromosome == "chr13":
            #     print("donor_idx: ", donor_idx)
            donor_idx += 1
            if base_num+donor_idx+1 >= end:
                no_donor = True
                print("\tBreak because of out of range.")
                break
            if donor_idx > 10000:
                no_donor = True
                print("\tBreak because of try too many times ~ donor_idx")
                break
        if no_donor:
            continue
        
        # print("Find donor dimers!")


        while not(base_num+donor_idx+acceptor_idx <= end and sequence[base_num+donor_idx+acceptor_idx-2] == acceptor_1 and sequence[base_num+donor_idx+acceptor_idx-1] == acceptor_2 and acceptor_skip_idx > acceptor_skip_idx):
            acceptor_idx += 1
            if base_num+donor_idx+acceptor_idx > end:
                no_acceptor = True
                break
            if acceptor_idx > 100000:
                no_donor = True
                break
        if no_acceptor:
            continue

        donor_s = base_num+donor_idx-QUATER_SEQ_LEN
        donor_e = base_num+donor_idx+QUATER_SEQ_LEN
        acceptor_s = base_num+donor_idx+acceptor_idx-QUATER_SEQ_LEN
        acceptor_e = base_num+donor_idx+acceptor_idx+QUATER_SEQ_LEN


        ######################################################
        # Check if the donors and acceptors are in range.
        ######################################################
        if acceptor_s >= end:
            continue
        
        print("Donor: ", donor_s, "-", donor_e, " ;  Acceptor: ", acceptor_s, "-", acceptor_e)
        ######################################################
        # Make sure there are no donors or acceptors inside the sequences.
        ######################################################
        # in_da_set = False
        # for d in range(donor_s, donor_e):
        #     if (chromosome, d) in D_A_POSITIONS:
        #         in_da_set = True
        # for a in range(acceptor_s, acceptor_e):
        #     if (chromosome, a) in D_A_POSITIONS:
        #         in_da_set = True
        # if in_da_set:
        #     continue

        # D_A_POSITIONS.add((chromosome, base_num+donor_idx))
        # D_A_POSITIONS.add((chromosome, base_num+donor_idx+acceptor_idx))

        fw_da.write(chromosome + "\t" + str(base_num+donor_idx) + "\t" + str(base_num+donor_idx+acceptor_idx+1) + "\t" + "JUNC_donor\t1\t"+strand_select+"\n")
        fw_donor.write(chromosome + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + "JUNC_donor\t1\t+\n")
        fw_acceptor.write(chromosome + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + "JUNC_donor\t1\t"+strand_select+"\n")
    
    
    
    
    


def main():
    # with open("../1_GET_BAM_JUNCS/BAM_junctions/"+SEQ_LENGTH+"bp/"+str(THRESHOLD)+"_juncs/d_a.bed", "r") as f:
    #     lines = f.read().splitlines()
    #     for line in lines:
    #         # print(line)
    #         eles = line.split("\t")
    #         # print("eles[0], eles[1]: ", eles[0], eles[1], eles[2])
    #         # Adding donor into the set.
    #         D_A_POSITIONS.add((eles[0], eles[1]))
    #         # Adding acceptor into the set.
    #         D_A_POSITIONS.add((eles[0], eles[2]))

    # with open("../2_GET_REF_JUNCS/REF_junctions/ref_d_a.sort.bed", "r") as f:
    #     lines = f.read().splitlines()
    #     for line in lines:
    #         # print(line)
    #         eles = line.split("\t")
    #         # print("eles[0], eles[1]: ", eles[0], eles[1], eles[2])
    #         # Adding donor into the set.
    #         D_A_POSITIONS.add((eles[0], eles[1]))
    #         # Adding acceptor into the set.
    #         D_A_POSITIONS.add((eles[0], eles[2]))

    # print("D_A_POSITIONS size: ", len(D_A_POSITIONS))

    with open(hg38_ref, 'r') as handle:
        Path("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/donor/").mkdir(parents=True, exist_ok=True)
        Path("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/acceptor/").mkdir(parents=True, exist_ok=True)
        Path("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/d_a/").mkdir(parents=True, exist_ok=True)

        record_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))

        with open("../6_GET_BAM_NEG_OPP_STRAND_JUNCS/BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.pos.only.bed") as f:
            lines = f.read().splitlines()
            for line in lines:
                eles = line.split("\t")
                print(eles)

                chr = eles[0]
                start = int(eles[1])
                end = int(eles[2])
                strand = eles[3]


                record = record_dict[chr]

                # Extract individual parts of the FASTA record
                identifier = record.id
                chromosome = record.description
                sequence = record.seq
                sequence = str(sequence).upper()
                # if (chromosome in targets.keys()):
                task(chromosome, sequence, start, end ,strand)

    fw_donor.close()
    fw_acceptor.close()
    fw_da.close()

if __name__ == "__main__":
    main()