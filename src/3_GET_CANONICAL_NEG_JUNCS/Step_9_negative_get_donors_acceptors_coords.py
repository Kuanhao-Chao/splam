from Bio import SeqIO
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

SAMPLE_NUM = 1261186
EACH_JUNC_PER_CHROM = 5000
MIN_JUNC = 400
hg38_ref = "../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa"
output_bed = "../NEG_junctions/neg_junctions.bed"
output_file = "../INPUTS/input_neg.fa"

def task(description, sequence):
    fw_donor = open("../NEG_junctions/donor/"+description+"_donor.bed", "w")
    fw_acceptor = open("../NEG_junctions/acceptor/"+description+"_acceptor.bed", "w")
    fw_da = open("../NEG_junctions/d_a/"+description+"_d_a.bed", "w")
    print("description: ", description)
    for i in range(EACH_JUNC_PER_CHROM):
        if i % 1000 == 0:
            print("\t", i)
        # print(description)
        # select_num = random.choice(range(0, chrs[description]-10000))
        select_num = random.randint(0, chrs[description]-10000)

        # Finding the 'GT'
        donor_idx = 0
        acceptor_idx = 0
        
        no_donor = False
        no_acceptor = False
        while not(sequence[select_num+donor_idx] == "G" and sequence[select_num+donor_idx+1] == "T"):
            donor_idx += 1

            if select_num+donor_idx+1 >= chrs[description]:
                no_donor = True
                i -= 1
                break
        if no_donor:
            continue

        while not(sequence[select_num+donor_idx+acceptor_idx-2] == "A" and sequence[select_num+donor_idx+acceptor_idx-1] == "G" and acceptor_idx > MIN_JUNC):
            acceptor_idx += 1

            if select_num+donor_idx+acceptor_idx >= chrs[description]:
                no_acceptor = True
                i -= 1
                break
        if no_acceptor:
            continue

        # print("Donor    (", select_num+donor_idx, "): ", sequence[select_num+donor_idx: select_num+donor_idx+2])
        # print("Acceptor (", select_num+donor_idx+acceptor_idx,"):", sequence[select_num+donor_idx+acceptor_idx-2: select_num+donor_idx+acceptor_idx])

        donor_s = select_num+donor_idx-200
        donor_e = select_num+donor_idx+200
        acceptor_s = select_num+donor_idx+acceptor_idx-200
        acceptor_e = select_num+donor_idx+acceptor_idx+200
        fw_da.write(description + "\t" + str(select_num+donor_idx) + "\t" + str(select_num+donor_idx+acceptor_idx+1) + "\t" + "JUNC_donor\t1\t+\n")

        if donor_e > chrs[description] or acceptor_e > chrs[description]:
            i -= 1
            continue
        fw_donor.write(description + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + "JUNC_donor\t1\t+\n")
        fw_acceptor.write(description + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + "JUNC_donor\t1\t+\n")


def main():
    with open(hg38_ref, 'r') as handle:
        # print("Len: ", len(SeqIO.parse(handle, 'fasta')))
        os.mkdir("NEG_junctions/donor/")
        os.mkdir("NEG_junctions/acceptor/")
        os.mkdir("NEG_junctions/d_a/")

        workers = 20
        with ThreadPoolExecutor(workers) as pool:
            for record in SeqIO.parse(handle, 'fasta'):            
                # Extract individual parts of the FASTA record
                identifier = record.id
                description = record.description
                sequence = record.seq
                sequence = str(sequence).upper()
                if (description in chrs.keys()):
                    processed = pool.submit(task, description, sequence)

if __name__ == "__main__":
    main()