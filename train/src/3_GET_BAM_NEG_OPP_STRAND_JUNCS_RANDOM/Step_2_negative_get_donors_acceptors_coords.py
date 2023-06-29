from Bio import SeqIO
import random
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

SEQ_LENGTH="800"
QUATER_SEQ_LEN = int(SEQ_LENGTH) // 4
EACH_JUNC_PER_LOCUS = 500
MIN_JUNC = 200
MAX_JUNC = 20000

THRESHOLD = "100"
hg38_ref = "../../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa"
output_bed = "./NEG_rev_junctions/"+SEQ_LENGTH+"bp/neg_junctions.bed"
output_file = "../INPUTS/"+SEQ_LENGTH+"bp/input_can_neg.fa"

os.makedirs("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/", exist_ok=True)
os.makedirs("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/", exist_ok=True)
os.makedirs("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/d_a/", exist_ok=True)
os.makedirs("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/donor/", exist_ok=True)
os.makedirs("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/acceptor/", exist_ok=True)


fw_donor = open("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/donor/donor.bed", "w")
fw_acceptor = open("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/acceptor/acceptor.bed", "w")
fw_da = open("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/d_a/d_a.bed", "w")

def task(chromosome, sequence, start, end ,strand):
    idx = 0
    JUNC_BASE = 0
    while True:
        idx += 1
        if idx > EACH_JUNC_PER_LOCUS:
            break
        base_num = random.randint(start, end)
        JUNC_BASE = random.randint(MIN_JUNC, MAX_JUNC)
        
        ################################
        # Finding the first and second dimers
        ################################
        donor_idx = 0
        acceptor_idx = JUNC_BASE
        
        has_donor = False
        has_acceptor = False

        donor_1 = ""
        donor_2 = ""
        acceptor_1 = ""
        acceptor_1 = ""

        strand_select = "."
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

        
        ################################
        # Finding the first dimer
        ################################
        while True:
            # (1) condition to break => Found the first dimer
            if (sequence[base_num+donor_idx] == donor_1 and sequence[base_num+donor_idx+1] == donor_2):
                has_donor = True
                break

            # (2) condition to break => donor is out of range.
            if base_num+donor_idx+1 >= end:
                has_donor = True
                print("\tBreak because of out of range.")
                break

            # (3) condition to break => Trying too many times.
            # if donor_idx > 10000:
            #     has_donor = True
            #     print("\tBreak because of try too many times ~ donor_idx")
            #     break
            donor_idx += 1
        if not has_donor:
            continue
        
        ################################
        # Finding the second dimer
        ################################
        while True:
            # (1) condition to break => Found the second dimer
            if (sequence[base_num+donor_idx+acceptor_idx-2] == acceptor_1 and sequence[base_num+donor_idx+acceptor_idx-1] == acceptor_2):
                has_acceptor = True
                break

            # (2) condition to break => acceptor is out of range.
            if base_num+donor_idx+acceptor_idx > end:
                has_acceptor = False
                break
            acceptor_idx += 1

        if not has_acceptor:
            continue

        if has_donor and has_acceptor:



            # print("donor_1   : ", sequence[base_num+donor_idx])
            # print("donor_2   : ", sequence[base_num+donor_idx+1])
            # print("acceptor_1: ", sequence[base_num+donor_idx+acceptor_idx-2])
            # print("acceptor_2: ", sequence[base_num+donor_idx+acceptor_idx-1])


            splice_junc_len = acceptor_idx + 1


            flanking_size = QUATER_SEQ_LEN
            if splice_junc_len < QUATER_SEQ_LEN:
                flanking_size = splice_junc_len

            donor = base_num+donor_idx
            acceptor = base_num+donor_idx+acceptor_idx
            if (strand_select == "+"):
                donor_s = donor - QUATER_SEQ_LEN
                donor_e = donor + flanking_size
                acceptor_s = acceptor - flanking_size
                acceptor_e = acceptor + QUATER_SEQ_LEN

            elif (strand_select == "-"):
                donor_s = donor - flanking_size
                donor_e = donor + QUATER_SEQ_LEN
                acceptor_s = acceptor - QUATER_SEQ_LEN
                acceptor_e = acceptor + flanking_size

                if (flanking_size < 200 ):
                    print(donor_s, donor, donor_e, acceptor_s, acceptor_e)

            ######################################################
            # Check if the donors and acceptors are in range.
            ######################################################
            if acceptor_s >= end:
                continue
            
            fw_da.write(chromosome + "\t" + str(base_num+donor_idx) + "\t" + str(base_num+donor_idx+acceptor_idx+1) + "\t" + "JUNC_donor\t1\t"+strand_select+"\n")
            

            if (strand_select == "+"):
                fw_donor.write(chromosome + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + "JUNC_donor\t1\t"+strand_select+"\n")
                fw_acceptor.write(chromosome + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + "JUNC_donor\t1\t"+strand_select+"\n")

            elif (strand_select == "-"):
                fw_acceptor.write(chromosome + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + "JUNC_donor\t1\t"+strand_select+"\n")
                fw_donor.write(chromosome + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + "JUNC_donor\t1\t"+strand_select+"\n")



def main():

    with open(hg38_ref, 'r') as handle:
        Path("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/donor/").mkdir(parents=True, exist_ok=True)
        Path("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/acceptor/").mkdir(parents=True, exist_ok=True)
        Path("./NEG_rev_junctions/"+SEQ_LENGTH+"bp/d_a/").mkdir(parents=True, exist_ok=True)

        record_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))

        with open("./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.merge.sort.bed") as f:
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