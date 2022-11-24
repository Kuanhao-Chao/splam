import pandas as pd
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

def main():
    threshold = 100
    fw_donor = open("../BAM_junctions/donor.bed", "w")
    fw_acceptor = open("../BAM_junctions/acceptor.bed", "w")

    d_a_bed = "../BAM_junctions/d_a.bed"
    fw_da = open(d_a_bed, "w")
    # fw_d = open("BAM_junctions/d.bed", "w")
    # fw_a = open("BAM_junctions/a.bed", "w")
    JUNCS = set()

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
                
            if chr == "chr22_KI270733v1_random" or chr == "chr22_KI270734v1_random":
                continue
            if donor_s > 0:
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