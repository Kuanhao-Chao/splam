import pandas as pd
import os 
import sys

def get_hg38_chrom_size():
    
    chrs = {}
    with open('GRCh38.p14_assembly_report.txt', 'r') as file:       
        # skip header
        next(file)

        # read the file line by line
        for line in file:  
            # split by tabs
            columns = line.strip().split('\t')
            ucsc_name = columns[9]
            chrom_size = int(columns[8])

            # store the key-value pair in the dictionary
            chrs[ucsc_name] = chrom_size
    
    return chrs

def get_chrom_size(path):
    chrs = {}
    with open(path, 'r') as file:       
        # skip header
        next(file)

        # read the file line by line
        for line in file:  
            # split by tabs
            columns = line.strip().split('\t')
            refseq_name = columns[6]
            chrom_size = int(columns[8])

            # store the key-value pair in the dictionary
            chrs[refseq_name] = chrom_size
    
    return chrs


def main(argv):
    #chrs = get_hg38_chrom_size()
    
    chrs = get_chrom_size('../Dataset/' + argv[1].split('_parsed.bed')[0] + '_assembly_report.txt')

    threshold = "100"
    SEQ_LEN="800"
    QUARTER_SEQ_LEN = int(SEQ_LEN) // 4

    #################################
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    #################################
    os.makedirs(argv[0]+"/juncs/", exist_ok=True)
    fw_donor = open(argv[0]+"/juncs/donor.bed", "w")
    fw_acceptor = open(argv[0]+"/juncs/acceptor.bed", "w")
    
    d_a_bed = argv[0]+"/juncs/d_a.bed"
    fw_da = open(d_a_bed, "w")
    JUNCS = set()

    # with open(argv[0]+"/junction.bed", "r") as f:
    # with open("../Dataset/"+argv[0]+"/"+argv[0]+".bed", "r") as f:
    with open(argv[0]+'/'+argv[1], 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            eles = line.split("\t")
            if len(eles) == 1:
                continue
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

            flanking_size = QUARTER_SEQ_LEN
            if splice_junc_len < QUARTER_SEQ_LEN:
                flanking_size = splice_junc_len

            if (strand == "+"):
                donor_s = donor - QUARTER_SEQ_LEN
                donor_e = donor + flanking_size
                acceptor_s = acceptor - flanking_size
                acceptor_e = acceptor + QUARTER_SEQ_LEN

            elif (strand == "-"):
                donor_s = donor - flanking_size
                donor_e = donor + QUARTER_SEQ_LEN
                acceptor_s = acceptor - QUARTER_SEQ_LEN
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
    main(sys.argv[1:])