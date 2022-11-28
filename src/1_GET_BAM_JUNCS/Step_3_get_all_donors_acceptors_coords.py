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
    threshold = "100"
    fw_da = open("../BAM_junctions/"+threshold+"_juncs/d_a.bed", "w")

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
                
            if donor_e >= chrs[chr] or acceptor_e >= chrs[chr]:
                continue
            if donor_s < 0 or acceptor_s < 0:
                continue
            if (strand == "+"):
                fw_da.write(chr + "\t" + str(donor) + "\t" + str(acceptor+1) + "\t" + junc_name + "\t" + score + "\t" + strand + "\n")
            elif (strand == "-"):
                fw_da.write(chr + "\t" + str(acceptor) + "\t" + str(donor+1) + "\t" + junc_name + "\t" + score + "\t" + strand + "\n")

    fw_da.close()

if __name__ == "__main__":
    main()