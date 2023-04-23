import sys

def main():
    output_dir = "./dataset/"
    output_files = [output_dir+"pos/", output_dir+"pos_MANE/", output_dir+"pos_ALTS/", output_dir+"neg_1/", output_dir+"neg_random/"]

    target = "splam"
    for output_file in output_files:
        for target in ["spliceai/spliceai.noN", "spliceai/spliceai.N", "splam/splam"]:
            print("####### target: ", target)
            fr_da = open(output_file+target+".juncs.seq.fa", "r")
            print("Target file: ", output_file+target+".juncs.seq.fa")
            lines_da = fr_da.read().splitlines()
            print("lines_da: ", len(lines_da))
            
            line_num = len(lines_da)

            canonical_d_count = 0
            noncanonical_d_count = 0
            canonical_a_count = 0
            noncanonical_a_count = 0

            donors = {}
            acceptors = {}
            chr_name = ""
            strand = ""
            for idx in range(line_num):
                if idx % 2 == 0:
                    chr_name = lines_da[idx]
                    strand = lines_da[idx][-2]
                else:
                    seq = lines_da[idx]
                    len_seq = len(seq)
                    x = seq.upper()
                    # if x[5100] == "N" or x[QUATER_SEQ_LEN+1] == "N" or x[QUATER_SEQ_LEN*3-1] == "N" or x[QUATER_SEQ_LEN*3] == "N":
                    #     continue

                    if (target == "spliceai/spliceai.noN"):
                        # donor_dimer = x[5100-1:5100+1]
                        donor_dimer = x[5200:5200+2]
                        acceptor_dimer = x[len_seq-5200-2:len_seq-5200]

                    elif (target == "spliceai/spliceai.N"):
                        donor_dimer = x[5200:5200+2]
                        acceptor_dimer = x[len_seq-5200-2:len_seq-5200]

                    elif (target == "splam/splam"):
                        donor_dimer = x[200:200+2]
                        acceptor_dimer = x[len_seq-200-2:len_seq-200]

                    # if strand == "-":
                    if donor_dimer not in donors.keys():
                        donors[donor_dimer] = 1
                    else:
                        donors[donor_dimer] += 1

                    if acceptor_dimer not in acceptors.keys():
                        acceptors[acceptor_dimer] = 1
                    else:
                        acceptors[acceptor_dimer] += 1

                    if (donor_dimer == "GT"):
                        canonical_d_count += 1
                    else:
                        noncanonical_d_count += 1
                    if (acceptor_dimer == "AG"):
                        canonical_a_count += 1
                    else:
                        noncanonical_a_count += 1
            print("Target: ", target)
            print("\tCanonical donor count: ", canonical_d_count)
            print("\tNoncanonical donor count: ", noncanonical_d_count)

            print("\tCanonical acceptor count: ", canonical_a_count)
            print("\tNoncanonical acceptor count: ", noncanonical_a_count)

            for key, value in donors.items():
                print("\tDonor   : ", key, " (", value, ")")
            for key, value in acceptors.items():
                print("\tAcceptor: ", key, " (", value, ")")
if __name__ == "__main__":
    main()