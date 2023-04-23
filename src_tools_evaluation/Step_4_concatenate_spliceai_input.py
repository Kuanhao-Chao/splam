import os

def main():
    SEQ_LEN = "800"
    HALF_SEQ_LEN = int(SEQ_LEN) // 2
    QUATER_SEQ_LEN = int(SEQ_LEN) // 4

    print("QUATER_SEQ_LEN: ", QUATER_SEQ_LEN)

    # output_files = ["./dataset/pos/", "./dataset/neg_can/", "./dataset/neg_noncan/", "./dataset/outlier_test/"]
    output_dir = "./dataset/"
    output_files = [output_dir+"pos/", output_dir+"pos_MANE/", output_dir+"pos_ALTS/", output_dir+"neg_1/", output_dir+"neg_random/"]

    for output_file in output_files:
        print(">> output_file")
        fr = open(output_file+"spliceai/spliceai.juncs.seq.fa", "r")
        fw_noN = open(output_file+"spliceai/spliceai.noN.juncs.seq.fa", "w")
        fw_N = open(output_file+"spliceai/spliceai.N.juncs.seq.fa", "w")

        lines = fr.read().splitlines()
        line_num = len(lines)

        donors_noN = {}
        donors_N = {}
        acceptors_noN = {}
        acceptors_N = {}

        for idx in range(line_num):
            if idx % 2 == 0:
                chr_name = lines[idx]
                strand = lines[idx][-2]
            else:
                seq = lines[idx]
                length = len(seq)
                
                for method in ["noN", "N"]:
                    if method == "noN":
                        x = seq
                    else:
                        x = 'N'*(5000) + seq[5000:-5000] + 'N'*(5000)
                    x = x.upper()
                    # if x[QUATER_SEQ_LEN + 5000] == "N" or x[QUATER_SEQ_LEN+1 + 5000] == "N" or x[len(x) - QUATER_SEQ_LEN-5000-2] == "N" or x[len(x) - QUATER_SEQ_LEN-5000-1] == "N":
                    # # or x[QUATER_SEQ_LEN*3-1 + 5000] == "N" or x[QUATER_SEQ_LEN*3 + 5000] == "N":
                    #     continue
                    
                    eles = lines[idx-1].split(":")
                    chr = eles[0]
                    strand = eles[-1][-2]
                    s_e = eles[1].split("-")
                    start = str(int(s_e[0])+5200)
                    end = str(int(s_e[1].split("(")[0])-5200)
                    # [:-3]
                    # print(">> start - end: ", start, " - ", end)

                    if method == "noN":
                        fw_noN.write(chr+";"+start+";"+end+";"+strand+"\n")
                        fw_noN.write(x + "\n")
                    else:
                        fw_N.write(chr+";"+start+";"+end+";"+strand+"\n")
                        fw_N.write(x + "\n")
                    

                    donor_dimer = x[QUATER_SEQ_LEN + 5000:QUATER_SEQ_LEN+2 + 5000]
                    acceptor_dimer = x[len(x) - QUATER_SEQ_LEN-5000-2:len(x)-QUATER_SEQ_LEN-5000]


                    if method == "noN":
                        if donor_dimer not in donors_noN.keys():
                            donors_noN[donor_dimer] = 1
                        else:
                            donors_noN[donor_dimer] += 1
                        if acceptor_dimer not in acceptors_noN.keys():
                            acceptors_noN[acceptor_dimer] = 1
                        else:
                            acceptors_noN[acceptor_dimer] += 1
                    
                    else:
                        if donor_dimer not in donors_N.keys():
                            donors_N[donor_dimer] = 1
                        else:
                            donors_N[donor_dimer] += 1
                        if acceptor_dimer not in acceptors_N.keys():
                            acceptors_N[acceptor_dimer] = 1
                        else:
                            acceptors_N[acceptor_dimer] += 1
                
        for key, value in donors_noN.items():
            print("\tDonor   : ", key, " (", value, ")")
        for key, value in acceptors_noN.items():
            print("\tAcceptor: ", key, " (", value, ")")
        print("\n")
        for key, value in donors_N.items():
            print("\tDonor   : ", key, " (", value, ")")
        for key, value in acceptors_N.items():
            print("\tAcceptor: ", key, " (", value, ")")
        print("\n")
        print("\n")
        fw_noN.close()
        fw_N.close()

if __name__ == "__main__":
    main()