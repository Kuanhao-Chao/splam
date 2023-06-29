import os

def main():
    SEQ_LEN = "800"
    HALF_SEQ_LEN = int(SEQ_LEN) // 2
    QUATER_SEQ_LEN = int(SEQ_LEN) // 4

    print("QUATER_SEQ_LEN: ", QUATER_SEQ_LEN)

    # output_files = ["./dataset/pos/", "./dataset/outlier_test/", "./dataset/neg_5/"]

    output_dir = "./dataset/"
    output_files = [output_dir+"pos/", output_dir+"pos_MANE/", output_dir+"pos_ALTS/", output_dir+"neg_1/", output_dir+"neg_random/"]

    for output_file in output_files:
        print("output_file: ", output_file)
        fw = open(output_file+"splam/splam.juncs.seq.fa", "w")
        fr_donor = open(output_file+"splam/splam.juncs.donor.seq.fa", "r")
        fr_acceptor = open(output_file+"splam/splam.juncs.acceptor.seq.fa", "r")

        lines_d = fr_donor.read().splitlines()
        lines_a = fr_acceptor.read().splitlines()

        line_num = len(lines_d)

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
                chr_name = lines_d[idx]
                strand = lines_d[idx][-2]
            else:
                seq_d = lines_d[idx]
                seq_a = lines_a[idx]
                len_d = len(seq_d)
                len_a = len(seq_a)

                if len_d != len_a:
                    print("seq_d: ", len_d)
                    print("seq_a: ", len_a)
                if len_d == HALF_SEQ_LEN and len_a == HALF_SEQ_LEN:
                    x = seq_d + seq_a
                    # y = (250, 750)
                else:
                    x = seq_d + (HALF_SEQ_LEN - len_d) * 'N' + (HALF_SEQ_LEN - len_a) * 'N' + seq_a
                    # y = (250, 750)
                x = x.upper()
                if len(x) != 800:
                    print(len(x))
                # if x[QUATER_SEQ_LEN] == "N" or x[QUATER_SEQ_LEN+1] == "N" or x[QUATER_SEQ_LEN*3-1] == "N" or x[QUATER_SEQ_LEN*3] == "N":
                #     print(chr_name)
                #     continue

                fw.write(lines_d[idx-1]+"_"+lines_a[idx-1][1:] + "\n")
                fw.write(x + "\n")

                donor_dimer = x[QUATER_SEQ_LEN:QUATER_SEQ_LEN+2]
                acceptor_dimer = x[QUATER_SEQ_LEN*3-2:QUATER_SEQ_LEN*3]

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
        print("Simple line counter: ", idx)
        print("Canonical donor count: ", canonical_d_count)
        print("Noncanonical donor count: ", noncanonical_d_count)

        print("Canonical acceptor count: ", canonical_a_count)
        print("Noncanonical acceptor count: ", noncanonical_a_count)

        for key, value in donors.items():
            print("\tDonor   : ", key, " (", value, ")")
        for key, value in acceptors.items():
            print("\tAcceptor: ", key, " (", value, ")")
        print("\n\n")
        fw.close()
if __name__ == "__main__":
    main()