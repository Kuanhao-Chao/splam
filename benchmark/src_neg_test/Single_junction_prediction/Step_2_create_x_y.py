def main():
    SEQ_LEN = 400
    fw = open("./input.fa", "w")
    fr_donor = open("./donor_seq.fa", "r")
    fr_acceptor = open("./acceptor_seq.fa", "r")

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
            fw.write(lines_d[idx]+"_"+lines_a[idx][1:] + "\n")
        else:
            seq_d = lines_d[idx]
            seq_a = lines_a[idx]
            len_d = len(seq_d)
            len_a = len(seq_a)

            if len_d != len_a:
                print("seq_d: ", len_d)
                print("seq_a: ", len_a)
            if len_d == SEQ_LEN and len_a == SEQ_LEN:
                x = seq_d + seq_a
                # y = (200, 602)
            else:
                x = seq_d + (SEQ_LEN - len_d) * 'N' + (SEQ_LEN - len_a) * 'N' + seq_a
                # y = (200, 602)
            
            # print("x: ", len(x))
            # print("x: ", len(x))

            x = x.upper()
            if x[200] == "N" or x[201] == "N" or x[599] == "N" or x[600] == "N":
                continue

            fw.write(x + "\n")

            donor_dimer = x[200:202]
            acceptor_dimer = x[598:600]


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
    print("Canonical donor count: ", canonical_d_count)
    print("Noncanonical donor count: ", noncanonical_d_count)

    print("Canonical acceptor count: ", canonical_a_count)
    print("Noncanonical acceptor count: ", noncanonical_a_count)

    for key, value in donors.items():
        print("Donor   : ", key, " (", value, ")")
    for key, value in acceptors.items():
        print("Acceptor: ", key, " (", value, ")")
    fw.close()
if __name__ == "__main__":
    main()