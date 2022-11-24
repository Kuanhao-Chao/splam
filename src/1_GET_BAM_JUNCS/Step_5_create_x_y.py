def main():
    fw = open("../INPUTS/input_pos.fa", "w")
    fr_donor = open("../BAM_junctions/donor_seq.fa", "r")
    fr_acceptor = open("../BAM_junctions/acceptor_seq.fa", "r")

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
    for idx in range(line_num):
        if idx % 2 == 0:
            chr_name = lines_d[idx]
            fw.write(lines_d[idx]+"_"+lines_a[idx][1:] + "\n")
        else:
            seq_d = lines_d[idx]
            seq_a = lines_a[idx]
            len_d = len(seq_d)
            len_a = len(seq_a)
            if len_d != len_a:
                print("seq_d: ", len_d)
                print("seq_a: ", len_a)
            if len_d == 400 and len_a == 400:
                x = seq_d + seq_a
                y = (200, 600)
            else:
                x = seq_d + (400 - len_d) * 'N' + (400 - len_a) * 'N' + seq_a
                # print("x: ", len(x))
                y = (200, 600)
            
            fw.write(x + "\n")
            x = x.upper()
            if x[200:202] not in donors.keys():
                donors[x[200:202]] = 1
            else:
                donors[x[200:202]] += 1

            if x[598:600] not in acceptors.keys():
                acceptors[x[598:600]] = 1
            else:
                acceptors[x[598:600]] += 1

            if (x[200:202] == "GT"):
                canonical_d_count += 1
            else:
                noncanonical_d_count += 1
            if (x[598:600] == "AG"):
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