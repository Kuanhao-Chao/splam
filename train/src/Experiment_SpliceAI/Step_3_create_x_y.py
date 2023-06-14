import os
import sys

def main(argv):
    SEQ_LEN="600"
    HALF_SEQ_LEN = int(SEQ_LEN)//2
    QUATER_SEQ_LEN = int(SEQ_LEN)//4
    os.makedirs("../../results/spliceAI/"+SEQ_LEN+"bp/"+argv[0]+"/INPUTS/", exist_ok=True)
    fw = open("../../results/spliceAI/"+SEQ_LEN+"bp/"+argv[0]+"/INPUTS/input.fa", "w")
    fr_donor = open("../../results/spliceAI/"+SEQ_LEN+"bp/"+argv[0]+"/juncs/donor_seq.fa", "r")
    fr_acceptor = open("../../results/spliceAI/"+SEQ_LEN+"bp/"+argv[0]+"/juncs/acceptor_seq.fa", "r")

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
            chromosome = lines_d[idx].split(":")[0]
            d_splits = lines_d[idx].split(":")[1].split("(")
            d_start, d_end = d_splits[0].split("-")
            d_strand = d_splits[1][0]

            a_splits = lines_a[idx].split(":")[1].split("(")
            a_start, a_end = a_splits[0].split("-")
            a_strand = a_splits[1][0]

            print("donor   : ", d_start, d_end)
            print("d_strand: ", d_strand)

            print("Acceptor   : ", a_start, a_end)
            print("a_strand: ", a_strand)

            donor_pos = (int(d_start) + int(d_end))//2
            acceptor_pos = (int(a_start) + int(a_end))//2

            fw.write(chromosome + ";" + str(donor_pos) +";"+ str(acceptor_pos) + ";" + d_strand + "\n")
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
                # y = (250, 602)
            else:
                x = seq_d + (HALF_SEQ_LEN - len_d) * 'N' + (HALF_SEQ_LEN - len_a) * 'N' + seq_a
                # y = (250, 602)
            
            print("x: ", len(x))
            print("x: ", len(x))

            x = x.upper()
            if x[QUATER_SEQ_LEN] == "N" or x[QUATER_SEQ_LEN+1] == "N" or x[QUATER_SEQ_LEN*3-1] == "N" or x[QUATER_SEQ_LEN*3] == "N":
                continue

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
    main(sys.argv[1:])