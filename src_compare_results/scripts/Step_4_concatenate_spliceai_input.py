import os
import sys

def main(argv):
    SEQ_LEN = "800"
    QUARTER_SEQ_LEN = int(SEQ_LEN) // 4

    # output_files = ["./dataset/pos/", "./dataset/neg_can/", "./dataset/neg_noncan/", "./dataset/outlier_test/"]

    output_filename = f'../output/{argv[0]}/{argv[1]}_spliceai_seq'
    fr = open(output_filename+'.fa', 'r')
    fw_noN = open(output_filename+'_noN.fa', 'w')
    fw_N = open(output_filename+'_N.fa', 'w')

    lines = fr.read().splitlines()
    line_num = len(lines)

    donors_noN = {}
    donors_N = {}
    acceptors_noN = {}
    acceptors_N = {}

    for idx in range(0, line_num, 2):
        # first line
        chr_name = lines[idx]
        strand = lines[idx][-2]
        chromosome = lines[idx].split(":")[0]
        splits = lines[idx].split(":")[1].split("(")
        start, end = splits[0].split("-")
        strand = splits[1][0]

        # reset start and end to original positions
        start = str(int(start) + 5200)
        end = str(int(end) - 5200)

        # second line
        idx += 1

        seq = lines[idx]
        length = len(seq)
        print(f'length: {length}')
        
        # noN method
        x_noN = seq
        x_noN= x_noN.upper()
        fw_noN.write(chromosome+";"+start+";"+end+";"+strand+"\n")
        fw_noN.write(x_noN + "\n")

        # N method
        x_N = 'N'*(5000) + seq[5000:-5000] + 'N'*(5000)
        x_N = x_N.upper()
        fw_N.write(chromosome+";"+start+";"+end+";"+strand+"\n")
        fw_N.write(x_N + "\n")

        # statistical comparison of dimers
        donor_dimer_noN = x_noN[QUARTER_SEQ_LEN+5000:QUARTER_SEQ_LEN+5000+2]
        acceptor_dimer_noN = x_noN[len(x_noN)-QUARTER_SEQ_LEN-5000-2:len(x_noN)-QUARTER_SEQ_LEN-5000]
        donor_dimer_N = x_N[QUARTER_SEQ_LEN+5000:QUARTER_SEQ_LEN+5000+2]
        acceptor_dimer_N = x_N[len(x_N)-QUARTER_SEQ_LEN-5000-2:len(x_N)-QUARTER_SEQ_LEN-5000]

        if donor_dimer_noN not in donors_noN.keys():
            donors_noN[donor_dimer_noN] = 1
        else:
            donors_noN[donor_dimer_noN] += 1
        if acceptor_dimer_noN not in acceptors_noN.keys():
            acceptors_noN[acceptor_dimer_noN] = 1
        else:
            acceptors_noN[acceptor_dimer_noN] += 1
        
        if donor_dimer_N not in donors_N.keys():
            donors_N[donor_dimer_N] = 1
        else:
            donors_N[donor_dimer_N] += 1
        if acceptor_dimer_N not in acceptors_N.keys():
            acceptors_N[acceptor_dimer_N] = 1
        else:
            acceptors_N[acceptor_dimer_N] += 1

    print('--- noN ---')
    for key, value in donors_noN.items():
        print("\tDonor   : ", key, " (", value, ")")
    for key, value in acceptors_noN.items():
        print("\tAcceptor: ", key, " (", value, ")")
    print("\n---- N ----")
    for key, value in donors_N.items():
        print("\tDonor   : ", key, " (", value, ")")
    for key, value in acceptors_N.items():
        print("\tAcceptor: ", key, " (", value, ")")
    print("\n")
    print("\n")
    fw_noN.close()
    fw_N.close()

if __name__ == "__main__":
    main(sys.argv[1:]) 