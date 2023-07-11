import os

SEQ_LEN = "800"
HALF_SEQ_LEN = int(SEQ_LEN) // 2
QUARTER_SEQ_LEN = int(SEQ_LEN) // 4

def main(db):
    os.makedirs(f'./4_output/{SEQ_LEN}bp/{db}/', exist_ok=True)
    fw = open(f'./4_output/{SEQ_LEN}bp/{db}/input_neg_random.fa', "w")
    fw_d_a = open(f"./4_output/{SEQ_LEN}bp/{db}/d_a_type.tsv", "w")
    fr_donor = open(f"./3_output/{SEQ_LEN}bp/{db}/donor_seq.fa", "r")
    fr_acceptor = open(f"./3_output/{SEQ_LEN}bp/{db}/acceptor_seq.fa", "r")

    lines_d = fr_donor.read().splitlines()
    lines_a = fr_acceptor.read().splitlines()

    line_num = min(len(lines_d), len(lines_a))

    canonical_d_count = 0
    noncanonical_d_count = 0
    canonical_a_count = 0
    noncanonical_a_count = 0

    donors = {}
    acceptors = {}
    chr_name = ""
    for idx in range(line_num):
        if idx % 2 == 0:
            # getting the header line
            d_e = lines_d[idx][1:].split(':')
            d_l = d_e[1].split('(')
            n1 = d_l[0].split('-')
            n2 = lines_a[idx].split(':')[1].split('(')[0].split('-')

            d_n = int((int(n1[0]) + int(n1[1])) / 2)
            a_n = int((int(n2[0]) + int(n2[1])) / 2)
            chr = d_e[0]
            strand = d_l[1][0]

            # final header format
            chr_name = f'>{chr};{d_n};{a_n};{strand};0\n'
        else:
            # getting the full sequence from both donor and acceptor fa files
            seq_d = lines_d[idx]
            seq_a = lines_a[idx]
            len_d = len(seq_d)
            len_a = len(seq_a)

            if len_d != len_a:
                print("seq_d: ", len_d)
                print("seq_a: ", len_a)
            if len_d == HALF_SEQ_LEN and len_a == HALF_SEQ_LEN:
                x = seq_d + seq_a
                # y = (200, 602)
            else:
                x = seq_d + (HALF_SEQ_LEN - len_d) * 'N' + (HALF_SEQ_LEN - len_a) * 'N' + seq_a
                # y = (200, 602)
            
            x = x.upper()
            if x[QUARTER_SEQ_LEN] == "N" or x[QUARTER_SEQ_LEN+1] == "N" or x[QUARTER_SEQ_LEN*3-1] == "N" or x[QUARTER_SEQ_LEN*3] == "N":
                continue

            fw.write(chr_name)
            fw.write(x + "\n")

            donor_dimer = x[QUARTER_SEQ_LEN:QUARTER_SEQ_LEN+2]
            acceptor_dimer = x[QUARTER_SEQ_LEN*3-2:QUARTER_SEQ_LEN*3]

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

    donors = sorted(donors.items(), key=lambda x: x[1], reverse=True)
    acceptors = sorted(acceptors.items(), key=lambda x: x[1], reverse=True)
    
    
    for key, value in donors[:5]:
        print("Donor   : ", key, " (", str(value), ")")
        fw_d_a.write("Donor   : " + key + "\t" + str(value) + "\n")
    for key, value in acceptors[:5]:
        print("Acceptor: ", key, " (", str(value), ")")
        fw_d_a.write("Acceptor: " + key + "\t" + str(value) + "\n")

    fw_d_a.close()
    fw.close()
    fr_donor.close()
    fr_acceptor.close()

if __name__ == "__main__":

    if os.getcwd() != 'src_neg_test':
        os.chdir('/home/smao10/SPLAM/benchmark/src_neg_test')

    datasets = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    idxs = [0,1,2] #CHANGEME

    for idx in idxs:
        main(datasets[idx])