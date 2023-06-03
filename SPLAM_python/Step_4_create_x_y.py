import os
import sys

def main(argv):
    # defining flanking size
    SEQ_LEN="800"
    HALF_SEQ_LEN = int(SEQ_LEN)//2
    QUARTER_SEQ_LEN = int(SEQ_LEN)//4

    # opening files
    os.makedirs(argv[0]+"/INPUTS/", exist_ok=True)    
    fw = open(argv[0]+"/INPUTS/input.fa", "w")
    fr_donor = open(argv[0]+"/juncs/donor_seq.fa", "r")
    fr_acceptor = open(argv[0]+"/juncs/acceptor_seq.fa", "r")

    # parsing donor and acceptor fa files
    lines_d = fr_donor.read().splitlines()
    lines_a = fr_acceptor.read().splitlines()
    line_num = len(lines_d)

    # initializing stats
    canonical_d_count = 0 # GT
    noncanonical_d_count = 0
    canonical_a_count = 0 # AG
    noncanonical_a_count = 0
    donors = {}
    acceptors = {}
    num_skipped = 0
    for idx in range(0, line_num, 2):
        # PARSE FIRST LINE
        # >chr1:10000-20000(+)
        chr_name = lines_d[idx]
        strand = lines_d[idx][-2]
        chromosome = lines_d[idx].split(":")[0]

        d_splits = lines_d[idx].split(":")[1].split("(")
        d_start, d_end = d_splits[0].split("-")
        d_strand = d_splits[1][0]

        a_splits = lines_a[idx].split(":")[1].split("(")
        a_start, a_end = a_splits[0].split("-")
        a_strand = a_splits[1][0]

        # print("donor   : ", d_start, d_end)
        # print("d_strand: ", d_strand)

        # print("Acceptor   : ", a_start, a_end)
        # print("a_strand: ", a_strand)

        if strand == "+":
            donor_pos = int(d_start) + QUARTER_SEQ_LEN
            acceptor_pos = int(a_end) - QUARTER_SEQ_LEN
        elif strand == "-":
            donor_pos = int(d_end) - QUARTER_SEQ_LEN
            acceptor_pos = int(a_start) + QUARTER_SEQ_LEN

        # PARSE SECOND LINE
        idx += 1

        seq_d = lines_d[idx]
        seq_a = lines_a[idx]
        len_d = len(seq_d)
        len_a = len(seq_a)

        if len_d != len_a:
            print(f'Unequal lengths: seq_d {len_d}, seq_a {len_a}')

        if len_d == HALF_SEQ_LEN and len_a == HALF_SEQ_LEN:
            # combine normally
            x = seq_d + seq_a
        else:
            # pad with repeating N
            x = seq_d + (HALF_SEQ_LEN - len_d) * 'N' + (HALF_SEQ_LEN - len_a) * 'N' + seq_a

        x = x.upper()
        if (len(x) != int(SEQ_LEN)): 
            print("x: ", len(x))
       
        # skip sequence if there are Ns in the sequence
        if x[QUARTER_SEQ_LEN] == "N" or x[QUARTER_SEQ_LEN+1] == "N" or x[QUARTER_SEQ_LEN*3-2] == "N" or x[QUARTER_SEQ_LEN*3-1] == "N":
            num_skipped += 1
            continue
        
        # write the final fasta entry as two lines
        fw.write(chromosome + ";" + str(donor_pos) +";"+ str(acceptor_pos) + ";" + d_strand + ";1\n")
        fw.write(x + "\n")

        # get stats on the dimers 
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

    # output stats
    print("Number of skips due to N in dimer: ", num_skipped)

    print("Canonical donor count: ", canonical_d_count)
    print("Noncanonical donor count: ", noncanonical_d_count)
    print("Canonical acceptor count: ", canonical_a_count)
    print("Noncanonical acceptor count: ", noncanonical_a_count)
    for key, value in donors.items():
        print("Donor   : ", key, " (", value, ")")
    for key, value in acceptors.items():
        print("Acceptor: ", key, " (", value, ")")

    fw.close()
    fr_acceptor.close()
    fr_donor.close()

if __name__ == "__main__":
    main(sys.argv[1:])