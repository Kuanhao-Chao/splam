import sys

def main():
    output_files = ["./OUTPUT/pos/", "./OUTPUT/neg_can/", "./OUTPUT/neg_noncan/", "./OUTPUT/neg_1/"]

    prefix = "../src/INPUTS/800bp/input_"
    output_testset_files = ["pos", "neg_can", "neg_noncan", "neg_1"]

    target = "splam"

    for idx in range(1,4):
        output_set = set()
        output_testset_set = set()

        output_file = output_files[idx]
        output_testset_file = output_testset_files[idx]

        ######################################################
        # Reading the SPLAM testing dataset.
        ######################################################
        fr_da = open(prefix+output_testset_file+"/test_"+output_testset_file+".shuffle.fa", "r")
        print("Target file: ", prefix+output_testset_file+"/test_"+output_testset_file+".shuffle.fa")
        lines_da = fr_da.read().splitlines()
        print("lines_da: ", len(lines_da))
        line_num = len(lines_da)
        
        for idx in range(line_num):
            if idx % 2 == 0:
                chr_name = lines_da[idx]
                strand = lines_da[idx][-2]
            else:
                seq = lines_da[idx]
                len_seq = len(seq)
                x = seq.upper()
                output_testset_set.add(x)
        print("\tlen(output_testset_set): ", len(output_testset_set))
        ######################################################
        # Reading the SPLAM junction for comparison dataset
        ######################################################
        fr_da = open(output_file+target+".juncs.seq.fa", "r")
        print("Target file: ", output_file+target+".juncs.seq.fa")
        lines_da = fr_da.read().splitlines()
        print("lines_da: ", len(lines_da))
        line_num = len(lines_da)

        for idx in range(line_num):
            if idx % 2 == 0:
                chr_name = lines_da[idx]
                strand = lines_da[idx][-2]
            else:
                seq = lines_da[idx]
                len_seq = len(seq)
                x = seq.upper()
                output_set.add(x)
        print("\tlen(output_testset_set): ", len(output_set))

        z = output_testset_set.intersection(output_set)
        print("Intersection set size: ", len(z))
        print("\n\n")

if __name__ == "__main__":
    main()