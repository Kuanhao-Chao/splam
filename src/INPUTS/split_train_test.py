import os
import sys 

TEST_CHRS = ['chr1', 'chr9']

def split_seq_name(seq):
    return seq.split(":")[0][1:]

def main(arg):
    output_dir = ""
    if arg == "pos_100":
        output_dir = "input_pos_100"
    elif arg == "neg":
        output_dir = "input_neg"
    elif arg == "noncan_neg":
        output_dir = "input_noncan_neg"
    elif arg == "neg_1":
        output_dir = "input_neg_1"


    train_fn = output_dir + "/" + "train_"+arg+".shuffle.fa"
    test_fn = output_dir + "/" + "test_"+arg+".shuffle.fa"

    os.makedirs(output_dir, exist_ok=True)

    train_fw = open(train_fn, "w")
    test_fw = open(test_fn, "w")

    with open("input_"+arg+".shuffle.fa", "r") as f:
        lines = f.read().splitlines()

        idx = 0
        with open("./input_"+arg+".shuffle.fa", "r") as f:
            lines = f.read().splitlines()
            seq_name = ""
            seq = ""
            for line in lines:
                # print(line)
                if idx % 2 == 0:
                    seq_name = split_seq_name(line)
                    print("seq_name: ", seq_name)
                    seq = ">"+seq_name+"\n"
                elif idx % 2 == 1:
                    seq = line+"\n"

                if seq_name in TEST_CHRS:
                    test_fw.write(seq)
                else:
                    train_fw.write(seq)

                idx += 1

    train_fw.close()
    test_fw.close()

    print("idx: ", idx)

if __name__ == "__main__":
    main(sys.argv[1])