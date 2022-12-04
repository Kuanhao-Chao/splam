import os

chrs = ["chr1",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr2",
        "chr20",
        "chr21",
        "chr22",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chrX",
        "chrY"]
def main():
    fr_d = open('./donor.all.bed', 'r')
    fr_a = open('./acceptor.all.bed', 'r')

    fw_d = open('./donor.bed', 'w')
    fw_a = open('./acceptor.bed', 'w')

    read_idx = 0
    lines_d = fr_d.read().splitlines()
    lines_a = fr_a.read().splitlines()
    Counter = 0
    THRESHOLD = 5000
    for chr_idx in chrs:
        Counter = 0
        while Counter < THRESHOLD:
            # print("read_idx: ", read_idx)
            if lines_d[read_idx].split("\t")[0] != chr_idx:
                read_idx += 1
            else: 
                print("chr_idx: ", chr_idx, "; ", Counter)
                fw_d.write(lines_d[read_idx]+"\n")
                fw_a.write(lines_a[read_idx]+"\n")
                if chr_idx == "chrY":
                    read_idx += 8
                else:
                    read_idx += 30

                Counter += 1

    fr_d.close()
    fr_a.close()
    fw_d.close()
    fw_a.close()

if __name__ == "__main__":
    main()