from Bio import SeqIO

def main():
    # fw = open("tiebrush_chr22_5.input.fa", "w")
    hg38 = "../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa"
    fr_hg38 = open(hg38, "r")
    output_file = "hg38.fa"
    # fr_acceptor = open("tiebrush_chr22_5.junc_seq.acceptor.fa", "r")

    fasta_sequences = SeqIO.parse(open(hg38),'fasta')
    with open(output_file) as out_file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            print(name)
            print(sequence)
            # new_sequence = some_function(sequence)
            # write_fasta(out_file)


    # lines_d = fr_donor.read().splitlines()
    # lines_a = fr_acceptor.read().splitlines()

    # line_num = len(lines_d)

    # for idx in range(line_num):
    #     if idx % 2 == 0:
    #         pass
    #         fw.write(lines_d[idx] + "\n")
    #     else:
    #         seq_d = lines_d[idx]
    #         seq_a = lines_a[idx]
    #         len_d = len(seq_d)
    #         len_a = len(seq_a)
    #         print(len_d)
    #         if len_d == 400 and len_a == 400:
    #             x = seq_d + seq_a
    #             y = (200, 600)
    #             # print(x)
    #         else:
    #             x = seq_d + (400 - len_d) * 'N' + (400 - len_a) * 'N' + seq_a
    #             y = (200, 600)
    #         fw.write(x + "\n")
    #         # print(x)
    #         # print(len(x))
    #         # print(y)
    # fw.close()
if __name__ == "__main__":
    main()