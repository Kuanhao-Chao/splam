f1 = open("../SPLAM/hg38_chrom_size_refseq.tsv", "r")
f2 = open("../SPLAM/hg38_chrom_size.tsv", "r")


fw = open("chrom_size.txt", "w")


lines = f1.read().splitlines()
for line in lines:
    # print(len(line), line)
    if len(line) > 1:
        print(line)
        key, value = line.split(" ")
        fw.write("tmp = \"" + key + "\";\n")
        fw.write("CHRS[tmp] = " + value + ";\n")

lines = f2.read().splitlines()
for line in lines:
    # print(len(line), line)
    if len(line) > 1:
        print(line)
        key, value = line.split("\t")
        fw.write("tmp = \"" + key + "\";\n")
        fw.write("CHRS[tmp] = " + value + ";\n")