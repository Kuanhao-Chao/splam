def chr_name_convert():
    f_chrs = open("../../Dataset/Refseq_2_UCSU_chromosome_names.tsv", "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = eles[1]
    return chrs

# store a list of paralogs here:
paralogs = []
with open("show_coords_paralogs.bed", "r") as fr:
    lines = fr.read().splitlines()
    for line in lines:
        paralogs.append(line)
        print(line)
print(paralogs)


# Read the gtf file (do the filtration here)
fw = open("paralogs.bed", "w")

chrs = chr_name_convert()

with open("../../Dataset/refseq_all_transcripts.gff", "r") as fr:
    lines = fr.read().splitlines()
    for line in lines:
        eles = line.split("\t")
        trans_id = eles[8].split(";")[0][3:]
        if trans_id in paralogs:
            print(eles)
            # eles
            fw.write(chrs[eles[0]]+ "\t" + eles[3] + "\t" + eles[4] + "\t" + trans_id + "\t0\t" + eles[6] + "\n")