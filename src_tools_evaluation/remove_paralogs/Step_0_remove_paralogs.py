import re
# gene_id_ls = []
gene_id_dict = {}
with open("../../Dataset/mart_export.txt") as fr:
    lines = fr.read().splitlines()
    for line in lines:
        gene_id = line.split(",")[0]
        if len(gene_id_dict) == 0:
            gene_id_dict = {gene_id}
        else:
            gene_id_dict.add(gene_id)

print(gene_id_dict)
print(len(gene_id_dict))



# from BCBio import GFF

# in_file = "../../Dataset/refseq_GCF_000001405.40_GRCh38.p14_genomic.gff"

# # in_handle = open(in_file)
# # for rec in GFF.parse(in_handle, target_lines=1000):
# #     print(rec)
# # in_handle.close()

# limit_info = dict(gff_id=["NC_000001.11", "NC_000009.12"])

# in_handle = open(in_file)
# for rec in GFF.parse(in_handle, limit_info=limit_info):
#     print(rec.features[0])
# in_handle.close()



target_paralog_bed = "../../Dataset/paralog.bed"
fw = open(target_paralog_bed, "w")

def find_paralogs(gff_file):
    global Hit
    with open(gff_file, "r") as fr:
        lines = fr.read().splitlines()
        for line in lines:
            entry = line.split("\t") 
            if len(entry) > 5:
                if entry[2] == "gene" or entry[2] == "pseudogene":
                    # if re.search("ID=gene:", line):
                    #     print(line, end='')
                    entry_parse = entry[8].split(";")
                    gene_id = entry_parse[0][8:]
                    if gene_id in gene_id_dict:
                        Hit += 1
                        print("entry_parse: ", entry)
                        fw.write("chr"+entry[0] + "\t" + entry[3] + "\t" + entry[4] + "\t" + gene_id + "\t" + "0" + "\t" + entry[6] + "\n")

Hit = 0
find_paralogs("../../Dataset/Homo_sapiens.GRCh38.109.chromosome.1.gff3")

find_paralogs("../../Dataset/Homo_sapiens.GRCh38.109.chromosome.9.gff3")

print(Hit)
fw.close()