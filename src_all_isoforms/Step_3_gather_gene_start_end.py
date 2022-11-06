import json
def main():
    gene_dict = {}
    with open("../Dataset/all_protein_coding_genes.gtf", "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            line_split = line.split("\t")
            start = int(line_split[3])
            end = int(line_split[4])
            # print(start, end)
            info = line_split[8].split(";")
            gene_id = info[0][9:-1]
            gene_dict [gene_id] = {}
            gene_dict [gene_id]["start"] = start
            gene_dict [gene_id]["end"] = end
        # print(len(gene_dict))
    with open("../Dataset/all_protein_coding_genes.js", "w") as f_out:
        json.dump(gene_dict, f_out)

if __name__ == "__main__":
    main()