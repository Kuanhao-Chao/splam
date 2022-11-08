import json
import bisect

def main(): 
    chrs_dict = {}
    chrs_list = ["chr"+str(i) for i in range(1,23)]
    chrs_list.append("chrM")
    print("chrs_list: ", chrs_list)
    f_chrs = open("../Dataset/GRCh38_RefSeq2UCSC.txt", "r")
    chrs = f_chrs.read().splitlines()
    for chr in chrs:
        old_chr, new_chr = chr.split("\t")
        chrs_dict[old_chr] = new_chr
    f_chrs.close()

    f_genes = open("../Dataset/all_protein_coding_genes.js", "r")
    genes_s_e = json.load(f_genes)
    f_genes.close()

    gene_dict = {}
    COUNTER = 0
    with open("../Dataset/all_protein_coding_exons.gtf", "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            split_line = line.split("\t")
            if (split_line[0] in chrs_dict.keys()):
                COUNTER += 1
                chr = chrs_dict[split_line[0]]
                start = int(split_line[3])
                end = int(split_line[4])
                strand = split_line[6]
                info = split_line[8].split(";")
                gene_id = info[0][9:-1]
                tranx_id = info[1][16:-1]
                # donor = ""
                # acceptor = ""
                # if strand == "+":
                #     acceptor = start
                #     donor = end
                # elif strand == "-":
                #     acceptor = end
                #     donor = start
                if gene_id not in gene_dict.keys():
                    gene_dict[gene_id] ={}
                gene_dict[gene_id]["chr"] = chr
                gene_dict[gene_id]["strand"] = strand
                if "tranx" not in gene_dict[gene_id]:
                    gene_dict[gene_id]["tranx"] = []    
                if tranx_id not in gene_dict[gene_id]["tranx"]:
                    gene_dict[gene_id]["tranx"].append(tranx_id)


                if "start" not in gene_dict[gene_id]:
                    gene_dict[gene_id]["start"] = []
                # print("gene_dict[gene_id][donor]:" , donor not in gene_dict[gene_id]["donor"])
                if start not in gene_dict[gene_id]["start"]:
                    bisect.insort(gene_dict[gene_id]["start"], start)


                if "end" not in gene_dict[gene_id]:
                    gene_dict[gene_id]["end"] = []    
                if end not in gene_dict[gene_id]["end"]:
                    bisect.insort(gene_dict[gene_id]["end"], end)

            # if COUNTER > 40:
            #     break


    with open("./dataset.txt", "w") as f_out:
        for gene, info in gene_dict.items():
            # print("gene: ", gene)
            # print("info: ", info)

            if gene in genes_s_e.keys():
                if int(genes_s_e[gene]["start"]) > 5000:
                    if info["chr"] in chrs_list:
                        info["start"] = info["start"][1:]
                        info["end"] = info["end"][:-1]

                        # There are splice sites!!
                        if len(info["start"]) > 0 and len(info["end"]) > 0:
                            print_line = gene + "\t0\t" + info["chr"] + "\t" + info["strand"] + "\t" + str(genes_s_e[gene]["start"]) + "\t" + str(genes_s_e[gene]["end"]) + "\t" + ','.join(map(str, info["start"])) + "\t" + ','.join(map(str, info["end"])) + "\n"
                f_out.write(print_line)
        

if __name__ == "__main__":
    main()