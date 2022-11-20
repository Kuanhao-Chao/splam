import os 
import re
import json 
import torch

def main():
    mapping_f = open("../Dataset/mapping.json")
    mapping = json.load(mapping_f)
    mapping_f.close()
    print(mapping)
    gene_id_2_name = {}
    gene_id_2_features = {}
    with open("../Dataset/hg38c_protein_and_lncRNA.gtf", 'r') as f:
        lists = f.read().splitlines() 

        transcript_id = ""
        gene_id = ""
        prev_gene_id = ""
        gene_name = ""
        prev_gene_name = ""
        chr = ""
        prev_chr = ""
        # transcript_start = 0
        # transcript_end = 0
        transcript_strand = "."
        boundaries = set()

        WRITE = False
        encoding_ls = [0]
        for line in lists:
            line = line.split("\t")
            # print(len(line))
            if len(line) < 8:
                continue

            if (line[2] == "exon"):
            #     WRITE = True
            #     s = ''.join(str(x) for x in encoding_ls)
            #     if (transcript_id != ""):
            #         f_w = open("../Dataset/tranx_label/"+transcript_id+".txt", "w")
            #         f_w.write(s)
            #         f_w.close()

            #         f_w = open("../Dataset/tranx_label/pos_"+transcript_id+".txt", "w")
            #         f_w.write(transcript_chr + "\t" + str(transcript_start) + "\t" + str(transcript_end))
            #         f_w.close()
                # pattern = "transcript_id \"[\w]*\";"
                # x = re.search(pattern, line[8])
                # # prev_transcript_id = transcript_id
                # transcript_id = x.group(0)[15:-2]
                # # print(transcript_id)

                # pattern = "gene_id \"[\w]*\";"
                # x = re.search(pattern, line[8])
                # # prev_transcript_id = transcript_id
                # gene_id = x.group(0)[9:-2]
                # # print(gene_id)

                # pattern = "gene_name \"[\w]*\";"
                # x = re.search(pattern, line[8])
                # # prev_transcript_id = transcript_id
                # # gene_name = x.group(0)[7:-2]
                # # print(gene_name)

                features = line[8].split(";")
                transcript_id = features[0][15:-1]
                gene_id = features[1][10:-1]
                gene_name = features[2][12:-1]
                chr = line[0]
                # print(chr)
                exon_start = int(line[3])
                exon_end = int(line[4])
                # print(exon_start, " ", exon_end)

                if prev_gene_id != gene_id and prev_gene_id != "":

                    # time to write out the label encoding.
                    print(prev_gene_id)
                    boundaries_sort = sorted(list(boundaries))
                    gene_start = boundaries_sort[0]
                    gene_end = boundaries_sort[-1]
                    encoding_ls = [mapping["None"]]*(gene_end - gene_start + 1)
                    # print(encoding_ls)

                    for ele in boundaries_sort:
                        encoding_ls[ele - gene_start] = mapping["Boundary"]
                    boundaries.clear()

                    gene_id_2_name[prev_gene_id] = prev_gene_name
                    gene_features = {}
                    gene_features["bam"] ="../Dataset/tranx_bam/" + prev_gene_id + ".bam"
                    gene_features["feature"] ="../Dataset/tranx_feature/" + prev_gene_id + ".pt"
                    gene_features["label"] ="../Dataset/tranx_label/" + prev_gene_id + ".pt"
                    gene_features["chr"] = prev_chr
                    gene_features["gene_start"] = gene_start
                    gene_features["gene_end"] = gene_end
                    gene_features["region"] ="../Dataset/tranx_label_raw/pos_" + prev_gene_id + ".txt"
                    gene_features["raw_label"] ="../Dataset/tranx_label_raw/" +prev_gene_id + ".txt"
                    gene_id_2_features[prev_gene_id] = gene_features

                    f_w = open("../Dataset/tranx_label_raw/"+prev_gene_id+".txt", "w")
                    s = ''.join(str(x) for x in encoding_ls)
                    f_w.write(s)
                    f_w.close()

                    encoding_ls_tensor = torch.Tensor(encoding_ls)
                    torch.save(encoding_ls_tensor, "../Dataset/tranx_label/" + prev_gene_id + ".pt")

                    f_w = open("../Dataset/tranx_label_raw/pos_"+prev_gene_id+".txt", "w")
                    f_w.write(prev_chr + "\t" + str(gene_start) + "\t" + str(gene_end))
                    f_w.close()

                
                prev_gene_id = gene_id
                prev_gene_name = gene_name
                prev_chr = chr
                boundaries.add(exon_start)
                boundaries.add(exon_end)

        

        # time to write out the label encoding.
        print(prev_gene_id)
        boundaries_sort = sorted(list(boundaries))
        gene_start = boundaries_sort[0]
        gene_end = boundaries_sort[-1]
        encoding_ls = [mapping["None"]]*(gene_end - gene_start + 1)
        # print(encoding_ls)

        for ele in boundaries_sort:
            encoding_ls[ele - gene_start] = mapping["Boundary"]
        boundaries.clear()

        gene_id_2_name[prev_gene_id] = prev_gene_name
        gene_features = {}
        gene_features["bam"] ="../Dataset/tranx_bam/" + prev_gene_id + ".bam"
        gene_features["feature"] ="../Dataset/tranx_feature/" + prev_gene_id + ".pt"
        gene_features["label"] ="../Dataset/tranx_label/" + prev_gene_id + ".pt"
        gene_features["chr"] = prev_chr
        gene_features["gene_start"] = gene_start
        gene_features["gene_end"] = gene_end
        gene_features["region"] ="../Dataset/tranx_label_raw/pos_" + prev_gene_id + ".txt"
        gene_features["raw_label"] ="../Dataset/tranx_label_raw/" + prev_gene_id + ".txt"
        gene_id_2_features[prev_gene_id] = gene_features

        f_w = open("../Dataset/tranx_label_raw/"+prev_gene_id+".txt", "w")
        s = ''.join(str(x) for x in encoding_ls)
        f_w.write(s)
        f_w.close()

        encoding_ls_tensor = torch.Tensor(encoding_ls)
        torch.save(encoding_ls_tensor, "../Dataset/tranx_label/" + prev_gene_id + ".pt")

        f_w = open("../Dataset/tranx_label_raw/pos_"+prev_gene_id+".txt", "w")
        f_w.write(prev_chr + "\t" + str(gene_start) + "\t" + str(gene_end))
        f_w.close()

    with open('../Dataset/geneid_2_name.txt', 'w') as convert_file:
        convert_file.write(json.dumps(gene_id_2_name))

    with open('../Dataset/geneid_2_features.txt', 'w') as convert_file:
        convert_file.write(json.dumps(gene_id_2_features))
                # encoding_ls = [mapping["None"]]*(transcript_end-transcript_start+1)

            # if (line[2] != "mRNA" and line[2] != "exon"):
            #     WRITE = False

            # if (line[2] == "exon") and WRITE:
            #     pattern_e = ";transcript_id=[\w]*\.[\w]"
            #     x_e = re.search(pattern_e, line[8])
            #     transcript_id_e = x_e.group(0)[15:]
            #     if transcript_id == transcript_id_e:
            #         # print(line)
            #         # print("transcript_id_e: ", transcript_id_e)
            #         # print("transcript_id: ", transcript_id)
            #         exon_start = int(line[3])
            #         exon_end = int(line[4])

            #         if (transcript_strand == "+"):
            #             encoding_ls[exon_start - transcript_start] = mapping["Acceptor_pos"]
            #             encoding_ls[exon_end - transcript_start] = mapping["Doner_pos"]
            #         elif (transcript_strand == "-"):
            #             encoding_ls[exon_start - transcript_start] = mapping["Doner_neg"]
            #             encoding_ls[exon_end - transcript_start] = mapping["Acceptor_neg"]
            #     else:
            #         # Should not enter here!
            #         print(">>>>")

if __name__ == "__main__":
    main()