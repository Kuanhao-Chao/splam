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
    JUNC_COUNTER = 0
    fw = open("./ref_d_a.bed", 'w')
    with open("../Dataset/hg38c_protein_and_lncRNA.gtf", 'r') as f:
        lists = f.read().splitlines() 

        transcript_id = ""
        prev_transcript_id = ""
        gene_id = ""
        prev_gene_id = ""
        gene_name = ""
        prev_gene_name = ""
        chr = ""
        prev_chr = ""
        strand = "."
        boundaries = set()

        starts = []
        ends = []

        WRITE = False
        encoding_ls = [0]
        for line in lists:
            line = line.split("\t")
            if len(line) < 8:
                continue

            if (line[2] == "exon"):

                features = line[8].split(";")
                # print("features: ", features)
                transcript_id = features[0][15:-1]
                gene_id = features[1][10:-1]
                gene_name = features[2][12:-1]
                chr = line[0]
                strand = line[6]
                exon_start = int(line[3])
                exon_end = int(line[4])
                
                # print("transcript_id: ", transcript_id)
                # print("gene_id: ", gene_id)
                # print("gene_name: ", gene_name)
                # print("chr: ", chr)
                # print("exon_start: ", exon_start)
                # print("exon_end: ", exon_end)

                if prev_transcript_id != transcript_id:
                    # print(transcript_id)
                    # print("starts: ", starts)
                    # print("ends  : ",  ends)
                    starts = starts[1:]
                    ends = ends[:-1]
                    for idx in range(len(starts)):
                        JUNC_COUNTER += 1
                        # print(chr + "\t" + str(ends[idx]) + "\t" + str(starts[idx]) + "\t" + "JUNC"+str(JUNC_COUNTER) + "\t1\t" + strand + "\n")
                        fw.write(chr + "\t" + str(ends[idx]) + "\t" + str(starts[idx]) + "\t" + "JUNC" + "\t1\t" + strand + "\n")
                        
                    starts.clear()
                    ends.clear()

                starts.append(exon_start)
                ends.append(exon_end)

                prev_transcript_id = transcript_id
                prev_gene_id = gene_id
                prev_gene_name = gene_name
                prev_chr = chr
    fw.close()
if __name__ == "__main__":
    main()