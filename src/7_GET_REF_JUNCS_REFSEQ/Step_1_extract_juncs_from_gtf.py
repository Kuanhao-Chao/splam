import os 
import re
import json 
import torch

def chr_name_convert():
    f_chrs = open("../../Dataset/Refseq_2_UCSU_chromosome_names.tsv", "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = eles[1]
    return chrs


def main():
    JUNC_COUNTER = 0
    os.makedirs("./REF_junctions/", exist_ok=True)
    fw = open("./REF_junctions/ref_d_a.bed", 'w')
    chrs = chr_name_convert()
    with open("../../Dataset/refseq_GCF_000001405.40_GRCh38.p14_genomic.gff", 'r') as f:
    # with open("../../Dataset/MANE.GRCh38.v1.0.ensembl_genomic.gff", 'r') as f:
        lists = f.read().splitlines() 
        transcript_id = ""
        prev_transcript_id = ""
        chr = ""
        prev_chr = ""
        strand = "."
        prev_strand = "."

        starts = []
        ends = []

        is_protein_entries = False
        for ele in lists:
            line = ele.split("\t")
            if len(line) < 8:
                continue

            if (line[2] == "gene"):
                match_prot = re.search(r"gene_biotype=protein_coding", line[8])
                if match_prot is not None:
                    # Following entries are protein coding genes
                    is_protein_entries = True
                else:
                    is_protein_entries = False
            
            # if is_protein_entries:
            #     fw.write(ele+"\n")

            if (line[2] == "exon"):
                match = re.search(r"transcript_id=\w+", line[8])
                if match is not None:
                    transcript_id = match.group()[14:]
                    chr = line[0]
                    strand = line[6]
                    exon_start = int(line[3])
                    exon_end = int(line[4])
                    # print("\tchr: ", chr)
                    # print("\tstrand: ", strand)
                    # print("\texon_start: ", exon_start)
                    # print("\texon_end: ", exon_end)

                    print("transcript_id: ", transcript_id)
                    # print("gene_id: ", gene_id)
                    # print("gene_name: ", gene_name)
                    # print("chr: ", chr)
                    # print("exon_start: ", exon_start)
                    # print("exon_end: ", exon_end)

                    if prev_transcript_id != transcript_id:
                        # print(transcript_id)
                        # print("starts: ", starts)
                        # print("ends  : ",  ends)

                        if prev_strand == '-':
                            starts.reverse()    
                            ends.reverse()
                        
                        starts = starts[1:]
                        ends = ends[:-1]

                        # if (prev_chr in chrs.keys()):
                        for idx in range(len(starts)):
                            JUNC_COUNTER += 1
                            # print(chr + "\t" + str(ends[idx]) + "\t" + str(starts[idx]) + "\t" + "JUNC"+str(JUNC_COUNTER) + "\t1\t" + strand + "\n")
                            # fw.write(chr + "\t" + str(ends[idx]) + "\t" + str(starts[idx]) + "\t" + "JUNC" + "\t1\t" + strand + "\n")
                            fw.write(chrs[prev_chr] + "\t" + str(ends[idx]) + "\t" + str(starts[idx]) + "\t" + "JUNC" + "\t1\t" + prev_strand + "\n")
                            
                        starts.clear()
                        ends.clear()
                    starts.append(exon_start)
                    ends.append(exon_end)

                    prev_transcript_id = transcript_id
                    prev_chr = chr
                    prev_strand = strand
    fw.close()

if __name__ == "__main__":
    main()