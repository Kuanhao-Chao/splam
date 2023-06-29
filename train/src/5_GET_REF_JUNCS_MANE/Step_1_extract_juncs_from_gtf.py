import os 
import re
import json 
import torch

def main():
    JUNC_COUNTER = 0
    os.makedirs("./REF_junctions/", exist_ok=True)
    fw = open("./REF_junctions/ref_d_a.bed", 'w')
    with open("../../../Dataset/MANE.GRCh38.v1.0.ensembl_genomic.gff", 'r') as f:
        lists = f.read().splitlines() 
        transcript_id = ""
        prev_transcript_id = ""
        chr = ""
        prev_chr = ""
        strand = "."
        prev_strand = "."
        starts = []
        ends = []
        for line in lists:
            line = line.split("\t")
            if len(line) < 8:
                continue

            if (line[2] == "exon"):
                match = re.search(r"transcript_id=\w+", line[8])
                if match is not None:
                    transcript_id = match.group()[14:]
                    chr = line[0]
                    strand = line[6]
                    exon_start = int(line[3])
                    exon_end = int(line[4])
                    print("transcript_id: ", transcript_id)
                    if prev_transcript_id != transcript_id:

                        if prev_strand == '-':
                            starts.reverse()    
                            ends.reverse()
                        
                        starts = starts[1:]
                        ends = ends[:-1]

                        for idx in range(len(starts)):
                            JUNC_COUNTER += 1
                            fw.write(prev_chr + "\t" + str(ends[idx]) + "\t" + str(starts[idx]) + "\t" + "JUNC" + "\t1\t" + prev_strand + "\n")
                            
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