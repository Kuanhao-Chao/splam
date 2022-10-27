import os 
import re
import json 
def main():
    mapping_f = open("../Dataset/mapping.json")
    mapping = json.load(mapping_f)
    print(mapping)
    with open("../Dataset/annotations_with_tranx_id.gff", 'r') as f:
        lists = f.read().splitlines() 

        transcript_id = ""
        transcript_chr = ""
        transcript_start = 0
        transcript_end = 0
        transcript_strand = "."
        # prev_transcript_id = ""

        WRITE = False
        encoding_ls = [0]
        for line in lists:
            line = line.split("\t")
            # print(len(line))
            if len(line) < 8:
                continue

            if (line[2] == "gene"):
                print(line)
            #     WRITE = True
            #     s = ''.join(str(x) for x in encoding_ls)
            #     if (transcript_id != ""):
            #         f_w = open("../Dataset/tranx_label/"+transcript_id+".txt", "w")
            #         f_w.write(s)
            #         f_w.close()

            #         f_w = open("../Dataset/tranx_label/pos_"+transcript_id+".txt", "w")
            #         f_w.write(transcript_chr + "\t" + str(transcript_start) + "\t" + str(transcript_end))
            #         f_w.close()
            #     pattern = ";transcript_id=[\w]*\.[\w]"
            #     x = re.search(pattern, line[8])
            #     # prev_transcript_id = transcript_id
            #     transcript_id = x.group(0)[15:]
            #     transcript_chr = line[0]
            #     transcript_start = int(line[3])
            #     transcript_end = int(line[4])
            #     encoding_ls = [mapping["None"]]*(transcript_end-transcript_start+1)
            #     transcript_strand = line[6]

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