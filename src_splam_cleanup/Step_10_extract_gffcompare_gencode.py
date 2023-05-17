import re
import sys
import os

# polyA
colname = ["", "base_s", "base_p", "exon_s", "exon_p", "intron_s", "intron_p", "intron_chain_s", "intron_chain_p", "transcript_s", "transcript_p", "locus_s", "locus_p", "matching_intron_chain", "matching_transcripts", "matching_loci", "missed_exons", "novel_exons", "missed_introns", "novel_introns", "missed_loci", "novel_loci"]

for annotation in ["chess",  "gencode",  "refseq_ucsc"]:
    print(">> Processing annotation: ", annotation)
    base_num = 17
    if annotation == "chess":
        base_num = 18
    for condition in ["BEFORE", "AFTER"]:
        # PolyA result
        dir_name = "../results/polyA/assembly/"+annotation+"/"
        os.makedirs(dir_name, exist_ok=True)
        filename = dir_name + condition+".tsv"

        with open(filename, "w") as fw:
            fw.write("\t".join(colname) + "\n")

            for sample in ["R2826", "R2835",  "R2839",  "R2845",  "R2855",  "R2857",  "R2869",  "R2874",  "R2894",  "R2895"]:
                file_name = "../results/polyA/"+sample+"/gffcompare/"+annotation+"/"+condition+"/res.stats"

                print("file_name: ", file_name)
                with open(file_name, 'r') as file:
                    input_text = file.read()

                # Extract all numbers using regular expressions
                numbers = re.findall(r"\d+\.\d+|\d+", input_text)
                print("numbers: ", numbers)

                # Print the extracted numbers
                stats = numbers[base_num:base_num+15]
                stats.append(numbers[base_num+17])
                stats.append(numbers[base_num+20])
                stats.append(numbers[base_num+23])
                stats.append(numbers[base_num+26])
                stats.append(numbers[base_num+29])
                stats.append(numbers[base_num+32])
                stats.insert(0, sample)
                print("numbers: ", stats)
                print(len(stats))
                fw.write("\t".join(stats) + "\n")

        # ribozero result
        dir_name = "../results/ribozero/assembly/"+annotation+"/"
        os.makedirs(dir_name, exist_ok=True)
        filename = dir_name + condition+".tsv"

        with open(filename, "w") as fw:
            fw.write("\t".join(colname) + "\n")

            for sample in ["R12258", "R12260", "R12263", "R12265", "R12266", "R12277", "R12278", "R12280", "R12285", "R12287"]:
                file_name = "../results/ribozero/"+sample+"/gffcompare/"+annotation+"/"+condition+"/res.stats"

                print("file_name: ", file_name)
                with open(file_name, 'r') as file:
                    input_text = file.read()

                # Extract all numbers using regular expressions
                numbers = re.findall(r"\d+\.\d+|\d+", input_text)
                print("numbers: ", numbers)
                # Print the extracted numbers
                stats = numbers[base_num:base_num+15]
                stats.append(numbers[base_num+17])
                stats.append(numbers[base_num+20])
                stats.append(numbers[base_num+23])
                stats.append(numbers[base_num+26])
                stats.append(numbers[base_num+29])
                stats.append(numbers[base_num+32])
                stats.insert(0, sample)
                print("numbers: ", stats)
                print(len(stats))
                fw.write("\t".join(stats) + "\n")
    print("\n\n")