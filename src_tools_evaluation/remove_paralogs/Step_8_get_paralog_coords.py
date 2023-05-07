with open("../../Dataset/refseq_all_transcripts.gff", "r") as fr:
    lines = fr.read().splitlines()
    for line in lines:
        print(line)