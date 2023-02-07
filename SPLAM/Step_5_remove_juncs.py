import pysam
import sys
import os

def main(argv):
    SEQ_LEN="800"
    ##############################
    # Read-in invalid junctions
    ##############################
    MODEL_OUTPUT_BASE = "../../results/"+SEQ_LEN+"bp/"+argv[0]+"/OUTPUT/"+argv[1]+"/"
    os.makedirs(MODEL_OUTPUT_BASE, exist_ok=True)
    os.makedirs(MODEL_OUTPUT_BASE+"/BAM/", exist_ok=True)

    threshold = 0.3
    invalid_juncs = set()
    fr = open("../../results/"+SEQ_LEN+"bp/"+argv[0]+"/OUTPUT/"+argv[1]+"/junc_scores.bed", 'r')
    lines = fr.read().splitlines()
    for line in lines:
        chr, start, end, name, score, strand, j_score = line.split("\t")
        if float(j_score) < threshold:
            # print("Detection!!!!!")
            print(">> chr, start, end, strand: ", chr, start, end, strand)
            invalid_juncs.add((chr, start, end, strand))

    # print("invalid_juncs: ", invalid_juncs)
    ##############################
    # Process BAM file
    ##############################
    samfile = pysam.AlignmentFile("../../Dataset/"+argv[0]+"/"+argv[0]+".bam", "rb")
    keep_file = pysam.AlignmentFile(MODEL_OUTPUT_BASE + "/BAM/" + argv[0] + ".cleaned.bam", "wb", template=samfile)
    discard_file = pysam.AlignmentFile(MODEL_OUTPUT_BASE + "/BAM/" + argv[0] + ".discard.bam", "wb", template=samfile)

    for read in samfile.fetch():
        # print(read)
        strand = "?"
        if read.is_forward:
            strand = "+"
        else:
            strand = "-"

        # print(read.reference_start, " - ", read.reference_end)
        # print("read.reference: ", read.reference_name)
        accum_len = read.reference_start
        keep_read = True
        if read.cigartuples == None:
            continue
        for operation, length in read.cigartuples:
            if operation == 3:
                # print("final_junc: ", accum_len, " - ", accum_len+length)
                # print((read.reference_name , accum_len, accum_len+length, strand))
                if ((read.reference_name , str(accum_len), str(accum_len+length), strand) in invalid_juncs):
                    keep_read = False
                    print("Detection!!!!!")
                accum_len += length
            accum_len += length
        if keep_read:
            keep_file.write(read)
        else:
            discard_file.write(read)

            
if __name__ == "__main__":
    main(sys.argv[1:])