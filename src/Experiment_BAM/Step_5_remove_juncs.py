import pysam

def main():
    ##############################
    # Read-in invalid junctions
    ##############################
    MODEL_OUTPUT_BASE = "./OUTPUT/SpliceAI_6_RB_p1_n1_nn1_TB_all_samples_thr_100_v8/"
    invalid_juncs = set()
    fr = open("./OUTPUT/SpliceAI_6_RB_p1_n1_nn1_TB_all_samples_thr_100_v8/removed_junc.bed", 'r')
    lines = fr.read().splitlines()
    for line in lines:
        chr, start, end, name, score, strand = line.split("\t")
        invalid_juncs.add((chr, start, end, strand))

    print("invalid_juncs: ", invalid_juncs)
    ##############################
    # Process BAM file
    ##############################
    samfile = pysam.AlignmentFile("../../Dataset/SRR1352129_chr22.bam", "rb")
    keep_file = pysam.AlignmentFile(MODEL_OUTPUT_BASE + "SRR1352129_chr22.keep.bam", "wb", template=samfile)
    discard_file = pysam.AlignmentFile(MODEL_OUTPUT_BASE + "SRR1352129_chr22.discard.bam", "wb", template=samfile)

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
    main()