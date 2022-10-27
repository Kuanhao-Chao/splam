import os 
import pysam
import numpy as np
import json
import torch
import sys, getopt


def main(argv):
    # files = os.listdir("../Dataset/tranx_bam/")
    # for file in files:
    #     print(file[-4:])
    #     if file[-4:] == ".bam":
    #         print(file)
    #         samfile = pysam.AlignmentFile("../Dataset/tranx_bam/"+file, "rb")
    #         print(samfile)

    gene_id = argv[0]
    samfile = pysam.AlignmentFile("../Dataset/tranx_bam/"+gene_id+".bam", "rb")
    # iter = samfile.pileup('chr1', 10, 110630211)110630211
    features_f = open("../Dataset/geneid_2_features.txt")
    features = json.load(features_f)
    features_f.close()

    input_tensor = []

    AVG_COV = 40
    AVG_QUA = 50

    print("gene_id: ", gene_id)
    # print("features[gene_id]['feature']: ", features[gene_id]['feature'])
    # print("features[gene_id]['gene_start']: ", features[gene_id]['gene_start'])
    # print("features[gene_id]['gene_end']: ", features[gene_id]['gene_end'])

    for pileupcolumn in samfile.pileup(features[gene_id]['chr'], features[gene_id]['gene_start'], features[gene_id]['gene_end']):
        if pileupcolumn.reference_pos < features[gene_id]['gene_start']:
            continue
        if pileupcolumn.reference_pos > features[gene_id]['gene_end']:
            continue

        # print("\ncoverage at base %s = %s" % (pileupcolumn.reference_pos, pileupcolumn.nsegments))
        column_encoding = np.zeros((4, 6))
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # query position is None if is_del or is_refskip is set.
                # print("pileupread.query_position: ", pileupread.query_position)
                # print('\tbase in read %s = %s (%s) %s' %
                #     (pileupread.alignment.query_name,
                #     pileupread.alignment.query_sequence[pileupread.query_position], pileupread.alignment.is_forward, 
                #     pileupread.alignment.query_qualities[pileupread.query_position]))

                encoding_idx_offset = 0
                col_seq_idx = 0
                if not pileupread.alignment.is_forward:
                    encoding_idx_offset = 3
                if pileupread.alignment.query_sequence[pileupread.query_position] == "A":
                    col_seq_idx = 0
                elif pileupread.alignment.query_sequence[pileupread.query_position] == "C":
                    col_seq_idx = 1
                elif pileupread.alignment.query_sequence[pileupread.query_position] == "G":
                    col_seq_idx = 2
                elif pileupread.alignment.query_sequence[pileupread.query_position] == "T":
                    col_seq_idx = 3
                else:
                    col_seq_idx = 4                    
                # print("col_seq_idx: ", col_seq_idx)

                # print("col_seq_idx: ", col_seq_idx,  " 0+encoding_idx_offset: ", 0+encoding_idx_offset)
                column_encoding[col_seq_idx][0+encoding_idx_offset] += 1
                # print(column_encoding)
                column_encoding[col_seq_idx][1+encoding_idx_offset] += 1
                # print("pileupread.alignment.query_alignment_qualities[pileupread.query_position]: ", pileupread.alignment.query_alignment_qualities[pileupread.query_position])
                column_encoding[col_seq_idx][2+encoding_idx_offset] += pileupread.alignment.query_qualities[pileupread.query_position]

        num_aligned_pos = sum(column_encoding[:,0])
        if num_aligned_pos == 0:
            num_aligned_pos = 1
        column_encoding[:,0] /= num_aligned_pos
        column_encoding[:,1] /= AVG_COV
        column_encoding[:,2] /= (num_aligned_pos * AVG_QUA)


        num_aligned_neg = sum(column_encoding[:,3])
        if num_aligned_neg == 0:
            num_aligned_neg = 1
        column_encoding[:,3] /= num_aligned_neg
        column_encoding[:,4] /= AVG_COV
        column_encoding[:,5] /= (num_aligned_neg * AVG_QUA)

        # print("num_aligned_pos: ", num_aligned_pos)
        # print("num_aligned_neg: ", num_aligned_neg)
        # print(column_encoding)
        input_tensor.append(column_encoding)

    samfile.close()

    input_tensor = np.array(input_tensor)
    input_tensor = torch.Tensor(input_tensor)
    torch.save(input_tensor, "../Dataset/tranx_feature/" + gene_id + ".pt")
    print("Length: ", features[gene_id]['gene_end'] - features[gene_id]['gene_start'] + 1)
    print("input_tensor: ", input_tensor.size())
        #     f = open("../Dataset/tranx_label/"+file, 'r')
        #     line = f.relsadline()
        #     line = line.split("\t")
        #     region = line[0]+":"+line[1]+"-"+line[2]
        #     print(region)

        #     output_bam = "../Dataset/tranx_feature/" + file[4:-4]+ ".bam"
        #     subprocess.run(["samtools", "view", "../Dataset/SRR1352415.bam", region, "-o", output_bam ])

# samtools view input.bam "Chr10:18000-45500" > output.bam




            


if __name__ == "__main__":
    main(sys.argv[1:])