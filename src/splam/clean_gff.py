import os
import json
import gffutils

def clean_gff(outdir, gff_db, threshold, bad_intron_num):
    threshold = float(threshold)
    bad_intron_num = int(bad_intron_num)
    print(f"[Info] Threshold: '{threshold}'")
    print(f"[Info] Bad intron number threshold: '{bad_intron_num}'\n")
    trans_intron_num_txt = os.path.join(outdir, "intron_num.txt")

    
    # reading the data from the file
    with open(trans_intron_num_txt) as f:
        intron_2_num = f.read()        
    # reconstructing the data as a dictionary
    intron_2_num = json.loads(intron_2_num)

    TOTAL_FILE = 50
    junc_score = outdir + "/junction_score.bed"
    with open(junc_score) as fr:
        lines = fr.read().splitlines()
        trans_dict = {}
        trans_ratio_dict = {}
        for line in lines:
            chrs, start, end, name, aln_num, strand, d_score, a_score, trans = line.split('\t')
            trans = trans.split(',')
            if float(d_score) < threshold and float(a_score) < threshold:
                for tran in trans:
                    if tran in trans_dict.keys():
                        trans_dict[tran] += 1
                    else:
                        trans_dict[tran] = 1
            else:
                for tran in trans:
                    if tran in trans_dict.keys():
                        pass
                    else:
                        trans_dict[tran] = 0
        transcript_num = 0
        ls = [0]*50
        trans_dir = outdir + "/trans/"
        # fw_list = []  # Empty list to store file contents
        os.makedirs(trans_dir, exist_ok=True)



        # for i in range(0, TOTAL_FILE):
        #     file_name = f"{trans_dir}trans_{i}.txt"  # Assuming the files are named as file_1.txt, file_2.txt, and so on
        #     try:
        #         fw = open(file_name, 'w')
        #         fw_list.append(fw)
        #     except FileNotFoundError:
        #         print(f"[ERROR] File {file_name} not found.")

        for tran, count in trans_dict.items():


            trans_ratio_dict[tran] = count / intron_2_num[tran]

            # if count > (TOTAL_FILE-1):
            #     count = (TOTAL_FILE-1)
            # fw_list[count].write(f'{tran}\n')

            # ls[count] += 1  
            transcript_num += 1

        # for id in range(len(ls)):
        #     print(f'[Info] {id} bad introns: {ls[id]} transcripts')
        #     fw_list[id].close()

        # print("trans_ratio_dict: ", trans_ratio_dict)
        print("[Info] Total transcripts: ", transcript_num)



        # Load the GFF database
        db = gffutils.FeatureDB(gff_db, keep_order=True)

        ################################################
        # Create a list of transcripts that should be removed
        ################################################
        # trans_removed = []
        # for i in range(bad_intron_num, 50):
        #     file_name = f"{trans_dir}trans_{i}.txt"  # Assuming the files are named as file_1.txt, file_2.txt, and so on
        #     try:
        #         fr = open(file_name, 'r')
        #         content = fr.read().splitlines()
        #         trans_removed = trans_removed + content
        #     except FileNotFoundError:
        #         print(f"File {file_name} not found.")

        trans_removed = []
        for tran, ratio in trans_ratio_dict.items():
            if ratio > 0.9 and intron_2_num[tran] >= 3:
                trans_removed.append(tran)

        print("trans_removed: ", len(trans_removed))
        print(f'[Info] Number of removed transcripts: {len(trans_removed)}\n')


        output_file = outdir + '/cleaned.gff'
        # Open the output file in write mode


        genes = db.features_of_type("gene")

        kept_transcript_count = 0
        with open(output_file, 'w') as fw:
            for gene in genes:
                # Retrieve child features
                children = db.children(gene, level=1)
                # Print child features
                for child in children:
                    kept_transcript_count += 1
                    if not child.id in trans_removed:
                        fw.write(str(child) + '\n')
                        exons = db.children(child, featuretype='exon', order_by='start', level=1)
                        for exon in exons:    
                            fw.write(str(exon) + '\n')

        print(f'[Info] {kept_transcript_count} transcripts processed\n')
        print(f"[Info] GFF file written to '{output_file}'")
              