import os
import gffutils

def clean_gff(outdir, gff_db, threshold, bad_intron_num):
    threshold = float(threshold)
    bad_intron_num = int(bad_intron_num)
    print(f"[Info] Threshold: '{threshold}'")
    print(f"[Info] Bad intron number threshold: '{bad_intron_num}'\n")

    TOTAL_FILE = 50
    junc_score = outdir + "/junction_score.bed"
    with open(junc_score) as fr:
        lines = fr.read().splitlines()
        trans_dict = {}

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
        fw_list = []  # Empty list to store file contents
        os.makedirs(trans_dir, exist_ok=True)
        for i in range(0, TOTAL_FILE):
            file_name = f"{trans_dir}trans_{i}.txt"  # Assuming the files are named as file_1.txt, file_2.txt, and so on
            try:
                fw = open(file_name, 'w')
                fw_list.append(fw)
            except FileNotFoundError:
                print(f"[ERROR] File {file_name} not found.")

        for tran, count in trans_dict.items():
            if count > (TOTAL_FILE-1):
                count = (TOTAL_FILE-1)
            fw_list[count].write(f'{tran}\n')

            ls[count] += 1  
            transcript_num += 1

        for id in range(len(ls)):
            print(f'[Info] {id} bad introns: {ls[id]} transcripts')
            fw_list[id].close()

        print("[Info] Total transcripts: ", transcript_num)


        # Load the GFF database
        db = gffutils.FeatureDB(gff_db, keep_order=True)

        ################################################
        # Create a list of transcripts that should be removed
        ################################################
        trans_removed = []
        for i in range(bad_intron_num, 50):
            file_name = f"{trans_dir}trans_{i}.txt"  # Assuming the files are named as file_1.txt, file_2.txt, and so on
            try:
                fr = open(file_name, 'r')
                content = fr.read().splitlines()
                trans_removed = trans_removed + content
            except FileNotFoundError:
                print(f"File {file_name} not found.")
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
                

        
        
        
        
        
        
        
        
        
        
        
        # feature_counter = 0
        # with open(output_file, 'w') as f:
        #     # Iterate over all features in the database
        #     for feature in db.all_features():
        #         feature_counter += 1
        #         # Skip itself and its children.                
        #         print("feature.id: ", feature.id)
        #         if feature.id in trans_removed or (not feature.attributes.get('Parent') is None and feature.attributes.get('Parent') and feature.attributes['Parent'][0] in trans_removed):
        #             # print("feature.attributes.get('Parent'): ", feature.attributes.get('Parent')[0])
        #             continue  # Skip this feature and its children
                    
        #         # Write the feature as a GFF formatted line
        #         f.write(str(feature) + '\n')
        #         if feature_counter % 100000 == 0:
        #             print(f'[Info] {feature_counter} features processed')
        # print(f'[Info] {feature_counter} features processed\n')
        # print(f"[Info] GFF file written to '{output_file}'")