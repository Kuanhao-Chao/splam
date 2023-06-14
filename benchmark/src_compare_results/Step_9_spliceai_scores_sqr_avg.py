import os
import pickle 
import sys
import numpy as np

# with open("./spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/d_scores_"+TYPE+".pkl", 'rb') as fw:
#     pickle.dump(donor_scores, fw)
# with open("./spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/a_scores_"+TYPE+".pkl", 'rb') as fw:
#     pickle.dump(acceptor_scores, fw)

# SPLICEAI_VERSION = sys.argv[1]
output_files = ["pos_MANE", "pos_ALTS"]#, "neg_1", "neg_random"]
for output_file in output_files:

    for TYPE in ['N', "noN"]:

        d_pred_global = np.array([])
        # a_label_global = pickle.load(fr)
        a_pred_global = np.array([])

        for SPLICEAI_VERSION in range(1,6, 1):
            for SPLICEAI_VERSION_2 in range(SPLICEAI_VERSION+1, 6, 1):
                print("@@ SPLICEAI_VERSION: ", SPLICEAI_VERSION)
                print("@@ SPLICEAI_VERSION: ", SPLICEAI_VERSION_2)
                SPLICEAI_VERSION = str(SPLICEAI_VERSION)
                SPLICEAI_VERSION_2 = str(SPLICEAI_VERSION_2)

                if (SPLICEAI_VERSION == SPLICEAI_VERSION_2):
                    continue

                print(">> SpliceAI, TYPE: ", TYPE)
                # spliceai
                print(">> output_file\t: ", output_file)


                with open("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da."+TYPE+".merged."+output_file+".pkl", "rb") as fr:
                    d_label_1 = pickle.load(fr)
                    d_pred_1 = pickle.load(fr)
                    a_label_1 = pickle.load(fr)
                    a_pred_1 = pickle.load(fr)

                with open("./spliceai_result_"+SPLICEAI_VERSION_2+"/spliceai.da."+TYPE+".merged."+output_file+".pkl", "rb") as fr:
                    d_label_2 = pickle.load(fr)
                    d_pred_2 = pickle.load(fr)
                    a_label_2 = pickle.load(fr)
                    a_pred_2 = pickle.load(fr)
                
                d_pred_diff = np.square(d_pred_1 - d_pred_2)
                a_pred_diff = np.square(a_pred_1 - a_pred_2)


                print("\t >>> d_pred_global : ", d_pred_global)
                print("\t >>> a_pred_global : ", a_pred_global)
                

                if (len(d_pred_global) == 0 or len(a_pred_global) == 0):
                    d_pred_global = d_pred_diff
                    a_pred_global = a_pred_diff    
                else:
                    d_pred_global = d_pred_global + d_pred_diff
                    a_pred_global = a_pred_global + a_pred_diff

                print("\t d_pred: ", d_pred_diff)
                print("\t a_pred: ", a_pred_diff)

                print("\t >>> after d_pred_global : ", d_pred_global)
                print("\t >>> after a_pred_global : ", a_pred_global)
                print("\n\n")

                # print("\td_label : ", len(d_label))
                # print("\td_pred: ", len(d_pred))
                # print("")
                # print("\ta_label : ", len(a_label))
                # print("\ta_pred: ", len(a_pred))
                # print("")

        d_pred_global /= 10
        a_pred_global /= 10

        print("## MAX d_pred_global: ", max(d_pred_global))
        print("## MAX a_pred_global: ", max(a_pred_global))

        os.makedirs("./spliceai_result_AVERAGE/", exist_ok=True)
        with open("./spliceai_result_AVERAGE/spliceai.da."+TYPE+".merged."+output_file+".distance.pkl", 'wb') as f: 
            pickle.dump(d_label_1, f)
            pickle.dump(d_pred_global, f)
            pickle.dump(d_label_1, f)
            pickle.dump(a_pred_global, f)
        print("\n")