import os
import pickle 
import sys
import numpy as np

# with open("./spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/d_scores_"+TYPE+".pkl", 'rb') as fw:
#     pickle.dump(donor_scores, fw)
# with open("./spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/a_scores_"+TYPE+".pkl", 'rb') as fw:
#     pickle.dump(acceptor_scores, fw)

# SPLICEAI_VERSION = sys.argv[1]
output_files = ["pos_MANE", "pos_ALTS", "BOTH", "FULL"]#, "neg_1", "neg_random"]
for output_file in output_files:

    for TYPE in ['N', "noN"]:

        d_pred_global = np.array([])
        # a_label_global = pickle.load(fr)
        a_pred_global = np.array([])

        for SPLICEAI_VERSION in ["1", "2", "3", "4", "5"]:
            print("@@ SPLICEAI_VERSION: ", SPLICEAI_VERSION)

            print(">> SpliceAI, TYPE: ", TYPE)
            # spliceai
            print(">> output_file\t: ", output_file)


            with open("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da."+TYPE+".merged."+output_file+".pkl", "rb") as fr:
                d_label = pickle.load(fr)
                d_pred = pickle.load(fr)
                a_label = pickle.load(fr)
                a_pred = pickle.load(fr)

            if (len(d_pred_global) == 0 or len(a_pred_global) == 0):
                d_pred_global = d_pred
                a_pred_global = a_pred    
            else:
                d_pred_global = d_pred_global + d_pred
                a_pred_global = a_pred_global + a_pred

            print("\t d_pred_global : ", d_pred_global)
            print("\td_pred: ", d_pred)
            print("")
            print("\ta_pred_global : ", a_pred_global)
            print("\ta_pred: ", a_pred)
            print("")

            print("\td_label : ", len(d_label))
            print("\td_pred: ", len(d_pred))
            print("")
            print("\ta_label : ", len(a_label))
            print("\ta_pred: ", len(a_pred))
            print("")
        d_pred_global /= 5
        a_pred_global /= 5

        os.makedirs("./spliceai_result_AVERAGE/", exist_ok=True)
        with open("./spliceai_result_AVERAGE/spliceai.da."+TYPE+".merged."+output_file+".pkl", 'wb') as f: 
            pickle.dump(d_label, f)
            pickle.dump(d_pred_global, f)
            pickle.dump(a_label, f)
            pickle.dump(a_pred_global, f)
        print("\n")


    # for MODEL_VERSION in ['SPLAM_v11']:#, "SPLAM_v12"]:
    #     print(">> SPLAM, MODEL_VERSION: ", MODEL_VERSION)
    #     # splam
    #     TYPE = "noshuffle"
    #     MODEL_VERSION = "SPLAM_v11"
    #     # splam
    #     print(">> output_file\t: ", output_file)

    #     with open("./splam_result/"+MODEL_VERSION+"/"+output_file+"/d_scores_"+TYPE+".pkl", "rb") as fr:
    #         donor_scores = pickle.load(fr)
    #         # print("donor_scores: ", donor_scores)
    #         print("donor_scores: ", len(donor_scores))

    #     with open("./splam_result/"+MODEL_VERSION+"/"+output_file+"/a_scores_"+TYPE+".pkl", "rb") as fr:
    #         acceptor_scores = pickle.load(fr)
    #         # print("acceptor_scores: ", acceptor_scores)
    #         print("acceptor_scores: ", len(acceptor_scores))
    # print("\n\n")



# # This is for both
# for TYPE in ['N', "noN"]:
#     print(">> SpliceAI, TYPE: ", TYPE)
#     # spliceai

#     # d_label_global = np.array()
#     d_pred_global = np.array([])
#     # a_label_global = pickle.load(fr)
#     a_pred_global = np.array([])

#     for SPLICEAI_VERSION in ["1", "2", "3", "4", "5"]:
#         print("@@ SPLICEAI_VERSION: ", SPLICEAI_VERSION)

#         with open("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da."+TYPE+".merged.BOTH.pkl", "rb") as fr:
#             d_label = pickle.load(fr)
#             d_pred = pickle.load(fr)
#             a_label = pickle.load(fr)
#             a_pred = pickle.load(fr)
#         if (len(d_pred_global) == 0 or len(a_pred_global) == 0):
#             d_pred_global = d_pred
#             a_pred_global = a_pred    
#         else:
#             d_pred_global = d_pred_global + d_pred
#             a_pred_global = a_pred_global + a_pred

#         # j_pred_prob = [x.numpy() for x in j_pred_prob]
#         # j_pred_prob = [x.numpy() for x in j_pred_prob]

#         print("\t d_pred_global : ", d_pred_global)
#         print("\td_pred: ", d_pred)
#         print("")
#         print("\ta_pred_global : ", a_pred_global)
#         print("\ta_pred: ", a_pred)
#         print("")

#         print("\td_label : ", len(d_label))
#         print("\td_pred: ", len(d_pred))
#         print("")
#         print("\ta_label : ", len(a_label))
#         print("\ta_pred: ", len(a_pred))
#         print("")
#     d_pred_global /= 5
#     a_pred_global /= 5
#     print("\t d_pred_global : ", d_pred_global)
#     print("\ta_pred_global : ", a_pred_global)
#     os.makedirs("./spliceai_result_AVERAGE/", exist_ok=True)
#     with open("./spliceai_result_AVERAGE/spliceai.da."+TYPE+".merged.BOTH.pkl", 'wb') as f: 
#         pickle.dump(d_label, f)
#         pickle.dump(d_pred_global, f)
#         pickle.dump(a_label, f)
#         pickle.dump(a_pred_global, f)
# print("\n")



# for MODEL_VERSION in ['SPLAM_v11', "SPLAM_v12"]:
#     print(">> SPLAM, MODEL_VERSION: ", MODEL_VERSION)
#     # splam
#     TYPE = "noshuffle"
#     MODEL_VERSION = "SPLAM_v11"
#     # splam
#     # print(">> output_file\t: ", output_file)

#     with open("./splam_result/"+MODEL_VERSION+"/splam.da."+TYPE+".merged.BOTH.pkl", "rb") as fr:
#         d_label = pickle.load(fr)
#         d_pred = pickle.load(fr)
#         a_label = pickle.load(fr)
#         a_pred = pickle.load(fr)

#         # j_pred_prob = [x.numpy() for x in j_pred_prob]
#         # j_pred_prob = [x.numpy() for x in j_pred_prob]

#         print("\td_label : ", len(d_label))
#         print("\td_pred: ", len(d_pred))
#         print("")
#         print("\ta_label : ", len(a_label))
#         print("\ta_pred: ", len(a_pred))
#         print("")

# print("\n\n")
