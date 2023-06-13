import os
import pickle 
import sys

# with open("./spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/d_scores_"+TYPE+".pkl", 'rb') as fw:
#     pickle.dump(donor_scores, fw)
# with open("./spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/a_scores_"+TYPE+".pkl", 'rb') as fw:
#     pickle.dump(acceptor_scores, fw)

SPLICEAI_VERSION = sys.argv[1]

output_files = ["pos_MANE", "pos_ALTS", "neg_1", "neg_random"]
for output_file in output_files:

    for TYPE in ['N', "noN"]:
        print(">> SpliceAI, TYPE: ", TYPE)
        # spliceai
        print(">> output_file\t: ", output_file)

        with open("./spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/d_scores_"+TYPE+".pkl", "rb") as fr:
            donor_scores = pickle.load(fr)
            # print("donor_scores: ", donor_scores)
            print("donor_scores: ", len(donor_scores))

        with open("./spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/a_scores_"+TYPE+".pkl", "rb") as fr:
            acceptor_scores = pickle.load(fr)
            # print("acceptor_scores: ", acceptor_scores)
            print("acceptor_scores: ", len(acceptor_scores))
    
    print("\n")

    for MODEL_VERSION in ['SPLAM_v11']:#, "SPLAM_v12"]:
        print(">> SPLAM, MODEL_VERSION: ", MODEL_VERSION)
        # splam
        TYPE = "noshuffle"
        MODEL_VERSION = "SPLAM_v11"
        # splam
        print(">> output_file\t: ", output_file)

        with open("./splam_result/"+MODEL_VERSION+"/"+output_file+"/d_scores_"+TYPE+".pkl", "rb") as fr:
            donor_scores = pickle.load(fr)
            # print("donor_scores: ", donor_scores)
            print("donor_scores: ", len(donor_scores))

        with open("./splam_result/"+MODEL_VERSION+"/"+output_file+"/a_scores_"+TYPE+".pkl", "rb") as fr:
            acceptor_scores = pickle.load(fr)
            # print("acceptor_scores: ", acceptor_scores)
            print("acceptor_scores: ", len(acceptor_scores))

    print("\n\n")



for TYPE in ['N', "noN"]:
    print(">> SpliceAI, TYPE: ", TYPE)
    # spliceai

    with open("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da."+TYPE+".merged.BOTH.pkl", "rb") as fr:
        d_label = pickle.load(fr)
        d_pred = pickle.load(fr)
        a_label = pickle.load(fr)
        a_pred = pickle.load(fr)

        # j_pred_prob = [x.numpy() for x in j_pred_prob]
        # j_pred_prob = [x.numpy() for x in j_pred_prob]

        print("\td_label : ", len(d_label))
        print("\td_pred: ", len(d_pred))
        print("")
        print("\ta_label : ", len(a_label))
        print("\ta_pred: ", len(a_pred))
        print("")

print("\n")

for MODEL_VERSION in ['SPLAM_v11', "SPLAM_v12"]:
    print(">> SPLAM, MODEL_VERSION: ", MODEL_VERSION)
    # splam
    TYPE = "noshuffle"
    MODEL_VERSION = "SPLAM_v11"
    # splam
    print(">> output_file\t: ", output_file)

    with open("./splam_result/"+MODEL_VERSION+"/splam.da."+TYPE+".merged.BOTH.pkl", "rb") as fr:
        d_label = pickle.load(fr)
        d_pred = pickle.load(fr)
        a_label = pickle.load(fr)
        a_pred = pickle.load(fr)

        # j_pred_prob = [x.numpy() for x in j_pred_prob]
        # j_pred_prob = [x.numpy() for x in j_pred_prob]

        print("\td_label : ", len(d_label))
        print("\td_pred: ", len(d_pred))
        print("")
        print("\ta_label : ", len(a_label))
        print("\ta_pred: ", len(a_pred))
        print("")

print("\n\n")
