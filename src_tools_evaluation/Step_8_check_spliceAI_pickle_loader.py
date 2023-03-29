import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
# from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve

def main():

    # os.makedirs("./IMG/spliceai/", exist_ok=True)
    # os.makedirs("./IMG/splam/", exist_ok=True)
    
    d_pred_neg = []
    d_label_neg = []
    a_pred_neg = []
    a_label_neg = []
    junc_name_neg = []

    d_pred_pos = []
    d_label_pos = []
    a_pred_pos = []
    a_label_pos = []
    junc_name_pos = []

    d_pred_pos_refseq_protein_isoforms = []
    d_label_pos_refseq_protein_isoforms = []
    a_pred_pos_refseq_protein_isoforms = []
    a_label_pos_refseq_protein_isoforms = []
    junc_name_pos_refseq_protein_isoforms = []

    d_pred = []
    d_label = []
    a_pred = []
    a_label = []
    junc_name = []

    splam_j_pred_prob = []
    splam_j_label_prob = []


    # SUBSET = 9900
    TARGETS = ["pos_refseq_protein_alts", "neg_1"]
    # , "neg_random"]
    # TARGETS = ["pos", "neg_1"]
    # TARGETS = ["pos", "neg_1", "neg_5"]
    # SUBSETS = [2000, 5000, 5000]
    # SUBSETS = [9900, 9900]
    # SUBSETS = [500, 8000]

    SUBSETS = [1000, 8000]
            #    , 9000]

    for TYPE in ["N", "noN"]:
    # for TYPE in ["noN"]:
        a_label = []
        d_label = []    
        a_pred = []
        d_pred = []     
        junc_name = []

        for idx in range(len(TARGETS)):
            target = TARGETS[idx]
            SUBSET = SUBSETS[idx]
            ###################################
            # Checking spliceai pkl file (pos)
            ###################################
            with open("./spliceai_result/spliceai."+TYPE+"."+target+".pkl",'rb') as f:
                print("Results of './spliceai_result/spliceai."+TYPE+"."+target+".pkl'")
                d_pred_lcl = pickle.load(f)
                d_label_lcl = pickle.load(f)
                d_label_lcl = np.array(d_label_lcl)
                d_pred_lcl = d_pred_lcl[:SUBSET]
                d_label_lcl = d_label_lcl[:SUBSET]

                # d_label_neg[:1000] = 1
                a_pred_lcl = pickle.load(f)
                a_label_lcl = pickle.load(f)
                a_label_lcl = np.array(a_label_lcl)
                a_pred_lcl = a_pred_lcl[:SUBSET]
                a_label_lcl = a_label_lcl[:SUBSET]
                junc_name_lcl = pickle.load(f)
                junc_name_lcl = junc_name_lcl[:SUBSET]
                
                if target == "pos" or target == "pos_refseq_protein_all" or target == "pos_refseq_protein_alts":
                    d_label_lcl = np.ones(SUBSET)
                    a_label_lcl = np.ones(SUBSET)

                # print("\td_pred_pos : ", len(d_pred_pos))
                print("\td_label_lcl: ", len(d_label_lcl))
                # print("\td_pred_neg: ", d_pred_neg)
                print("\td_label_lcl: ", d_label_lcl)
                # print("\ta_pred_pos : ", len(a_pred_pos))
                print("\ta_label_lcl: ", len(a_label_lcl))
                print("\ta_label_lcl: ", a_label_lcl)
                # print("\tjunc_name_pos: ", len(junc_name_pos))
                # print("\tjunc_name_pos: ", junc_name_pos)
                # print("\ta_pred_neg: ", a_pred_neg)
                # print("\ta_label_neg: ", a_label_neg)
                # print("")

            a_label = np.concatenate((a_label, a_label_lcl), axis=None)        
            d_label = np.concatenate((d_label, d_label_lcl), axis=None)        

            a_pred = np.concatenate((a_pred, a_pred_lcl), axis=None)        
            d_pred = np.concatenate((d_pred, d_pred_lcl), axis=None)        

            junc_name = (junc_name + junc_name_lcl)

            print("\td_pred : ", len(d_pred))
            print("\td_label: ", len(d_label))

            print("\td_pred: ", d_pred)
            print("\td_label: ", d_label)

            print("\ta_pred : ", len(a_pred))
            print("\ta_label: ", len(a_label))

            print("\ta_pred: ", a_pred)
            print("\ta_label: ", a_label)
            print("\n\tjunc_name: ", len(junc_name))
            # print("\tjunc_name: ", junc_name)
            print("")


        ###################################
        # Writing results to pickle
        ###################################
        with open("./spliceai_result/spliceai."+TYPE+".merged.pkl", 'wb') as f: 
            pickle.dump(d_pred, f)
            pickle.dump(d_label, f)
            pickle.dump(a_pred, f)
            pickle.dump(a_label, f)
            pickle.dump(junc_name, f)



    # ###################################
    # # Checking splam pkl file.
    # ###################################
    # for TYPE in ["shuffle", "noshuffle"]:
    #     with open("./splam_result/splam."+TYPE+".pkl",'rb') as f:
    #         print("Results of './splam_result/splam."+TYPE+".pkl'")
    #         j_pred_prob_min = pickle.load(f)
    #         j_label_prob_min = pickle.load(f)
    #         j_pred_prob_avg = pickle.load(f)
    #         j_label_prob_avg = pickle.load(f)

    #         # j_pred_prob = [x.numpy() for x in j_pred_prob]
    #         # j_pred_prob = [x.numpy() for x in j_pred_prob]

    #         print("\tj_pred_prob : ", len(j_pred_prob_min))
    #         print("\tj_label_prob: ", len(j_label_prob_min))
    #         print("")
    #         print("\tj_pred_prob : ", len(j_pred_prob_avg))
    #         print("\tj_label_prob: ", len(j_label_prob_avg))
    #         print("")

if __name__ == "__main__":
    main()