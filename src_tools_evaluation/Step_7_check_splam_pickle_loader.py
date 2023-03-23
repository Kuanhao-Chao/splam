import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
# from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve

def main():


    # SUBSET = 9900
    TARGETS = ["pos_refseq_protein_alts", "neg_1"]
    # , "neg_random"]
    # TARGETS = ["pos", "neg_1"]
    # TARGETS = ["pos", "neg_1", "neg_5"]
    # SUBSETS = [2000, 5000, 5000]
    # SUBSETS = [9900, 9900]
    # SUBSETS = [500, 8000]

    SUBSETS = [9000, 9000]
            #    , 9000]

    a_label = []
    d_label = []    
    a_pred = []
    d_pred = []     
    junc_name = []

    for idx in range(len(TARGETS)):
        target = TARGETS[idx]
        SUBSET = SUBSETS[idx]
        with open("./splam_result/splam.da.noshuffle."+target+".pkl",'rb') as f:
            print(">> Processing ", target)
            splam_noS_d_label_prob = pickle.load(f)[:SUBSET]
            splam_noS_d_pred_prob = pickle.load(f)[:SUBSET]
            splam_noS_a_label_prob = pickle.load(f)[:SUBSET]
            splam_noS_a_pred_prob = pickle.load(f)[:SUBSET]
            splam_noS_junc_name = pickle.load(f)[:SUBSET]

            # splam_noS_d_label_prob = np.concatenate((splam_noS_d_label_prob[:SUBSET], splam_noS_d_label_prob[:SUBSET]), axis=None)
            # splam_noS_d_pred_prob = np.concatenate((splam_noS_d_pred_prob[:SUBSET], splam_noS_d_pred_prob[:SUBSET]), axis=None)
            # splam_noS_a_label_prob = np.concatenate((splam_noS_a_label_prob[:SUBSET], splam_noS_a_label_prob[:SUBSET]), axis=None)
            # splam_noS_a_pred_prob = np.concatenate((splam_noS_a_pred_prob[:SUBSET], splam_noS_a_pred_prob[:SUBSET]), axis=None)
            # splam_noS_junc_name = splam_noS_junc_name[:SUBSET] + splam_noS_junc_name[:SUBSET]

            # print("splam_noS_d_label_prob : ", splam_noS_d_label_prob)
            # print("splam_noS_d_pred_prob  : ", splam_noS_d_pred_prob)

            # print("splam_noS_a_label_prob : ", splam_noS_a_label_prob)
            # print("splam_noS_a_pred_prob  : ", splam_noS_a_pred_prob)
            # print("splam_noS_junc_name    : ", splam_noS_junc_name)

            print("\tsplam_noS_d_label_prob : ", len(splam_noS_d_label_prob))
            print("\tsplam_noS_d_label_prob : ", splam_noS_d_label_prob)
            print("\tsplam_noS_d_pred_prob  : ", len(splam_noS_d_pred_prob))

            print("\tsplam_noS_a_label_prob : ", len(splam_noS_a_label_prob))
            print("\tsplam_noS_a_pred_prob  : ", len(splam_noS_a_pred_prob))
            print("\tsplam_noS_junc_name    : ", len(splam_noS_junc_name))

        a_label = np.concatenate((a_label, splam_noS_a_label_prob), axis=None)        
        d_label = np.concatenate((d_label, splam_noS_d_label_prob), axis=None)        
        a_pred = np.concatenate((a_pred, splam_noS_a_pred_prob), axis=None)        
        d_pred = np.concatenate((d_pred, splam_noS_d_pred_prob), axis=None)        

        junc_name = (junc_name + splam_noS_junc_name)

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
        # print("\tjunc_name: ", junc_name)
        print("")


        ###################################
        # Writing results to pickle
        ###################################
        with open("./splam_result/splam.da.noshuffle.merged.pkl", 'wb') as f: 
            pickle.dump(d_label, f)
            pickle.dump(d_pred, f)
            pickle.dump(a_label, f)
            pickle.dump(a_pred, f)
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