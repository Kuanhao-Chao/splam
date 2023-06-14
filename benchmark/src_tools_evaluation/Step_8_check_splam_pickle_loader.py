import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
import sys
# from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve

def main():


    # SUBSET = 9900
    # TARGETS = ["pos_MANE", "pos_ALTS", "neg_1", "neg_random"]
    
    for MANE_OR_ALTS in ["pos_MANE", "pos_ALTS", "BOTH", "FULL"]:
        TARGETS = [MANE_OR_ALTS, "neg_1", "neg_random"]
        SUBSETS = [2000, 10000, 10000]

        if MANE_OR_ALTS == "BOTH":
            TARGETS = ["pos_MANE", "pos_ALTS", "neg_1", "neg_random"]
            SUBSETS = [2000, 2000, 10000, 10000]

        if MANE_OR_ALTS == "FULL":
            TARGETS = ["pos_MANE", "pos_ALTS", "neg_1", "neg_random"]
            SUBSETS = [10000, 10000, 10000, 10000]

        for SPLAM_VERSION in ["SPLAM_v11"]:#, "SPLAM_v12"]:
            a_label = []
            d_label = []    
            a_pred = []
            d_pred = []     
            junc_name = []
            print(">> Processing ", SPLAM_VERSION)
            for idx in range(len(TARGETS)):
                TARGET = TARGETS[idx]
                SUBSET = SUBSETS[idx]

                ###########################
                # Donor scores
                ###########################
                print("./splam_result/"+SPLAM_VERSION+"/"+TARGET+"/d_scores_noshuffle.pkl")        
                with open("./splam_result/"+SPLAM_VERSION+"/"+TARGET+"/d_scores_noshuffle.pkl", 'rb') as f:
                    print(">> Processing ", TARGET)
                    splam_noS_d_pred_prob = pickle.load(f)[:SUBSET]
                    
                    if TARGET == "pos_MANE" or TARGET == "pos_ALTS":
                        splam_noS_d_label_prob = np.ones(SUBSET)
                    else:
                        splam_noS_d_label_prob = np.zeros(SUBSET)
                    print("len(splam_noS_d_label_prob): ", len(splam_noS_d_label_prob))
                    print("len(splam_noS_d_pred_prob): ", len(splam_noS_d_pred_prob))

                d_label = np.concatenate((d_label, splam_noS_d_label_prob), axis=None)        
                d_pred = np.concatenate((d_pred, splam_noS_d_pred_prob), axis=None)        


                ###########################
                # Acceptor scores
                ###########################
                print("./splam_result/"+SPLAM_VERSION+"/"+TARGET+"/a_scores_noshuffle.pkl")        
                with open("./splam_result/"+SPLAM_VERSION+"/"+TARGET+"/a_scores_noshuffle.pkl", 'rb') as f:
                    print(">> Processing ", TARGET)
                    splam_noS_a_pred_prob = pickle.load(f)[:SUBSET]
                    if TARGET == "pos_MANE" or TARGET == "pos_ALTS":
                        splam_noS_a_label_prob = np.ones(SUBSET)
                    else:
                        splam_noS_a_label_prob = np.zeros(SUBSET)
                    print("len(splam_noS_a_label_prob): ", len(splam_noS_a_label_prob))
                    print("len(splam_noS_a_pred_prob): ", len(splam_noS_a_pred_prob))

                a_label = np.concatenate((a_label, splam_noS_a_label_prob), axis=None)        
                a_pred = np.concatenate((a_pred, splam_noS_a_pred_prob), axis=None)        


                print("\td_pred : ", len(d_pred))
                print("\td_label: ", len(d_label))

                print("\td_pred: ", d_pred)
                print("\td_label: ", d_label)

                print("\ta_pred : ", len(a_pred))
                print("\ta_label: ", len(a_label))

                print("\ta_pred: ", a_pred)
                print("\ta_label: ", a_label)
                print("")


            ###################################
            # Writing results to pickle
            ###################################
            print(">> Final check before output.")
            print("\td_pred : ", len(d_pred))
            print("\td_label: ", len(d_label))

            print("\td_pred: ", d_pred)
            print("\td_pred > 0.5: ", len(d_pred[d_pred > 0.5]))
            print("\td_label: ", d_label)
            print("\td_label == 1: ", len(d_label[d_label == 1]))

            print("\ta_pred : ", len(a_pred))
            print("\ta_label: ", len(a_label))

            print("\ta_pred: ", a_pred)
            print("\ta_label: ", a_label)

            print("\ta_pred: ", a_pred)
            print("\ta_pred > 0.5: ", len(a_pred[a_pred > 0.5]))
            print("\ta_label: ", a_label)
            print("\ta_label == 1: ", len(a_label[a_label == 1]))
            print("")
            with open("./splam_result/"+SPLAM_VERSION+"/splam.da.noshuffle.merged."+MANE_OR_ALTS+".pkl", 'wb') as f: 
                pickle.dump(d_label, f)
                pickle.dump(d_pred, f)
                pickle.dump(a_label, f)
                pickle.dump(a_pred, f)

            ###################################
            # Checking splam pkl file.
            ###################################
            with open("./splam_result/"+SPLAM_VERSION+"/splam.da.noshuffle.merged."+MANE_OR_ALTS+".pkl",'rb') as f:
                print("./splam_result/"+SPLAM_VERSION+"/splam.da.noshuffle.merged."+MANE_OR_ALTS+".pkl")
                j_pred_prob_min = pickle.load(f)
                j_label_prob_min = pickle.load(f)
                j_pred_prob_avg = pickle.load(f)
                j_label_prob_avg = pickle.load(f)

                # j_pred_prob = [x.numpy() for x in j_pred_prob]
                # j_pred_prob = [x.numpy() for x in j_pred_prob]

                print("\tj_pred_prob : ", len(j_pred_prob_min))
                print("\tj_label_prob: ", len(j_label_prob_min))
                print("")
                print("\tj_pred_prob : ", len(j_pred_prob_avg))
                print("\tj_label_prob: ", len(j_label_prob_avg))
                print("")

if __name__ == "__main__":
    main()