import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
import sys
# from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve

def main():

    SPLICEAI_VERSION = sys.argv[1]
    
    for MANE_OR_ALTS in ["pos_MANE", "pos_ALTS", "BOTH", "FULL"]:
        TARGETS = [MANE_OR_ALTS, "neg_1", "neg_random"]
        SUBSETS = [2000, 10000, 10000]

        if MANE_OR_ALTS == "BOTH":
            TARGETS = ["pos_MANE", "pos_ALTS", "neg_1", "neg_random"]
            SUBSETS = [2000, 2000, 10000, 10000]

        if MANE_OR_ALTS == "FULL":
            TARGETS = ["pos_MANE", "pos_ALTS", "neg_1", "neg_random"]
            SUBSETS = [10000, 10000, 10000, 10000]

        for TYPE in ["noN", "N"]:
            a_label = []
            d_label = []    
            a_pred = []
            d_pred = []     
            junc_name = []
            for idx in range(len(TARGETS)):
                TARGET = TARGETS[idx]
                SUBSET = SUBSETS[idx]

                ###########################
                # Donor scores
                ###########################
                print("./spliceai_result_"+SPLICEAI_VERSION+"/"+TARGET+"/d_scores_"+TYPE+".pkl")        
                with open("./spliceai_result_"+SPLICEAI_VERSION+"/"+TARGET+"/d_scores_"+TYPE+".pkl", 'rb') as f:
                    print(">> Processing ", TARGET)
                    spliceai_noS_d_pred_prob = pickle.load(f)[:SUBSET]
                    if TARGET == "pos_MANE" or TARGET == "pos_ALTS":
                        spliceai_noS_d_label_prob = np.ones(SUBSET)
                    else:
                        spliceai_noS_d_label_prob = np.zeros(SUBSET)                    

                    print("len(spliceai_noS_d_label_prob): ", len(spliceai_noS_d_label_prob))
                    print("len(spliceai_noS_d_pred_prob): ", len(spliceai_noS_d_pred_prob))

                d_label = np.concatenate((d_label, spliceai_noS_d_label_prob), axis=None)
                d_pred = np.concatenate((d_pred, spliceai_noS_d_pred_prob), axis=None)        


                ###########################
                # Acceptor scores
                ###########################
                print("./spliceai_result_"+SPLICEAI_VERSION+"/"+TARGET+"/a_scores_"+TYPE+".pkl")        
                with open("./spliceai_result_"+SPLICEAI_VERSION+"/"+TARGET+"/a_scores_"+TYPE+".pkl", 'rb') as f:
                    print(">> Processing ", TARGET)
                    spliceai_noS_a_pred_prob = pickle.load(f)[:SUBSET]
                    if TARGET == "pos_MANE" or TARGET == "pos_ALTS":
                        spliceai_noS_a_label_prob = np.ones(SUBSET)
                    else:
                        spliceai_noS_a_label_prob = np.zeros(SUBSET)                    
                    print("len(spliceai_noS_a_label_prob): ", len(spliceai_noS_a_label_prob))
                    print("len(spliceai_noS_a_pred_prob): ", len(spliceai_noS_a_pred_prob))

                a_label = np.concatenate((a_label, spliceai_noS_a_label_prob), axis=None)
                a_pred = np.concatenate((a_pred, spliceai_noS_a_pred_prob), axis=None)        


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
            with open("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da."+TYPE+".merged."+MANE_OR_ALTS+".pkl", 'wb') as f: 
                pickle.dump(d_label, f)
                pickle.dump(d_pred, f)
                pickle.dump(a_label, f)
                pickle.dump(a_pred, f)

            ###################################
            # Checking spliceai pkl file.
            ###################################
            with open("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da."+TYPE+".merged."+MANE_OR_ALTS+".pkl",'rb') as f:
                print("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da."+TYPE+".merged."+MANE_OR_ALTS+".pkl")
                d_label = pickle.load(f)
                d_pred = pickle.load(f)
                a_label = pickle.load(f)
                a_pred = pickle.load(f)

                # j_pred_prob = [x.numpy() for x in j_pred_prob]
                # j_pred_prob = [x.numpy() for x in j_pred_prob]

                print("\td_label : ", len(d_label))
                print("\td_pred: ", len(d_pred))
                print("")
                print("\ta_label : ", len(a_label))
                print("\ta_pred: ", len(a_pred))
                print("")

if __name__ == "__main__":
    main()