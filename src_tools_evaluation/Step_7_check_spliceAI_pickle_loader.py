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

    d_pred = []
    d_label = []
    a_pred = []
    a_label = []
    junc_name = []

    splam_j_pred_prob = []
    splam_j_label_prob = []


    SUBSET = 10000

    for TYPE in ["N", "noN"]:
    # for TYPE in ["noN"]:
        ###################################
        # Checking spliceai pkl file (pos)
        ###################################
        # with open("./spliceai_result_pos_only/spliceai."+TYPE+".pkl",'rb') as f:
        #     print("Results of './spliceai_result_pos_only/spliceai."+TYPE+".pkl'")
        #     d_pred_pos_prob = pickle.load(f)
        #     d_label_pos_prob = pickle.load(f)
        #     d_pred_pos_prob = d_pred_pos_prob[:SUBSET]
        #     d_label_pos_prob = np.array(d_label_pos_prob)
        #     d_label_pos_prob[:SUBSET] = 1
        #     d_label_pos_prob = d_label_pos_prob[:SUBSET]

        #     a_pred_pos_prob = pickle.load(f)
        #     a_label_pos_prob = pickle.load(f)
        #     a_pred_pos_prob = a_pred_pos_prob[:SUBSET]
        #     a_label_pos_prob = np.array(a_label_pos_prob)
        #     a_label_pos_prob[:SUBSET] = 1
        #     a_label_pos_prob = a_label_pos_prob[:SUBSET]


        #     print("\td_pred_pos_prob : ", len(d_pred_pos_prob))
        #     print("\td_label_pos_prob: ", len(d_label_pos_prob))
        #     # print("\td_pred_pos_prob: ", d_pred_pos_prob)
        #     # print("\td_label_pos_prob: ", d_label_pos_prob)
        #     print("\ta_pred_pos_prob : ", len(a_pred_pos_prob))
        #     print("\ta_label_pos_prob: ", len(a_label_pos_prob))
        #     # print("\ta_pred_pos_prob: ", a_pred_pos_prob)
        #     # print("\ta_label_pos_prob: ", a_label_pos_prob)
        #     # print("")
        
        ###################################
        # Checking spliceai pkl file (pos)
        ###################################
        with open("./spliceai_result/spliceai."+TYPE+".pos.pkl",'rb') as f:
            print("Results of './spliceai_result/spliceai."+TYPE+".pos.pkl'")
            d_pred_pos = pickle.load(f)
            d_label_pos = pickle.load(f)
            d_label_pos = np.array(d_label_pos)
            d_pred_pos = d_pred_pos[:SUBSET]
            d_label_pos = d_label_pos[:SUBSET]
            # d_label_neg[:1000] = 1
            a_pred_pos = pickle.load(f)
            a_label_pos = pickle.load(f)
            a_label_pos = np.array(a_label_pos)
            a_pred_pos = a_pred_pos[:SUBSET]
            a_label_pos = a_label_pos[:SUBSET]
            junc_name_pos= pickle.load(f)

            # print("\td_pred_pos : ", len(d_pred_pos))
            # print("\td_label_pos: ", len(d_label_pos))
            # # print("\td_pred_neg: ", d_pred_neg)
            # # print("\td_label_neg: ", d_label_neg)
            # print("\ta_pred_pos : ", len(a_pred_pos))
            # print("\ta_label_pos: ", len(a_label_pos))
            # print("\tjunc_name_pos: ", len(junc_name_pos))
            # print("\tjunc_name_pos: ", junc_name_pos)
            # print("\ta_pred_neg: ", a_pred_neg)
            # print("\ta_label_neg: ", a_label_neg)
            # print("")

        ###################################
        # Checking spliceai pkl file (neg)
        ###################################
        with open("./spliceai_result/spliceai."+TYPE+".neg_1.pkl",'rb') as f:
            print("Results of './spliceai_result/spliceai."+TYPE+".pkl'")
            d_pred_neg = pickle.load(f)
            d_label_neg = pickle.load(f)
            d_label_neg = np.array(d_label_neg)
            d_pred_neg = d_pred_neg[:SUBSET]
            d_label_neg = d_label_neg[:SUBSET]
            # d_label_neg[:1000] = 1
            a_pred_neg = pickle.load(f)
            a_label_neg = pickle.load(f)
            a_label_neg = np.array(a_label_neg)
            a_pred_neg = a_pred_neg[:SUBSET]
            a_label_neg = a_label_neg[:SUBSET]
            junc_name_neg = pickle.load(f)
            # a_label_neg[:1000] = 1


            # print("\td_pred_neg : ", len(d_pred_neg))
            # print("\td_label_neg: ", len(d_label_neg))
            # # print("\td_pred_neg: ", d_pred_neg)
            # # print("\td_label_neg: ", d_label_neg)
            # print("\ta_pred_neg : ", len(a_pred_neg))
            # print("\ta_label_neg: ", len(a_label_neg))
            # print("\tjunc_name_neg: ", len(junc_name_neg))
            # print("\tjunc_name_neg: ", junc_name_neg)
            # print("\ta_pred_neg: ", a_pred_neg)
            # print("\ta_label_neg: ", a_label_neg)
            # print("")

        a_label = np.concatenate((a_label_pos, a_label_neg), axis=None)        
        d_label = np.concatenate((d_label_pos, d_label_neg), axis=None)        

        a_pred = np.concatenate((a_pred_pos, a_pred_neg), axis=None)        
        d_pred = np.concatenate((d_pred_pos, d_pred_neg), axis=None)        

        junc_name = (junc_name_pos + junc_name_neg)

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