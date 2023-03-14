import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
# from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve

def main():

    SUBSET = 1000
    boundary = 1000
    # for TYPE in ["N", "noN"]:
    # with open("./splam_result/splam.da.shuffle.pkl",'rb') as f:
    #     splam_S_d_label_prob = pickle.load(f)
    #     splam_S_d_pred_prob = pickle.load(f)
    #     splam_S_a_label_prob = pickle.load(f)
    #     splam_S_a_pred_prob = pickle.load(f)
    #     splam_S_junc_name = pickle.load(f)


    #     splam_S_d_label_prob = np.concatenate((splam_S_d_label_prob[:SUBSET], splam_S_d_label_prob[boundary:boundary+SUBSET]), axis=None)
    #     splam_S_d_pred_prob = np.concatenate((splam_S_d_pred_prob[:SUBSET], splam_S_d_pred_prob[boundary:boundary+SUBSET]), axis=None)
    #     splam_S_a_label_prob = np.concatenate((splam_S_a_label_prob[:SUBSET], splam_S_a_label_prob[boundary:boundary+SUBSET]), axis=None)
    #     splam_S_a_pred_prob = np.concatenate((splam_S_a_pred_prob[:SUBSET], splam_S_a_pred_prob[boundary:boundary+SUBSET]), axis=None)
    #     splam_S_junc_name = splam_S_junc_name[:SUBSET] + splam_S_junc_name[boundary:boundary+SUBSET]

    #     print("splam_S_d_label_prob : ", splam_S_d_label_prob)
    #     print("splam_S_d_pred_prob  : ", splam_S_d_pred_prob)

    #     print("splam_S_a_label_prob : ", splam_S_a_label_prob)
    #     print("splam_S_a_pred_prob  : ", splam_S_a_pred_prob)
    #     print("splam_S_junc_name    : ", splam_S_junc_name)

    
    with open("./splam_result/splam.da.noshuffle.pkl",'rb') as f:
        splam_noS_d_label_prob = pickle.load(f)
        splam_noS_d_pred_prob = pickle.load(f)
        splam_noS_a_label_prob = pickle.load(f)
        splam_noS_a_pred_prob = pickle.load(f)
        splam_noS_junc_name = pickle.load(f)

        splam_noS_d_label_prob = np.concatenate((splam_noS_d_label_prob[:SUBSET], splam_noS_d_label_prob[boundary:boundary+SUBSET]), axis=None)
        splam_noS_d_pred_prob = np.concatenate((splam_noS_d_pred_prob[:SUBSET], splam_noS_d_pred_prob[boundary:boundary+SUBSET]), axis=None)
        splam_noS_a_label_prob = np.concatenate((splam_noS_a_label_prob[:SUBSET], splam_noS_a_label_prob[boundary:boundary+SUBSET]), axis=None)
        splam_noS_a_pred_prob = np.concatenate((splam_noS_a_pred_prob[:SUBSET], splam_noS_a_pred_prob[boundary:boundary+SUBSET]), axis=None)
        splam_noS_junc_name = splam_noS_junc_name[:SUBSET] + splam_noS_junc_name[boundary:boundary+SUBSET]

        # print("splam_noS_d_label_prob : ", splam_noS_d_label_prob)
        # print("splam_noS_d_pred_prob  : ", splam_noS_d_pred_prob)

        # print("splam_noS_a_label_prob : ", splam_noS_a_label_prob)
        # print("splam_noS_a_pred_prob  : ", splam_noS_a_pred_prob)
        # print("splam_noS_junc_name    : ", splam_noS_junc_name)

        print("splam_noS_d_label_prob : ", len(splam_noS_d_label_prob))
        print("splam_noS_d_pred_prob  : ", len(splam_noS_d_pred_prob))

        print("splam_noS_a_label_prob : ", len(splam_noS_a_label_prob))
        print("splam_noS_a_pred_prob  : ", len(splam_noS_a_pred_prob))
        print("splam_noS_junc_name    : ", len(splam_noS_junc_name))
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