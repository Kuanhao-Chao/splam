import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
from util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, PrecisionRecallDisplay

def plot_pr_curve(true_y, y_prob, label, option):
    """
    plots the roc curve based of the probabilities
    """
    # precision, recall, thresholds = precision_recall_curve_self(true_y, y_prob, label, option)
    precision, recall, thresholds = precision_recall_curve(true_y, y_prob)
    plt.plot(recall, precision)
    # , label=label, marker='o')
    plt.legend()
    plt.xlabel('Recall')
    plt.ylabel('Precision')

def plot_pr_curve_J(true_y, y_prob_d, y_prob_a, label, option, choice):
    """
    plots the roc curve based of the probabilities
    """
    # precision, recall, thresholds = precision_recall_curve_J_level_self(true_y, y_prob_d, y_prob_a, label, option)
    precision, recall, thresholds = precision_recall_curve_J_level_self(true_y, y_prob_d, y_prob_a, label, option, choice)
    plt.plot(recall, precision, label=label)
    # , marker='o')
    plt.legend()
    plt.xlabel('Recall')
    plt.ylabel('Precision')

def plot_roc_curve(true_y, y_prob, label, option):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = roc_curve(true_y, y_prob)
    # fpr, tpr, thresholds = ROC_curve_self(true_y, y_prob, label, option)
    plt.plot(fpr, tpr, label=label)
    # , marker='o' )
    plt.legend()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

def plot_roc_curve_J(true_y, y_prob_d, y_prob_a, label, option, choice):
    """
    plots the roc curve based of the probabilities
    """
    precision, recall, thresholds = roc_curve_J_level_self(true_y, y_prob_d, y_prob_a, label, option, choice)
    plt.plot(precision, recall, label=label)
    # , marker='o')
    plt.legend()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')


def plot_thresholds(true_y, y_prob_d, y_prob_a, type):
    # 0.1 - 0.9 (9)
    D_TP = [0]*1001
    D_TN = [0]*1001
    D_FP = [0]*1001
    D_FN = [0]*1001

    A_TP = [0]*1001
    A_TN = [0]*1001
    A_FP = [0]*1001
    A_FN = [0]*1001

    J_TP = [0]*1001
    J_TN = [0]*1001
    J_FP = [0]*1001
    J_FN = [0]*1001
    """
    plots TP, FP, TN, FN
    """
    thresholds = []
    for idx in range(0, 1001, 1):
        threshold = (idx)/1000
        thresholds.append(threshold)
        print("threshold: ", threshold)
        # print("Len(true_y)  : ", len(true_y))
        # print("Len(y_prob_d): ", len(y_prob_d))
        # print("Len(y_prob_a): ", len(y_prob_a))
        # print("Len(true_y): ", true_y)

        labels_1 = np.where(true_y == 1)
        labels_0 = np.where(true_y == 0)
        ####################
        # Donor 
        ####################
        thre_d = np.where(y_prob_d >= threshold)
        thre_d_0 = np.where(y_prob_d < threshold)
        # print("thre_d: ", thre_d)
        # print("thre_a: ", thre_a)

        TPs = len(np.intersect1d(labels_1, thre_d))
        FNs = len(np.setdiff1d(labels_1, thre_d))
        FPs = len(np.setdiff1d(thre_d, labels_1))
        TNs = len(np.intersect1d(labels_0, thre_d_0))
        # TNs = len(true_y) - TPs - FNs - FPs
        print("\tDonor TPs: ", TPs)
        print("\tDonor FNs: ", FNs)
        print("\tDonor FPs: ", FPs)
        print("\tDonor TNs: ", TNs)
        D_TP[idx] = TPs
        D_TN[idx] = TNs
        D_FP[idx] = FPs
        D_FN[idx] = FNs

        ####################
        # Acceptor 
        ####################
        thre_a = np.where(y_prob_a >= threshold)
        thre_a_0 = np.where(y_prob_a < threshold)
        TPs = len(np.intersect1d(labels_1, thre_a))
        FNs = len(np.setdiff1d(labels_1, thre_a))
        FPs = len(np.setdiff1d(thre_a, labels_1))
        TNs = len(np.intersect1d(labels_0, thre_a_0))
        # TNs = len(true_y) - TPs - FNs - FPs
        print("\tAcceptor TPs: ", TPs)
        print("\tAcceptor FNs: ", FNs)
        print("\tAcceptor FPs: ", FPs)
        print("\tAcceptor TNs: ", TNs)
        A_TP[idx] = TPs
        A_TN[idx] = TNs
        A_FP[idx] = FPs
        A_FN[idx] = FNs

        ####################
        # junction 
        ####################
        thre_j = np.intersect1d(thre_d, thre_a)
        TPs = len(np.intersect1d(labels_1, thre_j))
        FNs = len(np.setdiff1d(labels_1, thre_j))
        FPs = len(np.setdiff1d(thre_j, labels_1))
        TNs = len(true_y) - TPs - FNs - FPs
        print("\tJunction TPs: ", TPs)
        print("\tJunction FNs: ", FNs)
        print("\tJunction FPs: ", FPs)
        print("\tJunction TNs: ", TNs)
        J_TP[idx] = TPs
        J_TN[idx] = TNs
        J_FP[idx] = FPs
        J_FN[idx] = FNs

    ####################
    # Donor 
    ####################
    D_TP = np.array(D_TP)
    D_TN = np.array(D_TN)
    D_FP = np.array(D_FP)
    D_FN = np.array(D_FN)

    plt.bar(thresholds, D_TP, color='r', width=0.001)
    plt.bar(thresholds, D_FN, bottom=D_TP, width=0.001, color='b')
    plt.bar(thresholds, D_TN, bottom=D_TP+D_FN, width=0.001, color='y')
    plt.bar(thresholds, D_FP, bottom=D_TP+D_FN+D_TN, width=0.001, color='g')

    # # , color='b')
    # plt.bar(thresholds, D_FP, bottom=D_TP+D_TN)
    # # , color='b')
    # plt.bar(thresholds, D_FN, bottom=D_TP+D_TN+D_FP)
    # , color='b')
    plt.xlabel("Thresholds")
    plt.ylabel("Count")
    plt.legend(["TP", "FN", "TN", "FP"])
    # plt.title("Scores by Teams in 4 Rounds")
    plt.savefig("./IMG/T_histogram/"+type+"_donor_thresholds.png")
    plt.close()

    ####################
    # Acceptor 
    ####################
    A_TP = np.array(A_TP)
    A_TN = np.array(A_TN)
    A_FP = np.array(A_FP)
    A_FN = np.array(A_FN)

    plt.bar(thresholds, A_TP, color='r', width=0.001)
    plt.bar(thresholds, A_FN, bottom=A_TP, width=0.001, color='b')
    plt.bar(thresholds, A_TN, bottom=A_TP+A_FN, width=0.001, color='y')
    plt.bar(thresholds, A_FP, bottom=A_TP+A_FN+A_TN, width=0.001, color='g')

    # # , color='b')
    # plt.bar(thresholds, D_FP, bottom=D_TP+D_TN)
    # # , color='b')
    # plt.bar(thresholds, D_FN, bottom=D_TP+D_TN+D_FP)
    # , color='b')
    plt.xlabel("Thresholds")
    plt.ylabel("Count")
    plt.legend(["TP", "FN", "TN", "FP"])
    # plt.title("Scores by Teams in 4 Rounds")
    plt.savefig("./IMG/T_histogram/"+type+"_acceptor_thresholds.png")
    plt.close()

    ####################
    # junction 
    ####################
    J_TP = np.array(J_TP)
    J_TN = np.array(J_TN)
    J_FP = np.array(J_FP)
    J_FN = np.array(J_FN)

    plt.bar(thresholds, J_TP, color='r', width=0.001)
    plt.bar(thresholds, J_FN, bottom=J_TP, width=0.001, color='b')
    plt.bar(thresholds, J_TN, bottom=J_TP+J_FN, width=0.001, color='y')
    plt.bar(thresholds, J_FP, bottom=J_TP+J_FN+J_TN, width=0.001, color='g')

    plt.xlabel("Thresholds")
    plt.ylabel("Count")
    plt.legend(["TP", "FN", "TN", "FP"])
    # plt.title("Scores by Teams in 4 Rounds")
    plt.savefig("./IMG/T_histogram/"+type+"_junction_thresholds.png")
    plt.close()
        

def plot_thresholds_J(true_y, y_prob_j, type):
    J_TP = [0]*1001
    J_TN = [0]*1001
    J_FP = [0]*1001
    J_FN = [0]*1001
    """
    plots TP, FP, TN, FN
    """
    thresholds = []
    for idx in range(0, 1001, 1):
        threshold = (idx)/1000
        thresholds.append(threshold)
        print("threshold: ", threshold)
        
        labels_1 = np.where(true_y == 1)
        labels_0 = np.where(true_y == 0)

        ####################
        # Junction
        ####################
        thre_j = np.where(y_prob_j >= threshold)
        thre_j_0 = np.where(y_prob_j < threshold)
        # print("thre_d: ", thre_d)
        # print("thre_a: ", thre_a)

        TPs = len(np.intersect1d(labels_1, thre_j))
        FNs = len(np.setdiff1d(labels_1, thre_j))
        FPs = len(np.setdiff1d(thre_j, labels_1))
        TNs = len(np.intersect1d(labels_0, thre_j_0))
        # TNs = len(true_y) - TPs - FNs - FPs
        print("\tDonor TPs: ", TPs)
        print("\tDonor FNs: ", FNs)
        print("\tDonor FPs: ", FPs)
        print("\tDonor TNs: ", TNs)
        J_TP[idx] = TPs
        J_TN[idx] = TNs
        J_FP[idx] = FPs
        J_FN[idx] = FNs

    ####################
    # junction 
    ####################
    J_TP = np.array(J_TP)
    J_TN = np.array(J_TN)
    J_FP = np.array(J_FP)
    J_FN = np.array(J_FN)

    plt.bar(thresholds, J_TP, color='r', width=0.001)
    plt.bar(thresholds, J_FN, bottom=J_TP, width=0.001, color='b')
    plt.bar(thresholds, J_TN, bottom=J_TP+J_FN, width=0.001, color='y')
    plt.bar(thresholds, J_FP, bottom=J_TP+J_FN+J_TN, width=0.001, color='g')

    plt.xlabel("Thresholds")
    plt.ylabel("Count")
    plt.legend(["TP", "FN", "TN", "FP"])
    # plt.title("Scores by Teams in 4 Rounds")
    plt.savefig("./IMG/T_histogram/"+type+"_junction_thresholds.png")
    plt.close()


def main():
    #####################################
    # Creating directories for visualization.
    #####################################
    os.makedirs("./IMG/donor/", exist_ok=True)
    os.makedirs("./IMG/acceptor/", exist_ok=True)
    os.makedirs("./IMG/junction/", exist_ok=True)
    os.makedirs("./IMG/T_histogram/", exist_ok=True)

    #####################################
    # Declaring parameters for probability & prediction array
    #####################################
    spliceai_N_d_pred_prob = []
    spliceai_N_d_label_prob = []
    spliceai_N_a_pred_prob = []
    spliceai_N_a_label_prob = []

    spliceai_noN_d_pred_prob = []
    spliceai_noN_d_label_prob = []
    spliceai_noN_a_pred_prob = []
    spliceai_noN_a_label_prob = []

    splam_j_shuffle_pred_prob_min = []
    splam_j_shuffle_label_prob_min = []
    splam_j_noshuffle_pred_prob_min = []
    splam_j_noshuffle_label_prob_min = []

    splam_j_shuffle_pred_prob_avg = []
    splam_j_shuffle_label_prob_avg = []
    splam_j_noshuffle_pred_prob_avg = []
    splam_j_noshuffle_label_prob_avg = []


    # splam_j_nobatch_pred_prob = []
    # splam_j_nobatch_label_prob = []


    with open("./spliceai_result/spliceai.N.merged.pkl",'rb') as f:
        spliceai_N_d_pred_prob = pickle.load(f)
        spliceai_N_d_label_prob = pickle.load(f)
        spliceai_N_a_pred_prob = pickle.load(f)
        spliceai_N_a_label_prob = pickle.load(f)

        # spliceai_N_d_pred_prob = [x.numpy() for x in spliceai_N_d_pred_prob]
        # spliceai_N_a_pred_prob = [x.numpy() for x in spliceai_N_a_pred_prob]
        spliceai_N_d_pred_prob = np.array(spliceai_N_d_pred_prob)
        spliceai_N_a_pred_prob = np.array(spliceai_N_a_pred_prob)
        spliceai_N_d_label_prob = np.array(spliceai_N_d_label_prob)
        spliceai_N_a_label_prob = np.array(spliceai_N_a_label_prob)
        # spliceai_N_d_label_prob[:1000] = 1
        # spliceai_N_a_label_prob[:1000] = 1

        # spliceai_noN_d_label_prob = [float(i) for i in spliceai_noN_d_label_prob]
        # spliceai_noN_a_label_prob = [float(i) for i in spliceai_noN_a_label_prob]

        print("spliceai_N_d_pred_prob : ", spliceai_N_d_pred_prob)
        print("spliceai_N_d_label_prob: ", spliceai_N_d_label_prob)
        print("spliceai_N_a_pred_prob : ", spliceai_N_a_pred_prob)
        print("spliceai_N_a_label_prob: ", spliceai_N_a_label_prob)

        print("spliceai_N_d_pred_prob : ", len(spliceai_N_d_pred_prob))
        print("spliceai_N_d_label_prob: ", len(spliceai_N_d_label_prob))
        print("spliceai_N_a_pred_prob : ", len(spliceai_N_a_pred_prob))
        print("spliceai_N_a_label_prob: ", len(spliceai_N_a_label_prob))

    with open("./spliceai_result/spliceai.noN.merged.pkl",'rb') as f:
        spliceai_noN_d_pred_prob = pickle.load(f)
        spliceai_noN_d_label_prob = pickle.load(f)
        spliceai_noN_a_pred_prob = pickle.load(f)
        spliceai_noN_a_label_prob = pickle.load(f)

        # spliceai_noN_d_pred_prob = [x.numpy() for x in spliceai_noN_d_pred_prob]
        # spliceai_noN_a_pred_prob = [x.numpy() for x in spliceai_noN_a_pred_prob]
        spliceai_noN_d_pred_prob = np.array(spliceai_noN_d_pred_prob)
        spliceai_noN_a_pred_prob = np.array(spliceai_noN_a_pred_prob)
        spliceai_noN_d_label_prob = np.array(spliceai_noN_d_label_prob)
        spliceai_noN_a_label_prob = np.array(spliceai_noN_a_label_prob)
        # spliceai_noN_d_label_prob[:1000] = 1
        # spliceai_noN_a_label_prob[:1000] = 1


        # spliceai_noN_d_label_prob = [float(i) for i in spliceai_noN_d_label_prob]
        # spliceai_noN_a_label_prob = [float(i) for i in spliceai_noN_a_label_prob]

        print("spliceai_noN_d_pred_prob : ", spliceai_noN_d_pred_prob)
        print("spliceai_noN_d_label_prob: ", spliceai_noN_d_label_prob)
        print("spliceai_noN_a_pred_prob : ", spliceai_noN_a_pred_prob)
        print("spliceai_noN_a_label_prob: ", spliceai_noN_a_label_prob)

        print("spliceai_noN_d_pred_prob : ", len(spliceai_noN_d_pred_prob))
        print("spliceai_noN_d_label_prob: ", len(spliceai_noN_d_label_prob))
        print("spliceai_noN_a_pred_prob : ", len(spliceai_noN_a_pred_prob))
        print("spliceai_noN_a_label_prob: ", len(spliceai_noN_a_label_prob))

    with open("./splam_result/splam.shuffle.pkl",'rb') as f:
        splam_j_shuffle_pred_prob_min = pickle.load(f)
        splam_j_shuffle_label_prob_min = pickle.load(f)
        
        splam_j_shuffle_pred_prob_avg = pickle.load(f)
        splam_j_shuffle_label_prob_avg = pickle.load(f)

        print("splam_j_shuffle_pred_prob_min : ", splam_j_shuffle_pred_prob_min)
        print("splam_j_shuffle_label_prob_min: ", splam_j_shuffle_label_prob_min)

        print("splam_j_shuffle_pred_prob_min : ", len(splam_j_shuffle_pred_prob_min))
        print("splam_j_shuffle_label_prob_min: ", len(splam_j_shuffle_label_prob_min))

    with open("./splam_result/splam.noshuffle.pkl",'rb') as f:
        splam_j_noshuffle_pred_prob_min = pickle.load(f)
        splam_j_noshuffle_label_prob_min = pickle.load(f)

        splam_j_noshuffle_pred_prob_avg = pickle.load(f)
        splam_j_noshuffle_label_prob_avg = pickle.load(f)

        print("splam_j_noshuffle_pred_prob_min : ", splam_j_noshuffle_pred_prob_min)
        print("splam_j_noshuffle_label_prob_min: ", splam_j_noshuffle_label_prob_min)

        print("splam_j_noshuffle_pred_prob_min : ", len(splam_j_noshuffle_pred_prob_min))
        print("splam_j_noshuffle_label_prob_min: ", len(splam_j_noshuffle_label_prob_min))

    # with open("./spliceai_result/splam.nobatch.0.100.pkl",'rb') as f:
    #     splam_j_nobatch_pred_prob = pickle.load(f)
    #     splam_j_nobatch_label_prob = pickle.load(f)
    #     print("splam_j_nobatch_pred_prob : ", splam_j_nobatch_pred_prob)
    #     print("splam_j_nobatch_label_prob: ", splam_j_nobatch_label_prob)
    #     print("splam_j_nobatch_pred_prob : ", len(splam_j_nobatch_pred_prob))
    #     print("splam_j_nobatch_label_prob: ", len(splam_j_nobatch_label_prob))

    splam_j_shuffle_label_prob_min = np.array(splam_j_shuffle_label_prob_min)
    splam_j_shuffle_pred_prob_min = np.array(splam_j_shuffle_pred_prob_min)
    splam_j_noshuffle_label_prob_min = np.array(splam_j_noshuffle_label_prob_min)
    splam_j_noshuffle_pred_prob_min = np.array(splam_j_noshuffle_pred_prob_min)

    splam_j_shuffle_label_prob_avg = np.array(splam_j_shuffle_label_prob_avg)
    splam_j_shuffle_pred_prob_avg = np.array(splam_j_shuffle_pred_prob_avg)
    splam_j_noshuffle_label_prob_avg = np.array(splam_j_noshuffle_label_prob_avg)
    splam_j_noshuffle_pred_prob_avg = np.array(splam_j_noshuffle_pred_prob_avg)

    # splam_j_nobatch_label_prob = np.array(splam_j_nobatch_label_prob)
    # splam_j_nobatch_pred_prob = np.array(splam_j_nobatch_pred_prob)

    print("splam_j_noshuffle_label_prob_min: ", splam_j_noshuffle_label_prob_min)
    print("splam_j_noshuffle_pred_prob_min: ", splam_j_noshuffle_pred_prob_min)
    

    
    ###################################
    # Threshold plots
    ###################################
    plot_thresholds(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, spliceai_noN_a_pred_prob, "splcieai")
    plot_thresholds(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_pred_prob, "splcieai_N")
    plot_thresholds_J(splam_j_shuffle_label_prob_min, splam_j_shuffle_pred_prob_min, "splam_shuffle")
    plot_thresholds_J(splam_j_noshuffle_label_prob_min, splam_j_noshuffle_pred_prob_min, "splam_noshuffle")

    plot_thresholds_J(splam_j_shuffle_label_prob_avg, splam_j_shuffle_pred_prob_avg, "splam_shuffle_avg")
    plot_thresholds_J(splam_j_noshuffle_label_prob_avg, splam_j_noshuffle_pred_prob_avg, "splam_noshuffle_avg")
    # plot_thresholds_J(splam_j_nobatch_label_prob, splam_j_nobatch_pred_prob, "splam_nobatch")

    
    ###################################
    # PR of junction (min)
    ###################################
    plot_pr_curve_J(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, spliceai_noN_a_pred_prob, "spliceai_junc_sklearn", "sklearn", "min")
    # plot_pr_curve_J(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, spliceai_noN_a_pred_prob, "spliceai_junc_self", "self")

    plot_pr_curve_J(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_pred_prob, "spliceai_N_junc_sklearn", "sklearn", "min")
    # plot_pr_curve_J(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_pred_prob, "spliceai_N_junc_self", "self")

    plot_pr_curve(splam_j_shuffle_label_prob_min, splam_j_shuffle_pred_prob_min, "splam_junc_shuffle", "sklean")
    plot_pr_curve(splam_j_noshuffle_label_prob_min, splam_j_noshuffle_pred_prob_min, "splam_junc_noshuffle", "sklean")
    # plot_pr_curve(splam_j_nobatch_label_prob, splam_j_nobatch_pred_prob, "splam_junc_nobatch", "sklean")
    plt.savefig("./IMG/junction/junc_pr_min.png")
    plt.close()

    ###################################
    # PR of junction (avg)
    ###################################
    plot_pr_curve_J(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, spliceai_noN_a_pred_prob, "spliceai_junc_sklearn", "sklearn", "avg")
    # plot_pr_curve_J(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, spliceai_noN_a_pred_prob, "spliceai_junc_self", "self")

    plot_pr_curve_J(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_pred_prob, "spliceai_N_junc_sklearn", "sklearn", "avg")
    # plot_pr_curve_J(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_pred_prob, "spliceai_N_junc_self", "self")

    plot_pr_curve(splam_j_shuffle_label_prob_avg, splam_j_shuffle_pred_prob_avg, "splam_junc_shuffle", "sklean")
    plot_pr_curve(splam_j_noshuffle_label_prob_avg, splam_j_noshuffle_pred_prob_avg, "splam_junc_noshuffle", "sklean")
    # plot_pr_curve(splam_j_nobatch_label_prob, splam_j_nobatch_pred_prob, "splam_junc_nobatch", "sklean")
    plt.savefig("./IMG/junction/junc_pr_avg.png")
    plt.close()
    


    ###################################
    # ROC of junction
    ###################################
    plot_roc_curve_J(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, spliceai_noN_a_pred_prob, "spliceai_junc_sklearn", "sklearn", "min")
    # plot_roc_curve_J(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, spliceai_noN_a_pred_prob, "spliceai_junc_self", "self")

    plot_roc_curve_J(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_pred_prob, "spliceai_N_junc_sklearn", "sklearn", "min")
    # plot_roc_curve_J(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_pred_prob, "spliceai_N_junc_self", "self")

    plot_roc_curve(splam_j_shuffle_label_prob_min, splam_j_shuffle_pred_prob_min, "splam_junc_shuffle", "self")
    plot_roc_curve(splam_j_noshuffle_label_prob_min, splam_j_noshuffle_pred_prob_min, "splam_junc_noshuffle", "self")
    # plot_roc_curve(splam_j_nobatch_label_prob, splam_j_nobatch_pred_prob, "splam_junc_nobatch", "self")
    plt.savefig("./IMG/junction/junc_roc_min.png")
    plt.close()
    # print("spliceai_noN_d_label_prob: ", spliceai_noN_d_label_prob)
    # print("spliceai_noN_d_pred_prob: ", spliceai_noN_d_pred_prob)


    ################################### 
    # ROC of junction
    ###################################
    plot_roc_curve_J(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, spliceai_noN_a_pred_prob, "spliceai_junc_sklearn", "sklearn", "avg")
    # plot_roc_curve_J(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, spliceai_noN_a_pred_prob, "spliceai_junc_self", "self")

    plot_roc_curve_J(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_pred_prob, "spliceai_N_junc_sklearn", "sklearn", "avg")
    # plot_roc_curve_J(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_pred_prob, "spliceai_N_junc_self", "self")

    plot_roc_curve(splam_j_shuffle_label_prob_avg, splam_j_shuffle_pred_prob_avg, "splam_junc_shuffle", "self")
    plot_roc_curve(splam_j_noshuffle_label_prob_avg, splam_j_noshuffle_pred_prob_avg, "splam_junc_noshuffle", "self")
    # plot_roc_curve(splam_j_nobatch_label_prob, splam_j_nobatch_pred_prob, "splam_junc_nobatch", "self")
    plt.savefig("./IMG/junction/junc_roc_avg.png")
    plt.close()
    # print("spliceai_noN_d_label_prob: ", spliceai_noN_d_label_prob)
    # print("spliceai_noN_d_pred_prob: ", spliceai_noN_d_pred_prob)


    # ###################################
    # # PR / ROC of donor
    # ###################################
    # plot_pr_curve(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, "spliceai_donor", "self")
    # plot_pr_curve(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, "spliceai_N_donor", "self")
    # plot_pr_curve(splam_d_label_prob, splam_d_pred_prob, "splam_donor", "self")
    # # plt.show()
    # plt.savefig("./IMG/donor/donor_pr.png")
    # plt.close()

    # plot_roc_curve(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, "spliceai_donor", "self")
    # plot_roc_curve(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, "spliceai_N_donor", "self")
    # plot_roc_curve(splam_d_label_prob, splam_d_pred_prob, "splam_donor", "self")
    # # plt.show()
    # plt.savefig("./IMG/donor/donor_roc.png")
    # plt.close()

    # ###################################
    # # PR / ROC of acceptor
    # ###################################
    # plot_pr_curve(spliceai_noN_a_label_prob, spliceai_noN_a_pred_prob, "spliceai_acceptor", "self")
    # plot_pr_curve(spliceai_N_a_label_prob, spliceai_N_a_pred_prob, "spliceai_N_acceptor", "self")
    # plot_pr_curve(splam_a_label_prob, splam_a_pred_prob, "splam_acceptor", "self")
    # # plt.show()
    # plt.savefig("./IMG/acceptor/acceptor_pr.png")
    # plt.close()

    # plot_roc_curve(spliceai_noN_a_label_prob, spliceai_noN_a_pred_prob, "spliceai_acceptor", "self")
    # plot_roc_curve(spliceai_N_a_label_prob, spliceai_N_a_pred_prob, "spliceai_N_acceptor", "self")
    # plot_roc_curve(splam_a_label_prob, splam_a_pred_prob, "splam_acceptor", "self")
    # # plt.show()
    # plt.savefig("./IMG/acceptor/acceptor_roc.png")
    # plt.close()


if __name__ == "__main__":
    main()