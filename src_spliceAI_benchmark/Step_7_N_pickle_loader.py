import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, PrecisionRecallDisplay

def plot_pr_curve(true_y, y_prob, label):
    """
    plots the roc curve based of the probabilities
    """
    precision, recall, thresholds = precision_recall_curve_self(true_y, y_prob, label)
    # for threshold in thresholds:
    #     print("threshold: ", threshold)
    # prd = PrecisionRecallDisplay(precision, recall)
    # prd.plot()
    plt.plot(recall, precision , label=label)
    print("precision: ", precision)
    print("recall   : ", recall)
    plt.legend()
    plt.xlabel('Recall')
    plt.ylabel('Precision')

def plot_pr_curve_J(true_y, y_prob_d, y_prob_a, label):
    """
    plots the roc curve based of the probabilities
    """
    precision, recall, thresholds = precision_recall_curve_J_level_self(true_y, y_prob_d, y_prob_a, label)
    # for threshold in thresholds:
    #     print("threshold: ", threshold)
    plt.plot(recall, precision, label=label)
    print("precision: ", precision)
    print("recall   : ", recall)
    plt.legend()
    plt.xlabel('Recall')
    plt.ylabel('Precision')

def plot_roc_curve(true_y, y_prob, label):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = roc_curve(true_y, y_prob)
    plt.plot(fpr, tpr, label=label)
    plt.legend()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

def plot_roc_curve_J(true_y, y_prob_d, y_prob_a, label):
    """
    plots the roc curve based of the probabilities
    """
    precision, recall, thresholds = roc_curve_J_level_self(true_y, y_prob_d, y_prob_a, label)
    # for threshold in thresholds:
    #     print("threshold: ", threshold)
    plt.plot(precision, recall, label=label)
    print("precision: ", precision)
    print("recall   : ", recall)
    plt.legend()
    plt.xlabel('Recall')
    plt.ylabel('Precision')


def plot_thresholds(true_y, y_prob_d, y_prob_a, type):

    # 0.1 - 0.9 (9)
    D_TP = [0]*999
    D_TN = [0]*999
    D_FP = [0]*999
    D_FN = [0]*999

    A_TP = [0]*999
    A_TN = [0]*999
    A_FP = [0]*999
    A_FN = [0]*999

    J_TP = [0]*999
    J_TN = [0]*999
    J_FP = [0]*999
    J_FN = [0]*999
    """
    plots TP, FP, TN, FN
    """
    thresholds = []
    for idx in range(0, 999, 1):
        threshold = (idx+1)/1000
        thresholds.append(threshold)
        print("threshold: ", threshold)
        # print("Len(true_y)  : ", len(true_y))
        # print("Len(y_prob_d): ", len(y_prob_d))
        # print("Len(y_prob_a): ", len(y_prob_a))
        # print("Len(true_y): ", true_y)

        labels_1 = np.where(true_y == 1)

        ####################
        # Donor 
        ####################
        thre_d = np.where(y_prob_d >= threshold)
        # print("thre_d: ", thre_d)
        # print("thre_a: ", thre_a)

        TPs = len(np.intersect1d(labels_1, thre_d))
        FNs = len(np.setdiff1d(labels_1, thre_d))
        FPs = len(np.setdiff1d(thre_d, labels_1))
        TNs = len(true_y) - TPs - FNs - FPs
        print("Donor TPs: ", TPs)
        print("Donor FNs: ", FNs)
        print("Donor FPs: ", FPs)
        print("Donor TNs: ", TNs)
        D_TP[idx] = TPs
        D_TN[idx] = TNs
        D_FP[idx] = FPs
        D_FN[idx] = FNs

        ####################
        # Acceptor 
        ####################
        thre_a = np.where(y_prob_a >= threshold)
        TPs = len(np.intersect1d(labels_1, thre_a))
        FNs = len(np.setdiff1d(labels_1, thre_a))
        FPs = len(np.setdiff1d(thre_a, labels_1))
        TNs = len(true_y) - TPs - FNs - FPs
        print("Acceptor TPs: ", TPs)
        print("Acceptor FNs: ", FNs)
        print("Acceptor FPs: ", FPs)
        print("Acceptor TNs: ", TNs)
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
        print("Junction TPs: ", TPs)
        print("Junction FNs: ", FNs)
        print("Junction FPs: ", FPs)
        print("Junction TNs: ", TNs)
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

    plt.bar(thresholds, D_TP, color='r', width=0.01)
    plt.bar(thresholds, D_FN, bottom=D_TP, width=0.01, color='b')
    plt.bar(thresholds, D_TN, bottom=D_TP+D_FN, width=0.01, color='y')
    plt.bar(thresholds, D_FP, bottom=D_TP+D_FN+D_TN, width=0.01, color='g')

    # # , color='b')
    # plt.bar(thresholds, D_FP, bottom=D_TP+D_TN)
    # # , color='b')
    # plt.bar(thresholds, D_FN, bottom=D_TP+D_TN+D_FP)
    # , color='b')
    plt.xlabel("Thresholds")
    plt.ylabel("Count")
    plt.legend(["TP", "FN", "TN", "FP"])
    # plt.title("Scores by Teams in 4 Rounds")
    plt.savefig("./IMG/"+type+"_donor_thresholds.png")
    plt.close()



    ####################
    # Acceptor 
    ####################
    A_TP = np.array(A_TP)
    A_TN = np.array(A_TN)
    A_FP = np.array(A_FP)
    A_FN = np.array(A_FN)

    plt.bar(thresholds, A_TP, color='r', width=0.01)
    plt.bar(thresholds, A_FN, bottom=A_TP, width=0.01, color='b')
    plt.bar(thresholds, A_TN, bottom=A_TP+A_FN, width=0.01, color='y')
    plt.bar(thresholds, A_FP, bottom=A_TP+A_FN+A_TN, width=0.01, color='g')

    # # , color='b')
    # plt.bar(thresholds, D_FP, bottom=D_TP+D_TN)
    # # , color='b')
    # plt.bar(thresholds, D_FN, bottom=D_TP+D_TN+D_FP)
    # , color='b')
    plt.xlabel("Thresholds")
    plt.ylabel("Count")
    plt.legend(["TP", "FN", "TN", "FP"])
    # plt.title("Scores by Teams in 4 Rounds")
    plt.savefig("./IMG/"+type+"_acceptor_thresholds.png")
    plt.close()


    ####################
    # junction 
    ####################
    J_TP = np.array(J_TP)
    J_TN = np.array(J_TN)
    J_FP = np.array(J_FP)
    J_FN = np.array(J_FN)

    plt.bar(thresholds, J_TP, color='r', width=0.01)
    plt.bar(thresholds, J_FN, bottom=J_TP, width=0.01, color='b')
    plt.bar(thresholds, J_TN, bottom=J_TP+J_FN, width=0.01, color='y')
    plt.bar(thresholds, J_FP, bottom=J_TP+J_FN+J_TN, width=0.01, color='g')

    plt.xlabel("Thresholds")
    plt.ylabel("Count")
    plt.legend(["TP", "FN", "TN", "FP"])
    # plt.title("Scores by Teams in 4 Rounds")
    plt.savefig("./IMG/"+type+"_junction_thresholds.png")
    plt.close()
        


        # print("thre_array: ", thre_array)

        # true_y
        # if true_y[idx] == 1 and acceptor_l[idx] == 1:
        #     # print("INNNN!!POS_NUM")
        #     if donor_p > threshold:
        #         D_TP[idx]+=1
        #     else:
        #         D_FN[idx]+=1

        #     if acceptor_p > threshold:
        #         A_TP[idx]+=1
        #     else:
        #         A_FN[idx]+=1

        #     if donor_p > threshold and acceptor_p > threshold:
        #         J_TP[idx]+=1
        #     else:
        #         J_FN[idx]+=1
        # else:
        #     # print("OUTOUT!!POS_NUM")
        #     if donor_p > threshold:
        #         D_FP[idx]+=1
        #     else:
        #         D_TN[idx]+=1

        #     if acceptor_p > threshold:
        #         A_FP[idx]+=1
        #     else:
        #         A_TN[idx]+=1

        #     if donor_p > threshold and acceptor_p > threshold:
        #         J_FP[idx]+=1
        #     else:
        #         J_TN[idx]+=1
    # plt.plot(precision, recall, label=label)
    # print("precision: ", precision)
    # print("recall   : ", recall)
    # plt.legend()
    # plt.xlabel('Recall')
    # plt.ylabel('Precision')


# def threshold_assessment(donor_p, acceptor_p, donor_l, acceptor_l):
#     for idx in range(0,9,1):
#         threshold = (idx+1)*0.1
#         print("threshold:  ", threshold)
#         # print("Donor   : ", np.where(Y_pred[:, :, 2] > threshold))
#         # print("Acceptor: ", np.where(Y_pred[:, :, 1] > threshold))

#         if donor_l[idx] == 1 and acceptor_l[idx] == 1:
#             # print("INNNN!!POS_NUM")
#             if donor_p > threshold:
#                 D_TP[idx]+=1
#             else:
#                 D_FN[idx]+=1

#             if acceptor_p > threshold:
#                 A_TP[idx]+=1
#             else:
#                 A_FN[idx]+=1

#             if donor_p > threshold and acceptor_p > threshold:
#                 J_TP[idx]+=1
#             else:
#                 J_FN[idx]+=1
#         else:
#             # print("OUTOUT!!POS_NUM")
#             if donor_p > threshold:
#                 D_FP[idx]+=1
#             else:
#                 D_TN[idx]+=1

#             if acceptor_p > threshold:
#                 A_FP[idx]+=1
#             else:
#                 A_TN[idx]+=1

#             if donor_p > threshold and acceptor_p > threshold:
#                 J_FP[idx]+=1
#             else:
#                 J_TN[idx]+=1


def main():

    ratio = 12000
    os.makedirs("./IMG/spliceai/", exist_ok=True)
    os.makedirs("./IMG/splam/", exist_ok=True)
    

    spliceai_N_d_pred_prob = []
    spliceai_N_d_label_prob = []
    spliceai_N_a_pred_prob = []
    spliceai_N_a_label_prob = []

    spliceai_d_pred_prob = []
    spliceai_d_label_prob = []
    spliceai_a_pred_prob = []
    spliceai_a_label_prob = []

    splam_d_pred_prob = []
    splam_d_label_prob = []
    splam_a_pred_prob = []
    splam_a_label_prob = []
    with open("./INPUT_N/spliceai_11812.pkl",'rb') as f:
        spliceai_N_d_pred_prob = pickle.load(f)
        spliceai_N_d_label_prob = pickle.load(f)
        spliceai_N_a_pred_prob = pickle.load(f)
        spliceai_N_a_label_prob = pickle.load(f)

        spliceai_N_d_pred_prob = [x.numpy() for x in spliceai_N_d_pred_prob]
        spliceai_N_a_pred_prob = [x.numpy() for x in spliceai_N_a_pred_prob]
        spliceai_N_d_label_prob = [1]*3000+[0]*(len(spliceai_N_d_pred_prob)-3000)
        spliceai_N_a_label_prob = [1]*3000+[0]*(len(spliceai_N_a_pred_prob)-3000)

        # spliceai_d_label_prob = [float(i) for i in spliceai_d_label_prob]
        # spliceai_a_label_prob = [float(i) for i in spliceai_a_label_prob]

        print("spliceai_N_d_pred_prob : ", spliceai_N_d_pred_prob)
        print("spliceai_N_d_label_prob: ", spliceai_N_d_label_prob)
        print("spliceai_N_a_pred_prob : ", spliceai_N_a_pred_prob)
        print("spliceai_N_a_label_prob: ", spliceai_N_a_label_prob)

        print("spliceai_N_d_pred_prob : ", len(spliceai_N_d_pred_prob))
        print("spliceai_N_d_label_prob: ", len(spliceai_N_d_label_prob))
        print("spliceai_N_a_pred_prob : ", len(spliceai_N_a_pred_prob))
        print("spliceai_N_a_label_prob: ", len(spliceai_N_a_label_prob))

    with open("./INPUT/spliceai_12000.pkl",'rb') as f:
        spliceai_d_pred_prob = pickle.load(f)
        spliceai_d_label_prob = pickle.load(f)
        spliceai_a_pred_prob = pickle.load(f)
        spliceai_a_label_prob = pickle.load(f)

        spliceai_d_pred_prob = [x.numpy() for x in spliceai_d_pred_prob]
        spliceai_a_pred_prob = [x.numpy() for x in spliceai_a_pred_prob]

        # spliceai_d_label_prob = [float(i) for i in spliceai_d_label_prob]
        # spliceai_a_label_prob = [float(i) for i in spliceai_a_label_prob]

        print("spliceai_d_pred_prob : ", spliceai_d_pred_prob)
        print("spliceai_d_label_prob: ", spliceai_d_label_prob)
        print("spliceai_a_pred_prob : ", spliceai_a_pred_prob)
        print("spliceai_a_label_prob: ", spliceai_a_label_prob)

        print("spliceai_d_pred_prob : ", len(spliceai_d_pred_prob))
        print("spliceai_d_label_prob: ", len(spliceai_d_label_prob))
        print("spliceai_a_pred_prob : ", len(spliceai_a_pred_prob))
        print("spliceai_a_label_prob: ", len(spliceai_a_label_prob))


    with open("./INPUT/splam.pkl",'rb') as f:
        splam_d_pred_prob = pickle.load(f)
        splam_d_label_prob = pickle.load(f)
        splam_a_pred_prob = pickle.load(f)
        splam_a_label_prob = pickle.load(f)


        print("splam_d_pred_prob : ", splam_d_pred_prob)
        print("splam_d_label_prob: ", splam_d_label_prob)
        print("splam_a_pred_prob : ", splam_a_pred_prob)
        print("splam_a_label_prob: ", splam_a_label_prob)

        print("splam_d_pred_prob : ", len(splam_d_pred_prob))
        print("splam_d_label_prob: ", len(splam_d_label_prob))
        print("splam_a_pred_prob : ", len(splam_a_pred_prob))
        print("splam_a_label_prob: ", len(splam_a_label_prob))



    spliceai_d_label_prob = np.array(spliceai_d_label_prob)
    spliceai_d_pred_prob = np.array(spliceai_d_pred_prob)
    spliceai_a_pred_prob = np.array(spliceai_a_pred_prob)

    splam_d_label_prob = np.array(splam_d_label_prob)
    splam_d_pred_prob = np.array(splam_d_pred_prob)
    splam_a_pred_prob = np.array(splam_a_pred_prob)
    ###################################
    # Threshold plots
    ###################################
    plot_thresholds(spliceai_d_label_prob, spliceai_d_pred_prob, spliceai_a_pred_prob, "splcieai")

    plot_thresholds(splam_d_label_prob, splam_d_pred_prob, splam_a_pred_prob, "splam")




    # ###################################
    # # PR / ROC of junction
    # ###################################
    # plot_pr_curve_J(spliceai_d_label_prob, spliceai_d_pred_prob, spliceai_a_pred_prob, "spliceai_junc")
    # plot_pr_curve_J(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_pred_prob, "spliceai_N_junc")
    # plot_pr_curve_J(splam_d_label_prob, splam_d_pred_prob, splam_a_pred_prob, "splam_donor")
    # plt.savefig("./IMG/spliceai/junc_pr.png")
    # plt.close()

    # plot_roc_curve_J(spliceai_d_label_prob, spliceai_d_pred_prob, spliceai_a_pred_prob, "spliceai_junc")
    # plot_roc_curve_J(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_pred_prob, "spliceai_N_junc")
    # plot_roc_curve_J(splam_d_label_prob, splam_d_pred_prob, splam_a_pred_prob, "splam_donor")
    # plt.savefig("./IMG/spliceai/junc_roc.png")
    # plt.close()
    # # print("spliceai_d_label_prob: ", spliceai_d_label_prob)
    # # print("spliceai_d_pred_prob: ", spliceai_d_pred_prob)


    # ###################################
    # # PR / ROC of donor
    # ###################################
    # plot_pr_curve(spliceai_d_label_prob, spliceai_d_pred_prob, "spliceai_donor")
    # plot_pr_curve(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, "spliceai_N_donor")
    # plot_pr_curve(splam_d_label_prob, splam_d_pred_prob, "splam_donor")
    # # plt.show()
    # plt.savefig("./IMG/spliceai/donor_pr.png")
    # plt.close()

    # plot_roc_curve(spliceai_d_label_prob, spliceai_d_pred_prob, "spliceai_donor")
    # plot_roc_curve(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, "spliceai_N_donor")
    # plot_roc_curve(splam_d_label_prob, splam_d_pred_prob, "splam_donor")
    # # plt.show()
    # plt.savefig("./IMG/spliceai/donor_roc.png")
    # plt.close()


    # ###################################
    # # PR / ROC of acceptor
    # ###################################
    # plot_pr_curve(spliceai_a_label_prob, spliceai_a_pred_prob, "spliceai_acceptor")
    # plot_pr_curve(spliceai_N_a_label_prob, spliceai_N_a_pred_prob, "spliceai_N_acceptor")
    # plot_pr_curve(splam_a_label_prob, splam_a_pred_prob, "splam_acceptor")
    # # plt.show()
    # plt.savefig("./IMG/spliceai/acceptor_pr.png")
    # plt.close()

    # plot_roc_curve(spliceai_a_label_prob, spliceai_a_pred_prob, "spliceai_acceptor")
    # plot_roc_curve(spliceai_N_a_label_prob, spliceai_N_a_pred_prob, "spliceai_N_acceptor")
    # plot_roc_curve(splam_a_label_prob, splam_a_pred_prob, "splam_acceptor")
    # # plt.show()
    # plt.savefig("./IMG/spliceai/acceptor_roc.png")
    # plt.close()


if __name__ == "__main__":
    main()