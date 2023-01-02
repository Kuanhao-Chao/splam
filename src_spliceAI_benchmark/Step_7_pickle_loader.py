import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve

def plot_pr_curve(true_y, y_prob, label):
    """
    plots the roc curve based of the probabilities
    """
    precision, recall, thresholds = precision_recall_curve_self(true_y, y_prob, label)
    # for threshold in thresholds:
    #     print("threshold: ", threshold)
    plt.plot(precision, recall, label=label)
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
    plt.plot(precision, recall, label=label)
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

# 0.1 - 0.9 (9)
D_TP = [0]*9
D_TN = [0]*9
D_FP = [0]*9
D_FN = [0]*9

A_TP = [0]*9
A_TN = [0]*9
A_FP = [0]*9
A_FN = [0]*9

J_TP = [0]*9
J_TN = [0]*9
J_FP = [0]*9
J_FN = [0]*9


def threshold_assessment(donor_p, acceptor_p, donor_l, acceptor_l):
    for idx in range(0,9,1):
        threshold = (idx+1)*0.1
        print("threshold:  ", threshold)
        # print("Donor   : ", np.where(Y_pred[:, :, 2] > threshold))
        # print("Acceptor: ", np.where(Y_pred[:, :, 1] > threshold))

        if donor_l[idx] == 1 and acceptor_l[idx] == 1:
            # print("INNNN!!POS_NUM")
            if donor_p > threshold:
                D_TP[idx]+=1
            else:
                D_FN[idx]+=1

            if acceptor_p > threshold:
                A_TP[idx]+=1
            else:
                A_FN[idx]+=1

            if donor_p > threshold and acceptor_p > threshold:
                J_TP[idx]+=1
            else:
                J_FN[idx]+=1
        else:
            # print("OUTOUT!!POS_NUM")
            if donor_p > threshold:
                D_FP[idx]+=1
            else:
                D_TN[idx]+=1

            if acceptor_p > threshold:
                A_FP[idx]+=1
            else:
                A_TN[idx]+=1

            if donor_p > threshold and acceptor_p > threshold:
                J_FP[idx]+=1
            else:
                J_TN[idx]+=1


def main():

    ratio = 12000
    os.makedirs("./IMG/spliceai/", exist_ok=True)
    os.makedirs("./IMG/splam/", exist_ok=True)
    
    d_pred_prob = []
    d_label_prob = []
    a_pred_prob = []
    a_label_prob = []

    splam_d_pred_prob = []
    splam_d_label_prob = []
    splam_a_pred_prob = []
    splam_a_label_prob = []
    with open("./INPUT/spliceai_"+str(ratio)+".pkl",'rb') as f:
        d_pred_prob = pickle.load(f)
        d_label_prob = pickle.load(f)[:ratio//2]
        a_pred_prob = pickle.load(f)
        a_label_prob = pickle.load(f)[:ratio//2]

        d_pred_prob = [x.numpy() for x in d_pred_prob]
        a_pred_prob = [x.numpy() for x in a_pred_prob]

        # d_label_prob = [float(i) for i in d_label_prob]
        # a_label_prob = [float(i) for i in a_label_prob]

        print("d_pred_prob : ", d_pred_prob)
        print("d_label_prob: ", d_label_prob)
        print("a_pred_prob : ", a_pred_prob)
        print("a_label_prob: ", a_label_prob)

        print("d_pred_prob : ", len(d_pred_prob))
        print("d_label_prob: ", len(d_label_prob))
        print("a_pred_prob : ", len(a_pred_prob))
        print("a_label_prob: ", len(a_label_prob))
        

    with open("./INPUT/splam.pkl",'rb') as f:
        splam_d_pred_prob = pickle.load(f)
        splam_d_label_prob = pickle.load(f)
        splam_a_pred_prob = pickle.load(f)
        splam_a_label_prob = pickle.load(f)

        # splam_d_pred_prob = [x.numpy() for x in splam_d_pred_prob]
        # splam_a_pred_prob = [x.numpy() for x in splam_a_pred_prob]

        # d_label_prob = [float(i) for i in d_label_prob]
        # a_label_prob = [float(i) for i in a_label_prob]

        # print("d_pred_prob : ", splam_d_pred_prob)
        # print("d_label_prob: ", splam_d_label_prob)
        # print("a_pred_prob : ", splam_a_pred_prob)
        # print("a_label_prob: ", splam_a_label_prob)

        # print("d_pred_prob : ", len(splam_d_pred_prob))
        # print("d_label_prob: ", len(splam_d_label_prob))
        # print("a_pred_prob : ", len(splam_a_pred_prob))
        # print("a_label_prob: ", len(splam_a_label_prob))
        

    plot_pr_curve_J(d_label_prob, d_pred_prob, a_pred_prob, "spliceai_junc")
    plot_pr_curve_J(splam_d_label_prob, splam_d_pred_prob, splam_a_pred_prob, "splam_donor")
    plt.savefig("./IMG/spliceai/junc_pr.png")
    plt.close()

    plot_roc_curve_J(d_label_prob, d_pred_prob, a_pred_prob, "spliceai_junc")
    plot_roc_curve_J(splam_d_label_prob, splam_d_pred_prob, splam_a_pred_prob, "splam_donor")
    plt.savefig("./IMG/spliceai/junc_roc.png")
    plt.close()

    plot_pr_curve(d_label_prob, d_pred_prob, "spliceai_donor")
    plot_pr_curve(splam_d_label_prob, splam_d_pred_prob, "splam_donor")
    # plt.show()
    plt.savefig("./IMG/spliceai/donor_pr.png")
    plt.close()

    plot_roc_curve(d_label_prob, d_pred_prob, "spliceai_donor")
    plot_roc_curve(splam_d_label_prob, splam_d_pred_prob, "splam_donor")
    # plt.show()
    plt.savefig("./IMG/spliceai/donor_roc.png")
    plt.close()

    plot_pr_curve(a_label_prob, a_pred_prob, "spliceai_acceptor")
    plot_pr_curve(splam_a_label_prob, splam_a_pred_prob, "splam_acceptor")
    # plt.show()
    plt.savefig("./IMG/spliceai/acceptor_pr.png")
    plt.close()

    plot_roc_curve(a_label_prob, a_pred_prob, "spliceai_acceptor")
    plot_roc_curve(splam_a_label_prob, splam_a_pred_prob, "splam_acceptor")
    # plt.show()
    plt.savefig("./IMG/spliceai/acceptor_roc.png")
    plt.close()


    # for i in range(1,len(Y_ALL)):
    #     print("i: ", i)
    #     print(np.array(Y_ALL[i]).shape)
    #     # print(np.array(Y_ALL[i]))

    # threshold_assessment(d_pred_prob, a_pred_prob, d_label_prob, a_label_prob)

    # print("D_TP: ", D_TP)
    # print("D_TN: ", D_TN)
    # print("D_FP: ", D_FP)
    # print("D_FN: ", D_FN)

    # print("A_TP: ", A_TP)
    # print("A_TN: ", A_TN)
    # print("A_FP: ", A_FP)
    # print("A_FN: ", A_FN)

    # print("J_TP: ", J_TP)
    # print("J_TN: ", J_TN)
    # print("J_FP: ", J_FP)
    # print("J_FN: ", J_FN)



if __name__ == "__main__":
    main()