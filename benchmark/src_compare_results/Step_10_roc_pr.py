import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
from util import *
from sklearn.metrics import auc, accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, PrecisionRecallDisplay
from sklearn import svm

def plot_pr_curve(true_y, y_prob, label, option):
    """
    plots the roc curve based of the probabilities
    """
    # precision, recall, thresholds = precision_recall_curve_self(true_y, y_prob, label, option)
    precision, recall, thresholds = precision_recall_curve(true_y, y_prob)
    plt.plot(recall, precision, label=label)
    # , label=label, marker='o')
    plt.legend()
    plt.xlabel('Recall')
    plt.ylabel('Precision')

def plot_pr_curve_J(true_y, y_prob_d, y_prob_a, label, option, choice):
    """
    plots the roc curve based of the probabilities
    """
    # precision, recall, thresholds = precision_recall_curve_J_level_self(true_y, y_prob_d, y_prob_a, label, option)

    print("true_y: ", true_y)
    print("y_prob_d: ", y_prob_d)
    print("y_prob_a: ", y_prob_a)
    precision, recall, thresholds = precision_recall_curve_J_level_self(true_y, y_prob_d, y_prob_a, label, option, choice)
    auc_precision_recall = auc(recall, precision)
    print("auc_precision_recall: ", auc_precision_recall)
    plt.plot(recall, precision, label=label + "  ${\\bf  (" + str('%.4f' % auc_precision_recall) + ")}$")
    # plt.plot(recall, precision, label=f"{label:<50}{'%10.4f' % auc_precision_recall:>1}")

    # , marker='o')
    plt.legend(fontsize=8)
    plt.xlabel('Recall')
    plt.ylabel('Precision')

def plot_roc_curve(true_y, y_prob, label, option):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = roc_curve(true_y, y_prob)
    # fpr, tpr, thresholds = ROC_curve_self(true_y, y_prob, label, option)
    auc_roc = auc(fpr, tpr)
    plt.plot(fpr, tpr, label=label + ": "+str(auc_roc))
    # , marker='o' )
    plt.legend()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

def plot_roc_curve_J(true_y, y_prob_d, y_prob_a, label, option, choice):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = roc_curve_J_level_self(true_y, y_prob_d, y_prob_a, label, option, choice)
    auc_roc = auc(fpr, tpr)
    print("auc_roc: ", auc_roc)
    plt.plot(fpr, tpr, label=r""+label + "  ${\\bf  (" + str('%.4f' % auc_roc) + ")}$")
            #  "$\textbf{(}$" + str('%.4f' % auc_roc) + ")}")
    # , marker='o')
    plt.legend(fontsize=8)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')


def plot_thresholds(true_y, y_prob_d, y_prob_a, type, SPLAM_VERSION):
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
    plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/T_histogram/"+type+"_donor_thresholds.png")
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
    plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/T_histogram/"+type+"_acceptor_thresholds.png")
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
    plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/T_histogram/"+type+"_junction_thresholds.png")
    plt.close()
        

def plot_thresholds_J(true_y, y_prob_j, type, SPLAM_VERSION):
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
    plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/T_histogram/"+type+"_junction_thresholds.png")
    plt.close()


def main():
    # for SPLAM_VERSION in ["SPLAM_v11", "SPLAM_v12"]:

    #####################################
    # Declaring parameters for probability & prediction array
    #####################################

    for SPLICEAI_VERSION in ["1", "2", "3", "4", "5"]:
        for MANE_OR_ALTS in ["pos_MANE", "pos_ALTS"]:
            with open("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da.N.merged."+MANE_OR_ALTS+".pkl", "rb") as fr:
                spliceai_N_d_label = pickle.load(fr)
                spliceai_N_d_pred = pickle.load(fr)
                spliceai_N_a_label = pickle.load(fr)
                spliceai_N_a_pred = pickle.load(fr)
                # j_pred_prob = [x.numpy() for x in j_pred_prob]
                # j_pred_prob = [x.numpy() for x in j_pred_prob]
                print("\tspliceai_N_d_label : ", len(spliceai_N_d_label))
                print("\tspliceai_N_d_pred: ", len(spliceai_N_d_pred))
                print("\tspliceai_N_d_label : ", len(spliceai_N_d_label[spliceai_N_d_label == 1]))
                print("\tspliceai_N_d_pred: ", spliceai_N_d_pred)
                print("")
                print("\tspliceai_N_a_label : ", len(spliceai_N_a_label))
                print("\tspliceai_N_a_pred: ", len(spliceai_N_a_pred))
                print("\tspliceai_N_a_pred: ", len(spliceai_N_a_label[spliceai_N_a_label == 1]))
                print("\tspliceai_N_a_pred: ", spliceai_N_a_pred)
                print("")

            with open("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da.noN.merged."+MANE_OR_ALTS+".pkl", "rb") as fr:
                spliceai_noN_d_label = pickle.load(fr)
                spliceai_noN_d_pred = pickle.load(fr)
                spliceai_noN_a_label = pickle.load(fr)
                spliceai_noN_a_pred = pickle.load(fr)

                print("\tspliceai_noN_d_label : ", len(spliceai_noN_d_label))
                print("\tspliceai_noN_d_pred: ", len(spliceai_noN_d_pred))
                print("\tspliceai_noN_d_pred: ", spliceai_noN_d_pred)
                print("")
                print("\tspliceai_noN_a_label : ", len(spliceai_noN_a_label))
                print("\tspliceai_noN_a_pred: ", len(spliceai_noN_a_pred))
                print("\tspliceai_noN_a_pred: ", spliceai_noN_a_pred)
                print("")


            for SPLAM_VERSION in ["SPLAM_v11"]:#, "SPLAM_v12"]:
                #####################################
                # Creating directories for visualization.
                #####################################
                os.makedirs("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/junction/", exist_ok=True)
                os.makedirs("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/T_histogram/", exist_ok=True)

                with open("./splam_result/"+SPLAM_VERSION+"/splam.da.noshuffle.merged."+MANE_OR_ALTS+".pkl",'rb') as f:
                    splam_d_label = pickle.load(f)
                    splam_d_pred = pickle.load(f)
                    splam_a_label = pickle.load(f)
                    splam_a_pred = pickle.load(f)
                    print("\tsplam_d_label : ", len(splam_d_label))
                    print("\tsplam_d_pred: ", len(splam_d_pred))
                    print("")
                    print("\tsplam_a_label : ", len(splam_a_label))
                    print("\tsplam_a_pred: ", len(splam_a_pred))
                    print("")


                # with open("./splam_result/SPLAM_v12/splam.da.noshuffle.merged.pkl",'rb') as f:
                #     splam_v12_d_label_prob = pickle.load(f)
                #     splam_v12_d_pred_prob = pickle.load(f)
                #     splam_v12_a_label_prob = pickle.load(f)
                #     splam_v12_a_pred_prob = pickle.load(f)
                #     print("\tsplam_v12_d_label_prob : ", len(splam_v12_d_label_prob))
                #     print("\tsplam_v12_d_pred_prob: ", len(splam_v12_d_pred_prob))
                #     print("")
                #     print("\tsplam_v12_a_label_prob : ", len(splam_v12_a_label_prob))
                #     print("\tsplam_v12_a_pred_prob: ", len(splam_v12_a_pred_prob))
                #     print("")
                
                # splam_d_label = np.array(splam_d_label)
                # splam_d_pred = np.array(splam_d_pred)
                # splam_a_label = np.array(splam_a_label)
                # splam_a_pred = np.array(splam_a_pred)
                # splam_noS_junc_name = np.array(splam_noS_junc_name)

                # print("splam_d_label : ", (splam_d_label))
                # print("splam_d_pred: ", (splam_d_pred))
                # print("splam_a_label : ", (splam_a_label))
                # print("splam_a_pred: ", (splam_a_pred))
                # # print("splam_noS_junc_name  : ", (splam_noS_junc_name))

                # print("splam_d_label : ", len(splam_d_label))
                # print("splam_d_pred: ", len(splam_d_pred))
                # print("splam_a_label : ", len(splam_a_label))
                # print("splam_a_pred: ", len(splam_a_pred))
                # print("splam_noS_junc_name  : ", len(splam_noS_junc_name))



                ###################################
                # Threshold plots
                ###################################
                # plot_thresholds(spliceai_noN_d_label, spliceai_noN_d_pred, spliceai_noN_a_pred, "splcieai", SPLAM_VERSION)
                # plot_thresholds(spliceai_N_d_label, spliceai_N_d_pred, spliceai_N_a_pred, "splcieai_N", SPLAM_VERSION)

                # plot_thresholds(splam_d_label, splam_d_pred, splam_a_pred, "splam_noshuffle", SPLAM_VERSION)
                # plot_thresholds(splam_S_a_label_prob, splam_S_a_pred_prob, splam_S_d_pred_prob, "splam_shuffle", SPLAM_VERSION)

                # plot_thresholds_J(splam_j_shuffle_label_prob_min, splam_j_shuffle_pred_prob_min, "splam_shuffle", SPLAM_VERSION)
                # plot_thresholds_J(splam_j_noshuffle_label_prob_min, splam_j_noshuffle_pred_prob_min, "splam_noshuffle", SPLAM_VERSION)

                # plot_thresholds_J(splam_j_shuffle_label_prob_avg, splam_j_shuffle_pred_prob_avg, "splam_shuffle_avg", SPLAM_VERSION)
                # plot_thresholds_J(splam_j_noshuffle_label_prob_avg, splam_j_noshuffle_pred_prob_avg, "splam_noshuffle_avg", SPLAM_VERSION)
                # plot_thresholds_J(splam_j_nobatch_label_prob, splam_j_nobatch_pred_prob, "splam_nobatch", SPLAM_VERSION)


                ###################################
                # PR of junction (min)
                ###################################
                print("Plotting junction!!")
                fig, ax = plt.subplots()
                plot_pr_curve_J(splam_d_label, splam_d_pred, splam_a_pred, "SPLAM", "sklearn", "min")
                
                plot_pr_curve_J(spliceai_noN_d_label, spliceai_noN_d_pred,spliceai_noN_a_pred, "SpliceAI with 10k flanking sequence", "sklearn", "min")
                plot_pr_curve_J(spliceai_N_d_label, spliceai_N_d_pred,spliceai_N_a_pred, "SpliceAI with 10k N", "sklearn", "min")

                ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
                ax.plot([1, 0], [0, 1], transform=ax.transAxes, linestyle='dashed')
                plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/junction/junc_pr_min_ratio_2500-10000_"+MANE_OR_ALTS+".png", bbox_inches='tight', dpi=300)
                plt.close()






                # # ###################################
                # # # PR of junction (avg)
                # # ###################################
                # fig, ax = plt.subplots()
                # plot_pr_curve_J(splam_d_label, splam_d_pred, splam_a_pred, "SPLAM", "sklearn", "avg")
                # plot_pr_curve_J(spliceai_noN_d_label, spliceai_noN_d_pred, spliceai_noN_a_pred, "SpliceAI with 10k flanking sequence", "sklearn", "avg")
                # plot_pr_curve_J(spliceai_N_d_label, spliceai_N_d_pred, spliceai_N_a_pred, "SpliceAI with 10k N", "sklearn", "avg")

                # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
                # ax.plot([1, 0], [0, 1], transform=ax.transAxes, linestyle='dashed')
                # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/junction/junc_pr_avg.png", bbox_inches='tight', dpi=300)
                # plt.close()
                

                ###################################
                # ROC of junction (min)
                ###################################
                fig, ax = plt.subplots()
                plot_roc_curve_J(splam_d_label, splam_d_pred, splam_a_pred, "SPLAM", "sklearn", "min")
                plot_roc_curve_J(spliceai_noN_d_label, spliceai_noN_d_pred, spliceai_noN_a_pred, "SpliceAI with 10k flanking sequence", "sklearn", "min")
                plot_roc_curve_J(spliceai_N_d_label, spliceai_N_d_pred, spliceai_N_a_pred, "SpliceAI with 10k N", "sklearn", "min")

                ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
                ax.plot([0, 1], [0, 1], transform=ax.transAxes, linestyle='dashed')
                plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/junction/junc_roc_min_ratio_2500-10000_"+MANE_OR_ALTS+".png", bbox_inches='tight', dpi=300)
                plt.close()

                # ################################### 
                # # ROC of junction (avg)
                # ###################################
                # fig, ax = plt.subplots()
                # plot_roc_curve_J(splam_d_label, splam_d_pred, splam_a_pred, "SPLAM", "sklearn", "avg")
                # plot_roc_curve_J(spliceai_noN_d_label, spliceai_noN_d_pred, spliceai_noN_a_pred, "SpliceAI with 10k flanking sequence", "sklearn", "avg")
                # plot_roc_curve_J(spliceai_N_d_label, spliceai_N_d_pred, spliceai_N_a_pred, "SpliceAI with 10k N", "sklearn", "avg")

                # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
                # ax.plot([0, 1], [0, 1], transform=ax.transAxes, linestyle='dashed')
                # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/junction/junc_roc_avg.png", bbox_inches='tight', dpi=300)
                # plt.close()


                # # ###################################
                # # # PR / ROC of donor
                # # ###################################
                # # plot_pr_curve(spliceai_noN_d_label, spliceai_noN_d_pred, "spliceai_donor", "self")
                # # plot_pr_curve(spliceai_N_d_label, spliceai_N_d_pred, "spliceai_N_donor", "self")
                # # plot_pr_curve(splam_d_label, splam_d_pred, "splam_donor", "self")
                # # # plt.show()
                # # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/donor/donor_pr.png")
                # # plt.close()

                # # plot_roc_curve(spliceai_noN_d_label, spliceai_noN_d_pred, "spliceai_donor", "self")
                # # plot_roc_curve(spliceai_N_d_label, spliceai_N_d_pred, "spliceai_N_donor", "self")
                # # plot_roc_curve(splam_d_label, splam_d_pred, "splam_donor", "self")
                # # # plt.show()
                # # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/donor/donor_roc.png")
                # # plt.close()

                # # ###################################
                # # # PR / ROC of acceptor
                # # ###################################
                # # plot_pr_curve(spliceai_noN_a_label_prob, spliceai_noN_a_pred, "spliceai_acceptor", "self")
                # # plot_pr_curve(spliceai_N_a_label_prob, spliceai_N_a_pred, "spliceai_N_acceptor", "self")
                # # plot_pr_curve(splam_a_label, splam_a_pred, "splam_acceptor", "self")
                # # # plt.show()
                # # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/acceptor/acceptor_pr.png")
                # # plt.close()

                # # plot_roc_curve(spliceai_noN_a_label_prob, spliceai_noN_a_pred, "spliceai_acceptor", "self")
                # # plot_roc_curve(spliceai_N_a_label_prob, spliceai_N_a_pred, "spliceai_N_acceptor", "self")
                # # plot_roc_curve(splam_a_label, splam_a_pred, "splam_acceptor", "self")
                # # # plt.show()
                # # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/acceptor/acceptor_roc.png")
                # # plt.close()

if __name__ == "__main__":
    main()