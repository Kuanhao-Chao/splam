import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
import pandas as pd
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

def plot_pr_curve_J(true_y, y_prob_d, y_prob_a, label, option, choice, color):
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
    label_res = label + "  ${\\bf  (" + str('%.4f' % auc_precision_recall) + ")}$"
    plt_res = plt.plot(recall, precision, label=label_res, color=color)
    # plt.plot(recall, precision, label=f"{label:<50}{'%10.4f' % auc_precision_recall:>1}")

    # , marker='o')
    # plt.legend(fontsize=8)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    return plt_res, label_res


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


def plot_roc_curve_J(true_y, y_prob_d, y_prob_a, label, option, choice, color):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = roc_curve_J_level_self(true_y, y_prob_d, y_prob_a, label, option, choice)
    auc_roc = auc(fpr, tpr)
    print("auc_roc: ", auc_roc)
    label_res = r""+label + "  ${\\bf  (" + str('%.4f' % auc_roc) + ")}$"
    plt_res = plt.plot(fpr, tpr, label=label_res, color=color)
            #  "$\textbf{(}$" + str('%.4f' % auc_roc) + ")}")
    # , marker='o')
    # plt.legend(fontsize=8)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    return plt_res, label_res


def main(db):

    #####################################
    # Declaring parameters for probability & prediction array
    #####################################

    ### noN ###
    # positive
    noN_pos_df = pd.read_csv(f'../../src_pos_test/output/aggregate/avg_data.noN.{db}.csv')
    noN_pos_df['label'] = 1
    # negative
    noN_neg_df = pd.read_csv(f'../data/aggregate/avg_data.noN.{db}.csv')
    noN_neg_df['label'] = 0
    
    noN_merge_df = pd.concat([noN_pos_df, noN_neg_df], axis=1)

    ### N ###
    # positive
    N_pos_df = pd.read_csv(f'../../src_pos_test/output/aggregate/avg_data.N.{db}.csv')
    N_pos_df['label'] = 1
    # negative
    N_neg_df = pd.read_csv(f'../data/aggregate/avg_data.N.{db}.csv')
    N_neg_df['label'] = 0
    
    N_merge_df = pd.concat([N_pos_df, N_neg_df], axis=1)
    
    # load labels and predictions from dataframe
    spliceai_noN_d_label = noN_merge_df['label']
    spliceai_noN_a_label = noN_merge_df['label']
    spliceai_noN_d_pred = noN_merge_df['d_score_spliceai']
    spliceai_noN_a_pred = noN_merge_df['a_score_spliceai']

    spliceai_N_d_label = N_merge_df['label']
    spliceai_N_a_label = N_merge_df['label']
    spliceai_N_d_pred = N_merge_df['d_score_spliceai']
    spliceai_N_a_pred = N_merge_df['a_score_spliceai']

    splam_d_label = N_merge_df['label']
    splam_a_label = N_merge_df['label']
    splam_d_pred = N_merge_df['d_score_splam']
    splam_a_pred = N_merge_df['a_score_splam']

    print(f'SpliceAI noN_d labels:\n\t1:{len(spliceai_noN_d_label[spliceai_noN_d_label==1])}\t0:{len(spliceai_noN_d_label[spliceai_noN_d_label==0])}')
    print(f'SpliceAI noN_a labels:\n\t1:{len(spliceai_noN_a_label[spliceai_noN_a_label==1])}\t0:{len(spliceai_noN_a_label[spliceai_noN_a_label==0])}')
    print(f'SpliceAI N_d labels:\n\t1:{len(spliceai_N_d_label[spliceai_N_d_label==1])}\t0:{len(spliceai_N_d_label[spliceai_N_d_label==0])}')
    print(f'SpliceAI N_a labels:\n\t1:{len(spliceai_N_a_label[spliceai_N_a_label==1])}\t0:{len(spliceai_N_a_label[spliceai_N_a_label==0])}')
    print(f'SPLAM d labels:\n\t1:{len(splam_d_label[splam_d_label==1])}\t0:{len(splam_d_label[splam_d_label==0])}')
    print(f'SPLAM a labels:\n\t1:{len(splam_a_label[splam_a_label==1])}\t0:{len(splam_a_label[splam_a_label==0])}')
    
    

    for MANE_OR_ALTS in ["pos_MANE", "pos_ALTS"]:
        for SPLICEAI_VERSION in ["1", "2", "3", "4", "5", "AVERAGE"]:
            with open("../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da.N.merged."+MANE_OR_ALTS+".pkl", "rb") as fr:
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

            with open("../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da.noN.merged."+MANE_OR_ALTS+".pkl", "rb") as fr:
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

                with open("../../src_tools_evaluation/splam_result/"+SPLAM_VERSION+"/splam.da.noshuffle.merged."+MANE_OR_ALTS+".pkl",'rb') as f:
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

                # with open("../../src_tools_evaluation/splam_result/SPLAM_v12/splam.da.noshuffle.merged.pkl",'rb') as f:
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
                spliceainoN_plt, spliceainoN_label = plot_pr_curve_J(spliceai_noN_d_label, spliceai_noN_d_pred,spliceai_noN_a_pred, "SpliceAI", "sklearn", "min", color="#ff7f0e")
                spliceaiN_plt, spliceaiN_label = plot_pr_curve_J(spliceai_N_d_label, spliceai_N_d_pred,spliceai_N_a_pred, "SpliceAI-10k-Ns", "sklearn", "min", color="#1f77b4")
                splam_plt, splam_label = plot_pr_curve_J(splam_d_label, splam_d_pred, splam_a_pred, "SPLAM", "sklearn", "min", color="#2ca02c")

                ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
                ax.plot([1, 0], [0, 1], transform=ax.transAxes, linestyle='dashed', color="#d62728")
                #get handles and labels
                handles, labels = plt.gca().get_legend_handles_labels()

                #specify order of items in legend
                order = [2, 0 ,1]

                #add legend to plot
                plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=8)
                # plt.legend([splam_plt, spliceainoN_plt, spliceaiN_plt], [splam_label, spliceainoN_label, spliceaiN_label], fontsize=8)
                plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/junction/junc_pr_min_ratio_2000-10000_"+MANE_OR_ALTS+".png", bbox_inches='tight', dpi=300)
                plt.close()



                # ###################################
                # # PR of junction (avg)
                # ###################################
                # fig, ax = plt.subplots()
                # plot_pr_curve_J(splam_d_label, splam_d_pred, splam_a_pred, "SPLAM", "sklearn", "avg")
                # plot_pr_curve_J(spliceai_noN_d_label, spliceai_noN_d_pred, spliceai_noN_a_pred, "SpliceAI", "sklearn", "avg")
                # plot_pr_curve_J(spliceai_N_d_label, spliceai_N_d_pred, spliceai_N_a_pred, "SpliceAI-10k-Ns", "sklearn", "avg")

                # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
                # ax.plot([1, 0], [0, 1], transform=ax.transAxes, linestyle='dashed')
                # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/junction/junc_pr_avg.png", bbox_inches='tight', dpi=300)
                # plt.close()
                

                ###################################
                # ROC of junction (min)
                ###################################
                fig, ax = plt.subplots()
                spliceainoN_plt, spliceainoN_label = plot_roc_curve_J(spliceai_noN_d_label, spliceai_noN_d_pred, spliceai_noN_a_pred, "SpliceAI", "sklearn", "min", color="#ff7f0e")
                spliceaiN_plt, spliceaiN_label = plot_roc_curve_J(spliceai_N_d_label, spliceai_N_d_pred, spliceai_N_a_pred, "SpliceAI-10k-Ns", "sklearn", "min", color="#1f77b4")
                splam_plt, splam_label = plot_roc_curve_J(splam_d_label, splam_d_pred, splam_a_pred, "SPLAM", "sklearn", "min", color="#2ca02c")

                ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
                ax.plot([0, 1], [0, 1], transform=ax.transAxes, linestyle='dashed', color="#d62728")

                #get handles and labels
                handles, labels = plt.gca().get_legend_handles_labels()

                #specify order of items in legend
                order = [2, 0 ,1]

                #add legend to plot
                plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=8)
                # plt.legend([splam_plt, spliceainoN_plt, spliceaiN_plt], [splam_label, spliceainoN_label, spliceaiN_label], fontsize=8)
                plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/junction/junc_roc_min_ratio_2000-10000_"+MANE_OR_ALTS+".png", bbox_inches='tight', dpi=300)
                plt.close()

                # ################################### 
                # # ROC of junction (avg)
                # ###################################
                # fig, ax = plt.subplots()
                # plot_roc_curve_J(splam_d_label, splam_d_pred, splam_a_pred, "SPLAM", "sklearn", "avg")
                # plot_roc_curve_J(spliceai_noN_d_label, spliceai_noN_d_pred, spliceai_noN_a_pred, "SpliceAI", "sklearn", "avg")
                # plot_roc_curve_J(spliceai_N_d_label, spliceai_N_d_pred, spliceai_N_a_pred, "SpliceAI-10k-Ns", "sklearn", "avg")

                # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
                # ax.plot([0, 1], [0, 1], transform=ax.transAxes, linestyle='dashed')
                # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/junction/junc_roc_avg.png", bbox_inches='tight', dpi=300)
                # plt.close()


                # ###################################
                # # PR / ROC of donor
                # ###################################
                # plot_pr_curve(spliceai_noN_d_label, spliceai_noN_d_pred, "spliceai_donor", "self")
                # plot_pr_curve(spliceai_N_d_label, spliceai_N_d_pred, "spliceai_N_donor", "self")
                # plot_pr_curve(splam_d_label, splam_d_pred, "splam_donor", "self")
                # # plt.show()
                # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/donor/donor_pr.png")
                # plt.close()

                # plot_roc_curve(spliceai_noN_d_label, spliceai_noN_d_pred, "spliceai_donor", "self")
                # plot_roc_curve(spliceai_N_d_label, spliceai_N_d_pred, "spliceai_N_donor", "self")
                # plot_roc_curve(splam_d_label, splam_d_pred, "splam_donor", "self")
                # # plt.show()
                # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/donor/donor_roc.png")
                # plt.close()

                # ###################################
                # # PR / ROC of acceptor
                # ###################################
                # plot_pr_curve(spliceai_noN_a_label_prob, spliceai_noN_a_pred, "spliceai_acceptor", "self")
                # plot_pr_curve(spliceai_N_a_label_prob, spliceai_N_a_pred, "spliceai_N_acceptor", "self")
                # plot_pr_curve(splam_a_label, splam_a_pred, "splam_acceptor", "self")
                # # plt.show()
                # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/acceptor/acceptor_pr.png")
                # plt.close()

                # plot_roc_curve(spliceai_noN_a_label_prob, spliceai_noN_a_pred, "spliceai_acceptor", "self")
                # plot_roc_curve(spliceai_N_a_label_prob, spliceai_N_a_pred, "spliceai_N_acceptor", "self")
                # plot_roc_curve(splam_a_label, splam_a_pred, "splam_acceptor", "self")
                # # plt.show()
                # plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/acceptor/acceptor_roc.png")
                # plt.close()

if __name__ == "__main__":

    db = 'GRCm39'

    main(db)