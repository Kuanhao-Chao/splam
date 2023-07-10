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

    mkdir = lambda p : os.makedirs(os.path.dirname(p), exist_ok=True)

    POS_NUM = 2000
    NEG_NUM = 18000
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
    
    noN_pos_df = noN_pos_df.sample(n=POS_NUM, random_state=1)
    noN_neg_df = noN_neg_df.sample(n=NEG_NUM, random_state=1)

    noN_merge_df = pd.concat([noN_pos_df, noN_neg_df], axis=0)


    print("noN_pos_df: ", len(noN_pos_df))
    print("noN_neg_df: ", len(noN_neg_df))
    print("noN_merge_df: ", len(noN_merge_df))


    
    
    ### N ###
    # positive
    N_pos_df = pd.read_csv(f'../../src_pos_test/output/aggregate/avg_data.N.{db}.csv')
    N_pos_df['label'] = 1
    # negative
    N_neg_df = pd.read_csv(f'../data/aggregate/avg_data.N.{db}.csv')
    N_neg_df['label'] = 0
    
    N_pos_df = N_pos_df.sample(n=POS_NUM, random_state=1)
    N_neg_df = N_neg_df.sample(n=NEG_NUM, random_state=1)

    N_merge_df = pd.concat([N_pos_df, N_neg_df], axis=0)
    
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

    print(f'SpliceAI noN_d labels:\n\t1:{len(spliceai_noN_d_label[spliceai_noN_d_label==1])}\n\t0:{len(spliceai_noN_d_label[spliceai_noN_d_label==0])}')
    print(f'SpliceAI noN_a labels:\n\t1:{len(spliceai_noN_a_label[spliceai_noN_a_label==1])}\n\t0:{len(spliceai_noN_a_label[spliceai_noN_a_label==0])}')
    print(f'SpliceAI N_d labels:\n\t1:{len(spliceai_N_d_label[spliceai_N_d_label==1])}\n\t0:{len(spliceai_N_d_label[spliceai_N_d_label==0])}')
    print(f'SpliceAI N_a labels:\n\t1:{len(spliceai_N_a_label[spliceai_N_a_label==1])}\n\t0:{len(spliceai_N_a_label[spliceai_N_a_label==0])}')
    print(f'SPLAM d labels:\n\t1:{len(splam_d_label[splam_d_label==1])}\n\t0:{len(splam_d_label[splam_d_label==0])}')
    print(f'SPLAM a labels:\n\t1:{len(splam_a_label[splam_a_label==1])}\n\t0:{len(splam_a_label[splam_a_label==0])}')

    print(f'SpliceAI noN_d pred {len(spliceai_noN_d_pred)}')
    print(f'SpliceAI noN_a pred {len(spliceai_noN_a_pred)}')
    print(f'SpliceAI N_d pred {len(spliceai_N_d_pred)}')
    print(f'SpliceAI N_d pred {len(spliceai_N_a_pred)}')
    print(f'SPLAM d pred {len(splam_d_pred)}')
    print(f'SPLAM a pred {len(splam_a_pred)}')  


    ###################################
    # PR of junction (min)
    ###################################
    print("Plotting junction!!")
    fig, ax = plt.subplots()
    spliceainoN_plt, spliceainoN_label = plot_pr_curve_J(spliceai_noN_d_label, spliceai_noN_d_pred, spliceai_noN_a_pred, "SpliceAI", "sklearn", "min", color="#ff7f0e")
    spliceaiN_plt, spliceaiN_label = plot_pr_curve_J(spliceai_N_d_label, spliceai_N_d_pred, spliceai_N_a_pred, "SpliceAI-10k-Ns", "sklearn", "min", color="#1f77b4")
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
    path = f'../figures/pr/{db}_junc_pr_min_ratio_2000-10000.png'
    mkdir(path)
    # plt.savefig(path, bbox_inches='tight', dpi=300)
    # plt.close()
    plt.show()



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
    path = f'../figures/roc/{db}_junc_roc_min_ratio_2000-10000.png'
    mkdir(path)
    
    plt.show()


    
    # plt.savefig(path, bbox_inches='tight', dpi=300)
    # plt.close()

if __name__ == "__main__":

    databases = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    idxs = [1,2,3]
    
    for idx in idxs:
        main(databases[idx])