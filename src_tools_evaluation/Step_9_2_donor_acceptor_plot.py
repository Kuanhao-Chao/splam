import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
from util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, PrecisionRecallDisplay
import plotly.express as px 
import plotly.graph_objects as go
import pandas as pd

def plot_scatter_plot(label_d, score_d, label_a, score_a, filename):
    junc_prob = label_d.astype(bool)
    non_junc_prob = (1-label_d).astype(bool)
    print("\tspliceai_junc_prob: ", junc_prob)
    print("\tspliceai_non_junc_prob: ", non_junc_prob)

    # fig, ax = plt.subplots()
    # ax.set_title(filename)
    # ax.set_xlabel("Donor score")
    # ax.set_ylabel("Acceptor score")

    df = pd.DataFrame({'donor': score_d, 'acceptor': score_a, 'junciton': junc_prob})
    fig = px.scatter(df, x='donor',y='acceptor', color='junciton')
    fig.show()


# def plot_min_linear_plot(label_d, score_d, label_a, score_a, filename):
#     junction_score = np.minimum(score_d, score_a)

#     junc_prob = label_d.astype(bool)
#     non_junc_prob = (1-label_d).astype(bool)
#     print("\tspliceai_junc_prob: ", junc_prob)
#     print("\tspliceai_non_junc_prob: ", non_junc_prob)

#     fig, ax = plt.subplots()
#     ax.set_title(filename)
#     # ax.set_xlabel("Donor score")
#     # ax.set_ylabel("Acceptor score")

#     ar = np.arange(10) # just as an example array
#     junc_legend = ax.plot(junction_score[junc_prob], np.zeros_like(junction_score[junc_prob]) + 0.)
#                         #   , s = 1)
#     non_junc_legend = ax.plot(junction_score[non_junc_prob], np.zeros_like(junction_score[non_junc_prob]) + 0.)
#                             #   , s = 1)
#     # junc_legend = ax.scatter(score_d[junc_prob], score_a[junc_prob], s = 1)
#     # non_junc_legend = ax.scatter(score_d[non_junc_prob], score_a[non_junc_prob], s = 1)
#     ax.legend([junc_legend, non_junc_legend], ['Junction', 'Non Junction'])
#     fig.savefig(filename)
    

def main():
    #####################################
    # Creating directories for visualization.
    #####################################
    figure_root = "./IMG/d_a/"
    os.makedirs(figure_root, exist_ok=True)

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


    splam_S_d_pred_prob = []
    splam_S_d_label_prob = []
    splam_S_a_pred_prob = []
    splam_S_a_label_prob = []

    splam_noS_d_pred_prob = []
    splam_noS_d_label_prob = []
    splam_noS_a_pred_prob = []
    splam_noS_a_label_prob = []

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



    with open("./splam_result/splam.da.shuffle.pkl",'rb') as f:
        splam_S_d_label_prob = pickle.load(f)
        splam_S_d_pred_prob = pickle.load(f)
        
        splam_S_a_label_prob = pickle.load(f)
        splam_S_a_pred_prob = pickle.load(f)

        print("splam_S_d_label_prob : ", splam_S_d_label_prob)
        print("splam_S_d_pred_prob: ", splam_S_d_pred_prob)

        print("splam_S_a_label_prob : ", splam_S_a_label_prob)
        print("splam_S_a_pred_prob: ", splam_S_a_pred_prob)

    
    with open("./splam_result/splam.da.noshuffle.pkl",'rb') as f:
        splam_noS_d_label_prob = pickle.load(f)
        splam_noS_d_pred_prob = pickle.load(f)
        
        splam_noS_a_label_prob = pickle.load(f)
        splam_noS_a_pred_prob = pickle.load(f)

        print("splam_noS_d_label_prob : ", splam_noS_d_label_prob)
        print("splam_noS_d_pred_prob: ", splam_noS_d_pred_prob)

        print("splam_noS_a_label_prob : ", splam_noS_a_label_prob)
        print("splam_noS_a_pred_prob: ", splam_noS_a_pred_prob)

    # with open("./spliceai_result/splam.nobatch.0.100.pkl",'rb') as f:
    #     splam_j_nobatch_pred_prob = pickle.load(f)
    #     splam_j_nobatch_label_prob = pickle.load(f)
    #     print("splam_j_nobatch_pred_prob : ", splam_j_nobatch_pred_prob)
    #     print("splam_j_nobatch_label_prob: ", splam_j_nobatch_label_prob)
    #     print("splam_j_nobatch_pred_prob : ", len(splam_j_nobatch_pred_prob))
    #     print("splam_j_nobatch_label_prob: ", len(splam_j_nobatch_label_prob))

    splam_S_d_label_prob = np.array(splam_S_d_label_prob)
    splam_S_d_pred_prob = np.array(splam_S_d_pred_prob)
    splam_S_a_label_prob = np.array(splam_S_a_label_prob)
    splam_S_a_pred_prob = np.array(splam_S_a_pred_prob)

    splam_noS_d_label_prob = np.array(splam_noS_d_label_prob)
    splam_noS_d_pred_prob = np.array(splam_noS_d_pred_prob)
    splam_noS_a_label_prob = np.array(splam_noS_a_label_prob)
    splam_noS_a_pred_prob = np.array(splam_noS_a_pred_prob)





    ################################################
    # spliceai scatter plot
    ################################################
    # plot_scatter_plot(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, spliceai_noN_a_label_prob, spliceai_noN_a_pred_prob, figure_root+"spliceai_noN.png")
    # plot_min_linear_plot(spliceai_noN_d_label_prob, spliceai_noN_d_pred_prob, spliceai_noN_a_label_prob, spliceai_noN_a_pred_prob, figure_root+"spliceai_noN_1d.png")
    # plot_scatter_plot(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_label_prob, spliceai_N_a_pred_prob, figure_root+"spliceai_N.png")
    # # plot_min_linear_plot(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_label_prob, spliceai_N_a_pred_prob, figure_root+"spliceai_N_1d.png")

    plot_scatter_plot(splam_S_d_label_prob, splam_S_d_pred_prob, splam_S_a_label_prob, splam_S_a_pred_prob, figure_root+"splam_shuffle.png")
    # # plot_min_linear_plot(splam_S_d_label_prob, splam_S_d_pred_prob, splam_S_a_label_prob, splam_S_a_pred_prob, figure_root+"splam_shuffle_1d.png")

    # plot_scatter_plot(splam_noS_d_label_prob, splam_noS_d_pred_prob, splam_noS_a_label_prob, splam_noS_a_pred_prob, figure_root+"splam_noshuffle.png")
    # # plot_min_linear_plot(splam_noS_d_label_prob, splam_noS_d_pred_prob, splam_noS_a_label_prob, splam_noS_a_pred_prob, figure_root+"splam_noshuffle_1d.png")


if __name__ == "__main__":
    main()


