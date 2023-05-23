import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
from util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, PrecisionRecallDisplay
from scipy.stats import norm
from scipy.stats import gaussian_kde


def plot_scatter_plot(label_d, score_d, label_a, score_a, filename):
    junc_prob = label_d.astype(bool)
    non_junc_prob = (1-label_d).astype(bool)
    # print("\tspliceai_junc_prob: ", junc_prob)
    # print("\tspliceai_junc_prob: ", non_junc_prob)

    fig, ax = plt.subplots()
    ax.set_title(filename)
    ax.set_xlabel("Donor site score")
    ax.set_ylabel("Acceptor site score")


    junc_legend = ax.scatter(score_d[junc_prob], score_a[junc_prob], s = 0.3)
    non_junc_legend = ax.scatter(score_d[non_junc_prob], score_a[non_junc_prob], s = 0.3)    
    ax.legend([junc_legend, non_junc_legend], ['Junction', 'Non Junction'])

    fig.savefig(filename)
    # fig.close()


def plot_min_linear_plot(label_d, score_d, label_a, score_a, filename):
    junction_score = np.minimum(score_d, score_a)

    junc_prob = label_d.astype(bool)
    non_junc_prob = (1-label_d).astype(bool)
    print("\tspliceai_junc_prob: ", junc_prob)
    print("\tspliceai_junc_prob: ", non_junc_prob)

    fig, ax = plt.subplots()
    ax.set_title(filename)
    # ax.set_xlabel("Donor score")
    # ax.set_ylabel("Acceptor score")

    ar = np.arange(10) # just as an example array
    junc_legend = ax.plot(junction_score[junc_prob], np.zeros_like(junction_score[junc_prob]) + 0.)
                        #   , s = 0.3)
    non_junc_legend = ax.plot(junction_score[non_junc_prob], np.zeros_like(junction_score[non_junc_prob]) + 0.)
                            #   , s = 0.3)
    # junc_legend = ax.scatter(score_d[junc_prob], score_a[junc_prob], s = 0.3)
    # non_junc_legend = ax.scatter(score_d[non_junc_prob], score_a[non_junc_prob], s = 0.3)
    ax.legend([junc_legend, non_junc_legend], ['Junction', 'Non Junction'])
    fig.savefig(filename)
    

THRESHOLDS = [0.1, 0.01]



def main():

    MANE_OR_ALTS = ""
    for threshold in THRESHOLDS:
        #####################################
        # Creating directories for visualization.
        #####################################
        # figure_root = "./IMG/d_a/"

        #####################################
        # Declaring parameters for probability & prediction array
        #####################################
        spliceai_N_d_pred_prob = []
        spliceai_N_d_label_prob = []
        spliceai_N_a_pred_prob = []
        spliceai_N_a_label_prob = []

        spliceai_d_pred = []
        spliceai_d_label = []
        spliceai_a_pred = []
        spliceai_a_label = []


        splam_S_d_pred_prob = []
        splam_S_d_label_prob = []
        splam_S_a_pred_prob = []
        splam_S_a_label_prob = []

        splam_d_pred = []
        splam_d_pred = []
        splam_a_pred = []
        splam_noS_a_label_prob = []


        for SPLAM_VERSION in ["SPLAM_v11"]:#, "SPLAM_v12"]:
            for TYPE in ["pos_MANE", "pos_ALTS"]:
                #####################################
                # Creating directories for visualization.
                #####################################


                for TARGET in ["noN", "N"]:
                    figure_root = "./IMG/"+SPLAM_VERSION+"/distance/"
                    target_figure_root = figure_root+TARGET+"/"
                    os.makedirs(target_figure_root, exist_ok=True)
                    

                    with open("../../src_tools_evaluation/spliceai_result_AVERAGE/spliceai.da."+TARGET+".merged."+TYPE+".distance.pkl", "rb") as fr:
                        spliceai_d_label = pickle.load(fr)
                        spliceai_d_pred = pickle.load(fr)
                        spliceai_a_label = pickle.load(fr)
                        spliceai_a_pred = pickle.load(fr)

                        print("\tspliceai_d_label : ", len(spliceai_d_label))
                        print("\tspliceai_d_pred: ", len(spliceai_d_pred))
                        print("\tspliceai_d_pred: ", spliceai_d_pred)
                        print("")
                        print("\tspliceai_a_label : ", len(spliceai_a_label))
                        print("\tspliceai_a_pred: ", len(spliceai_a_pred))
                        print("\tspliceai_a_pred: ", spliceai_a_pred)
                        print("")



                    plt.rcParams['font.size'] = 8
                    # plt.subplots_adjust(wspace=0.4, hspace=-0.5)



                    selected_indices_combined = np.vstack([spliceai_d_pred, spliceai_a_pred])

                    density_combined = gaussian_kde(selected_indices_combined)

                    density_values = density_combined(selected_indices_combined)

                    vmin = -2000 #min(density_1_values.min(), density_2_values.min())-100
                    vmax = density_values.max()
                    print("vmax: ", vmax)
                    norm = plt.Normalize(vmin=vmin, vmax=vmax)
                    cmap = plt.get_cmap('Reds')  # You can choose a different colormap if desired

                    plt.scatter(spliceai_d_pred, spliceai_a_pred, s=0.6, c=density_values, cmap=cmap, alpha=1.0, norm=norm)

                    # Add a colorbar to show the density scale
                    cbar = plt.colorbar()
                    cbar.set_label('Density')

                    plt.xlabel("Donor difference score")
                    plt.ylabel("Acceptor difference score")

                    # Save the plot
                    plt.savefig(target_figure_root + TYPE + ".png", dpi=600, bbox_inches='tight')

                    plt.close()



                        # # fig, axes = plt.subplot_mosaic(
                        # #     [
                        # #         ["main", "TPFN", "FPFN"],
                        # #         ["main", "FPTN", "BLANK"],
                        # #     ],
                        # #     empty_sentinel="BLANK",
                        # #     width_ratios=[3, 1.3, 1.3],
                        # # )

                        # # fig, axes = plt.subplot_mosaic(
                        # #     [
                        # #         ["main", "TPFN", "FPTN", "FPFN"],
                        # #     ],
                        # #     empty_sentinel="BLANK",
                        # #     width_ratios=[3, 1.9, 1.9, 1.9],
                        # # )
                        # # # identify_axes(axes, fontsize=36)
                        
                        
                        # # junc_legend = 
                        # for axis in axes:
                        #     axes[axis].hlines(y=threshold, xmin=0, xmax=1, linewidth=1, color='#650021', linestyles="dashed")
                        #     axes[axis].vlines(x=threshold, ymin=0, ymax=1, linewidth=1, color='#650021', linestyles="dashed")
                        #     axes[axis].set_xlabel("Donor score", labelpad=1.5)
                        #     axes[axis].set_ylabel("Acceptor score", labelpad=1.5)

                        #     if axis == "main":
                        #         title = ""#"All"
                        #     elif axis == "TPFN":
                        #         title = ""#"TP and FN"
                        #     elif axis == "FPTN":
                        #         title = ""#"FP and TN"
                        #     elif axis == "FPFN":
                        #         title = ""#"FP and FN"
                        #     axes[axis].set_title(title)
                            


                        # spliceai_TP__splam_TP = axes["main"].scatter(d_score[spliceai_TP__splam_TP_idices], a_score[spliceai_TP__splam_TP_idices], s = 0, color="pink", alpha=1.0)
                        # spliceai_TP__splam_TP_len = len(d_score[spliceai_TP__splam_TP_idices])

                        # spliceai_TP__splam_FN = axes["main"].scatter(d_score[spliceai_TP__splam_FN_idices], a_score[spliceai_TP__splam_FN_idices], s = 0.3, color="#13E8EC", alpha=1.0)
                        # spliceai_TP__splam_FN_len = len(d_score[spliceai_TP__splam_FN_idices])

                        # spliceai_FN__splam_TP = axes["main"].scatter(d_score[spliceai_FN__splam_TP_idices], a_score[spliceai_FN__splam_TP_idices], s = 0.3, color="#EC1713", alpha=1.0)
                        # spliceai_FN__splam_TP_len = len(d_score[spliceai_FN__splam_TP_idices])

                        # spliceai_FN__splam_FN = axes["main"].scatter(d_score[spliceai_FN__splam_FN_idices], a_score[spliceai_FN__splam_FN_idices], s = 0.3, color="#118F14", alpha=1.0)
                        # spliceai_FN__splam_FN_len = len(d_score[spliceai_FN__splam_FN_idices])


                        # spliceai_TN__splam_TN = axes["main"].scatter(d_score[spliceai_TN__splam_TN_idices], a_score[spliceai_TN__splam_TN_idices], s = 0, color="gray", alpha=1.0)
                        # spliceai_TN__splam_TN_len = len(d_score[spliceai_TN__splam_TN_idices])

                        # spliceai_TN__splam_FP = axes["main"].scatter(d_score[spliceai_TN__splam_FP_idices], a_score[spliceai_TN__splam_FP_idices], s = 0.3, color="#FFA300", alpha=1.0)
                        # spliceai_TN__splam_FP_len = len(d_score[spliceai_TN__splam_FP_idices])

                        # spliceai_FP__splam_TN = axes["main"].scatter(d_score[spliceai_FP__splam_TN_idices], a_score[spliceai_FP__splam_TN_idices], s = 0.3, color="#005CFF", alpha=1.0)
                        # spliceai_FP__splam_TN_len = len(d_score[spliceai_FP__splam_TN_idices])

                        # spliceai_FP__splam_FP = axes["main"].scatter(d_score[spliceai_FP__splam_FP_idices], a_score[spliceai_FP__splam_FP_idices], s = 0.3, color="#8F118C", alpha=1.0)
                        # spliceai_FP__splam_FP_len = len(d_score[spliceai_FP__splam_FP_idices])


                        # spliceai_TP__splam_FN = axes["TPFN"].scatter(d_score[spliceai_TP__splam_FN_idices], a_score[spliceai_TP__splam_FN_idices], s = 0.3, color="#13E8EC", alpha=1.0)
                        # spliceai_TP__splam_FN_len = len(d_score[spliceai_TP__splam_FN_idices])
                        # spliceai_FN__splam_TP = axes["TPFN"].scatter(d_score[spliceai_FN__splam_TP_idices], a_score[spliceai_FN__splam_TP_idices], s = 0.3, color="#EC1713", alpha=1.0)
                        # spliceai_FN__splam_TP_len = len(d_score[spliceai_FN__splam_TP_idices])

                        # spliceai_TN__splam_FP = axes["FPTN"].scatter(d_score[spliceai_TN__splam_FP_idices], a_score[spliceai_TN__splam_FP_idices], s = 0.3, color="#FFA300", alpha=1.0)
                        # spliceai_TN__splam_FP_len = len(d_score[spliceai_TN__splam_FP_idices])
                        # spliceai_FP__splam_TN = axes["FPTN"].scatter(d_score[spliceai_FP__splam_TN_idices], a_score[spliceai_FP__splam_TN_idices], s = 0.3, color="#005CFF", alpha=1.0)
                        # spliceai_FP__splam_TN_len = len(d_score[spliceai_FP__splam_TN_idices])

                        # spliceai_FP__splam_FP = axes["FPFN"].scatter(d_score[spliceai_FP__splam_FP_idices], a_score[spliceai_FP__splam_FP_idices], s = 0.3, color="#8F118C", alpha=1.0)
                        # spliceai_FP__splam_FP_len = len(d_score[spliceai_FP__splam_FP_idices])
                        # spliceai_FN__splam_FN = axes["FPFN"].scatter(d_score[spliceai_FN__splam_FN_idices], a_score[spliceai_FN__splam_FN_idices], s = 0.3, color="#118F14", alpha=1.0)
                        # spliceai_FN__splam_FN_len = len(d_score[spliceai_FN__splam_FN_idices])


                        # # non_junc_legend = ax.scatter(score_d[non_junc_prob], score_a[non_junc_prob], s = 0.3)
                        # # lgd = fig.legend([spliceai_TP__splam_TP, spliceai_TP__splam_FN, spliceai_FN__splam_TP, spliceai_FN__splam_FN, spliceai_TN__splam_TN, spliceai_TN__splam_FP, spliceai_FP__splam_TN, spliceai_FP__splam_FP], ["spliceai_TP__splam_TP ("+str(spliceai_TP__splam_TP_len)+")", "spliceai_TP__splam_FN ("+str(spliceai_TP__splam_FN_len)+")", "spliceai_FN__splam_TP ("+str(spliceai_FN__splam_TP_len)+")", "spliceai_FN__splam_FN ("+str(spliceai_FN__splam_FN_len)+")", "spliceai_TN__splam_TN ("+str(spliceai_TN__splam_TN_len)+")", "spliceai_TN__splam_FP ("+str(spliceai_TN__splam_FP_len)+")", "spliceai_FP__splam_TN ("+str(spliceai_FP__splam_TN_len)+")", "spliceai_FP__splam_FP ("+str(spliceai_FP__splam_FP_len)+")"], loc='lower right', bbox_to_anchor=(0.9, 0.27))

                        # print("Threshold:  ", threshold, "; str(spliceai_TP__splam_TP_len): ", str(spliceai_TP__splam_TP_len))
                        # print("Threshold:  ", threshold, "; str(spliceai_TN__splam_TN_len): ", str(spliceai_TN__splam_TN_len))


                        # lgd = fig.legend([spliceai_TP__splam_FN, spliceai_FN__splam_TP, spliceai_TN__splam_FP, spliceai_FP__splam_TN, spliceai_FN__splam_FN, spliceai_FP__splam_FP], ["spliceai_TP__splam_FN ("+str(spliceai_TP__splam_FN_len)+")", "spliceai_FN & splam_TP ("+str(spliceai_FN__splam_TP_len)+")", "spliceai_TN & splam_FP ("+str(spliceai_TN__splam_FP_len)+")", "spliceai_FP & splam_TN ("+str(spliceai_FP__splam_TN_len)+")", "spliceai_FN & splam_FN ("+str(spliceai_FN__splam_FN_len)+")", "spliceai_FP & splam_FP ("+str(spliceai_FP__splam_FP_len)+")"], loc='lower right', bbox_to_anchor=(0.9, 0.27), ncol=3)


                        # # lgd = ax.legend([spliceai_TP__splam_TP, spliceai_TP__splam_FN, spliceai_FN__splam_TP], ["spliceai_TP__splam_TP ("+str(spliceai_TP__splam_TP_len)+")", "spliceai_TP__splam_FN ("+str(spliceai_TP__splam_FN_len)+")", "spliceai_FN__splam_TP ("+str(spliceai_FN__splam_TP_len)+")"], loc='center left', bbox_to_anchor=(1, 0.5))

                        # # lgd = ax.legend([spliceai_FN__splam_FN, spliceai_FP__splam_FP], ["spliceai_FN__splam_FN ("+str(spliceai_FN__splam_FN_len)+")", "spliceai_FP__splam_FP ("+str(spliceai_FP__splam_FP_len)+")"], loc='center left', bbox_to_anchor=(1, 0.5))

                        # # lgd = axes.legend([spliceai_TN__splam_FP, spliceai_FP__splam_TN], ["spliceai_TN__splam_FP ("+str(spliceai_TN__splam_FP_len)+")", "spliceai_FP__splam_TN ("+str(spliceai_FP__splam_TN_len)+")"], loc='center left', bbox_to_anchor=(1, 0.5))

                        # axes["main"].set_aspect('equal', adjustable='box')
                        # axes["TPFN"].set_aspect('equal', adjustable='box')
                        # axes["FPTN"].set_aspect('equal', adjustable='box')
                        # axes["FPFN"].set_aspect('equal', adjustable='box')

                        # # handles, labels = .get_legend_handles_labels()
                        # # fig.legend(handles, labels, loc='lower center')
                        # plt.savefig(target_figure_root + TOOL + "_" + str(threshold) + "_" + MANE_OR_ALTS + ".png", dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')
                        # # fig.legend(labels, loc='lower right', bbox_to_anchor=(1,-0.1), ncol=len(labels), bbox_transform=fig.transFigure)
                        # # plt.savefig(target_figure_root + TOOL + "_" + str(threshold) + ".png", dpi=300, bbox_inches='tight')
                        # plt.close()


if __name__ == "__main__":
    main()


