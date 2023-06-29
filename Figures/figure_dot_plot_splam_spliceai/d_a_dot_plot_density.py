import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
from util import *
from scipy.stats import gaussian_kde
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, PrecisionRecallDisplay



def plot_scatter_plot(d_score, a_score, indices_1, indices_2, axis, cmap_1, cmap_2):
    selected_indices_combined = np.vstack([d_score, a_score])

    density_combined = gaussian_kde(selected_indices_combined)
    density_all_values = density_combined(np.vstack([d_score, a_score]))

    density_1_values = density_combined(np.vstack([d_score[indices_1], a_score[indices_1]]))
    density_2_values = density_combined(np.vstack([d_score[indices_2], a_score[indices_2]]))

    vmin = -10 #min(density_1_values.min(), density_2_values.min())-100
    vmax = density_all_values.max()
    norm = plt.Normalize(vmin=vmin, vmax=vmax)

    plt_res_1 = axis.scatter(d_score[indices_1], a_score[indices_1], cmap=cmap_1, s=0.3, c=density_1_values, norm=norm)
    plt_res_2 = axis.scatter(d_score[indices_2], a_score[indices_2], cmap=cmap_2, s=0.3, c=density_2_values, norm=norm)



    # new_2d_array = np.vstack([np.concatenate([d_score[indices_1], d_score[indices_2]]), np.concatenate([a_score[indices_1], a_score[indices_2]])])
    # print("new_2d_array: ", new_2d_array)
    # density = gaussian_kde(new_2d_array)
    # # Calculate density values for each data point
    # print("density: ", density)

    # density_values = density(new_2d_array)
    # print("density_values: ", density_values)

    # # Define colormap and normalize density values
    # # cmap = 'Reds'
    # norm = plt.Normalize(vmin=-5, vmax=max(density_values))
    # print("norm: ", norm)

    # plt_res_1 = axis.scatter(d_score[indices_1], a_score[indices_1], cmap=cmap_1, s = 0.3, c=density_values, norm = norm)
    # plt_res_2 = axis.scatter(d_score[indices_2], a_score[indices_2], cmap=cmap_2, s = 0.3, c=density_values, norm = norm)


    return plt_res_1, plt_res_2

    # junc_prob = label_d.astype(bool)
    # non_junc_prob = (1-label_d).astype(bool)
    # # print("\tspliceai_junc_prob: ", junc_prob)
    # # print("\tspliceai_junc_prob: ", non_junc_prob)

    # fig, ax = plt.subplots()
    # ax.set_title(filename)
    # ax.set_xlabel("Donor site score")
    # ax.set_ylabel("Acceptor site score")


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



# # Helper function used for visualization in the following examples
# def identify_axes(ax_dict, fontsize=48):
#     """
#     Helper to identify the Axes in the examples below.

#     Draws the label in a large font in the center of the Axes.

#     Parameters
#     ----------
#     ax_dict : dict[str, Axes]
#         Mapping between the title / label and the Axes.
#     fontsize : int, optional
#         How big the label should be.
#     """
#     kw = dict(ha="center", va="center", fontsize=fontsize, color="darkgrey")
#     for k, ax in ax_dict.items():
#         ax.text(0.5, 0.5, k, transform=ax.transAxes, **kw)

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
            #####################################
            # Creating directories for visualization.
            #####################################

            for SPLICEAI_VERSION in ["1", "2", "3", "4", "5", "AVERAGE"]:

                for TARGET in ["noN", "N"]:
                    figure_root = "./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/scatter_plot_density/"
                    target_figure_root = figure_root+TARGET+"/"
                    os.makedirs(target_figure_root+"tsv_"+str(threshold)+"/", exist_ok=True)
                    

                    with open("../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da."+TARGET+".merged.BOTH.pkl", "rb") as fr:
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


                    #     # spliceai_d_pred = [x.numpy() for x in spliceai_d_pred]
                    #     # spliceai_a_pred = [x.numpy() for x in spliceai_a_pred]
                    # spliceai_d_pred = np.array(spliceai_d_pred)
                    # spliceai_a_pred = np.array(spliceai_a_pred)
                    # spliceai_d_label = np.array(spliceai_d_label)
                    # spliceai_a_label = np.array(spliceai_a_label)
                    # spliceai_d_label = np.array(spliceai_d_label)

                    # print("spliceai_d_pred : ", spliceai_d_pred)
                    # print("spliceai_d_label: ", spliceai_d_label)
                    # print("spliceai_a_pred : ", spliceai_a_pred)
                    # print("spliceai_a_label: ", spliceai_a_label)
                    # print("spliceai_d_label: ", spliceai_d_label)

                    # print("spliceai_d_pred : ", len(spliceai_d_pred))
                    # print("spliceai_d_label: ", len(spliceai_d_label))
                    # print("spliceai_a_pred : ", len(spliceai_a_pred))
                    # print("spliceai_a_label: ", len(spliceai_a_label))

                    # print("spliceai_d_label: ", len(spliceai_d_label))


                    # with open("../../src_tools_evaluation/splam_result/splam.da.shuffle.pkl",'rb') as f:
                    #     splam_S_d_label_prob = pickle.load(f)
                    #     splam_S_d_pred_prob = pickle.load(f)
                        
                    #     splam_S_a_label_prob = pickle.load(f)
                    #     splam_S_a_pred_prob = pickle.load(f)

                    #     print("splam_S_d_label_prob : ", splam_S_d_label_prob)
                    #     print("splam_S_d_pred_prob: ", splam_S_d_pred_prob)

                    #     print("splam_S_a_label_prob : ", splam_S_a_label_prob)
                    #     print("splam_S_a_pred_prob: ", splam_S_a_pred_prob)

                    

                    with open("../../src_tools_evaluation/splam_result/"+SPLAM_VERSION+"/splam.da.noshuffle.merged.BOTH.pkl",'rb') as f:
                        splam_d_label = pickle.load(f)
                        splam_d_pred = pickle.load(f)
                        splam_a_label = pickle.load(f)
                        splam_a_pred = pickle.load(f)
                        print("\tsplam_d_pred : ", len(splam_d_pred))
                        print("\tsplam_d_pred: ", len(splam_d_pred))
                        print("")
                        print("\tsplam_a_label : ", len(splam_a_label))
                        print("\tsplam_a_pred: ", len(splam_a_pred))
                        print("")




                    # # print(spliceai_d_pred >= threshold)
                    # spliceai_N_TP_idices = (spliceai_N_d_pred_prob >= threshold) & (spliceai_N_a_pred_prob >= threshold) & (spliceai_N_d_label_prob == 1)
                    # spliceai_N_FN_idices = ((spliceai_N_d_pred_prob < threshold) | (spliceai_N_a_pred_prob < threshold)) & (spliceai_N_d_label_prob == 1)
                    # spliceai_N_FP_idices = (spliceai_N_d_pred_prob >= threshold) & (spliceai_N_a_pred_prob >= threshold) & (spliceai_N_d_label_prob == 0)
                    # spliceai_N_TN_idices = ((spliceai_N_d_pred_prob < threshold) | (spliceai_N_a_pred_prob < threshold)) & (spliceai_N_d_label_prob == 0)

                    # print("spliceai_TP_idices : ", len(spliceai_N_junc_name[spliceai_N_TP_idices]))
                    # print("spliceai_FP_idices : ", len(spliceai_N_junc_name[spliceai_N_FP_idices]))
                    # print("spliceai_FN_idices : ", len(spliceai_N_junc_name[spliceai_N_FN_idices]))
                    # print("spliceai_TN_idices : ", len(spliceai_N_junc_name[spliceai_N_TN_idices]))


                    # print(spliceai_d_pred >= threshold)
                    spliceai_TP_idices = (spliceai_d_pred >= threshold) & (spliceai_a_pred >= threshold) & (spliceai_d_label == 1)
                    spliceai_FN_idices = ((spliceai_d_pred < threshold) | (spliceai_a_pred < threshold)) & (spliceai_d_label == 1)
                    spliceai_FP_idices = (spliceai_d_pred >= threshold) & (spliceai_a_pred >= threshold) & (spliceai_d_label == 0)
                    spliceai_TN_idices = ((spliceai_d_pred < threshold) | (spliceai_a_pred < threshold)) & (spliceai_d_label == 0)

                    print("spliceai_TP_idices : ", len(spliceai_d_label[spliceai_TP_idices]))
                    print("spliceai_FP_idices : ", len(spliceai_d_label[spliceai_FP_idices]))
                    print("spliceai_FN_idices : ", len(spliceai_d_label[spliceai_FN_idices]))
                    print("spliceai_TN_idices : ", len(spliceai_d_label[spliceai_TN_idices]))


                    splam_TP_idices = (splam_d_pred >= threshold) & (splam_a_pred >= threshold) & (splam_d_label== 1)
                    splam_FN_idices = ((splam_d_pred < threshold) | (splam_a_pred < threshold)) & (splam_d_label == 1)
                    splam_FP_idices = (splam_d_pred >= threshold) & (splam_a_pred >= threshold) & (splam_d_label == 0)
                    splam_TN_idices = ((splam_d_pred < threshold) | (splam_a_pred < threshold)) & (splam_d_label == 0)

                    print("splam_TP_idices : ", len(splam_d_pred[splam_TP_idices]))
                    print("splam_FP_idices : ", len(splam_d_pred[splam_FP_idices]))
                    print("splam_FN_idices : ", len(splam_d_pred[splam_FN_idices]))
                    print("splam_TN_idices : ", len(splam_d_pred[splam_TN_idices]))

                    spliceai_TP__splam_TP_idices = spliceai_TP_idices & splam_TP_idices
                    spliceai_TP__splam_FN_idices = spliceai_TP_idices & splam_FN_idices
                    spliceai_FN__splam_TP_idices = spliceai_FN_idices & splam_TP_idices
                    spliceai_FN__splam_FN_idices = spliceai_FN_idices & splam_FN_idices

                    spliceai_TN__splam_TN_idices = spliceai_TN_idices & splam_TN_idices
                    spliceai_TN__splam_FP_idices = spliceai_TN_idices & splam_FP_idices
                    spliceai_FP__splam_TN_idices = spliceai_FP_idices & splam_TN_idices
                    spliceai_FP__splam_FP_idices = spliceai_FP_idices & splam_FP_idices

                    # print("spliceai_TP__splam_FN_idices : ", len(splam_d_pred[spliceai_TP__splam_FN_idices]))
                    # print("spliceai_FN__splam_TP_idices : ", len(splam_d_pred[spliceai_FN__splam_TP_idices]))
                    # print("spliceai_FN__splam_FN_idices : ", len(splam_d_pred[spliceai_FN__splam_FN_idices]))

                    # print("spliceai_TP__splam_FN_idices : ", len(splam_d_pred[spliceai_TN__splam_FP_idices]))
                    # print("spliceai_FN__splam_TP_idices : ", len(splam_d_pred[spliceai_FP__splam_TN_idices]))
                    # print("spliceai_FN__splam_FN_idices : ", len(splam_d_pred[spliceai_FP__splam_FP_idices]))

                    # junc_prob = label_d.astype(bool)
                    # non_junc_prob = (1-label_d).astype(bool)
                    # # print("\tspliceai_junc_prob: ", junc_prob)
                    # # print("\tspliceai_junc_prob: ", non_junc_prob)


                    for TYPE in ["TP_TP", "TP_FN", "FN_TP", "FN_FN", "TN_TN", "TN_FP", "FP_TN", "FP_FP"]:
                        if TYPE == "TP_TP":
                            indices = spliceai_TP__splam_TP_idices
                        elif TYPE == "TP_FN":
                            indices = spliceai_TP__splam_FN_idices
                        elif TYPE == "FN_TP":
                            indices = spliceai_FN__splam_TP_idices
                        elif TYPE == "FN_FN":
                            indices = spliceai_FN__splam_FN_idices
                        elif TYPE == "TN_TN":
                            indices = spliceai_TN__splam_TN_idices
                        elif TYPE == "TN_FP":
                            indices = spliceai_TN__splam_FP_idices
                        elif TYPE == "FP_TN":
                            indices = spliceai_FP__splam_TN_idices
                        elif TYPE == "FP_FP":
                            indices = spliceai_FP__splam_FP_idices

                        fw = open(target_figure_root+"tsv_"+str(threshold)+"/"+TYPE+"_junc_"+str(threshold)+".tsv", "w")
                        # for junc in splam_d_pred[indices]:
                        #     fw.write('\t'.join(junc)+"\n")
                        # fw.close()




                    plt.rcParams['font.size'] = 4
                    for TOOL in ["spliceai", "splam"]:
                        if TOOL == "spliceai":
                            d_score = spliceai_d_pred
                            a_score = spliceai_a_pred
                        elif TOOL == "splam":
                            d_score = splam_d_pred
                            a_score = splam_a_pred

                        # fig, axes = plt.subplot_mosaic(
                        #     [
                        #         ["main", "TPFN", "FPFN"],
                        #         ["main", "FPTN", "BLANK"],
                        #     ],
                        #     empty_sentinel="BLANK",
                        #     width_ratios=[3, 1.3, 1.3],
                        # )

                        fig, axes = plt.subplot_mosaic(
                            [
                                ["main", "TPFN", "FPTN", "FPFN"],
                            ],
                            empty_sentinel="BLANK",
                            width_ratios=[3, 1.9, 1.9, 1.9],
                        )
                        # identify_axes(axes, fontsize=36)
                        
                        plt.subplots_adjust(wspace=0.4, hspace=-0.5)

                        # junc_legend = 
                        for axis in axes:
                            axes[axis].hlines(y=threshold, xmin=0, xmax=1, linewidth=1, color='#650021', linestyles="dashed")
                            axes[axis].vlines(x=threshold, ymin=0, ymax=1, linewidth=1, color='#650021', linestyles="dashed")
                            axes[axis].set_xlabel("Donor score", labelpad=1.5)
                            axes[axis].set_ylabel("Acceptor score", labelpad=1.5)

                            if axis == "main":
                                title = ""#"All"
                            elif axis == "TPFN":
                                title = ""#"TP and FN"
                            elif axis == "FPTN":
                                title = ""#"FP and TN"
                            elif axis == "FPFN":
                                title = ""#"FP and FN"
                            axes[axis].set_title(title)

                        # spliceai_TP__splam_TP = axes["main"].scatter(d_score[spliceai_TP__splam_TP_idices], a_score[spliceai_TP__splam_TP_idices], cmap=cmap, s = 0, c=density_values, alpha=1.0, norm = norm)
                        # spliceai_TP__splam_TP = plot_scatter_plot(d_score, a_score, spliceai_TP__splam_TP_idices, axes["main"], 'Reds')
                        # spliceai_TP__splam_TP_len = len(d_score[spliceai_TP__splam_TP_idices])

                        # spliceai_TP__splam_FN_main = axes["main"].scatter(d_score[spliceai_TP__splam_FN_idices], a_score[spliceai_TP__splam_FN_idices], s = 0.3, color="#13E8EC", alpha=1.0)
                        # spliceai_FN__splam_TP_main = axes["main"].scatter(d_score[spliceai_FN__splam_TP_idices], a_score[spliceai_FN__splam_TP_idices], s = 0.3, color="#EC1713", alpha=1.0)
                        spliceai_TP__splam_FN_main, spliceai_FN__splam_TP_main= plot_scatter_plot(d_score, a_score, spliceai_TP__splam_FN_idices, spliceai_FN__splam_TP_idices, axes["main"], 'Blues', 'Reds')
                        spliceai_TP__splam_FN_len = len(d_score[spliceai_TP__splam_FN_idices])
                        spliceai_FN__splam_TP_len = len(d_score[spliceai_FN__splam_TP_idices])


                        # spliceai_TN__splam_TN = axes["main"].scatter(d_score[spliceai_TN__splam_TN_idices], a_score[spliceai_TN__splam_TN_idices], s = 0, color="gray", alpha=1.0)
                        # spliceai_TN__splam_TN_len = len(d_score[spliceai_TN__splam_TN_idices])

                        # spliceai_TN__splam_FP_main = axes["main"].scatter(d_score[spliceai_TN__splam_FP_idices], a_score[spliceai_TN__splam_FP_idices], s = 0.3, color="#FFA300", alpha=1.0)
                        # spliceai_FP__splam_TN_main = axes["main"].scatter(d_score[spliceai_FP__splam_TN_idices], a_score[spliceai_FP__splam_TN_idices], s = 0.3, color="#005CFF", alpha=1.0)

                        spliceai_TN__splam_FP_main, spliceai_FP__splam_TN_main= plot_scatter_plot(d_score, a_score, spliceai_TN__splam_FP_idices, spliceai_FP__splam_TN_idices, axes["main"], 'Purples', 'Greens')
                        spliceai_TN__splam_FP_len = len(d_score[spliceai_TN__splam_FP_idices])
                        spliceai_FP__splam_TN_len = len(d_score[spliceai_FP__splam_TN_idices])


                        # spliceai_FN__splam_FN_main = axes["main"].scatter(d_score[spliceai_FN__splam_FN_idices], a_score[spliceai_FN__splam_FN_idices], s = 0.3, color="#118F14", alpha=1.0)
                        # spliceai_FP__splam_FP_main = axes["main"].scatter(d_score[spliceai_FP__splam_FP_idices], a_score[spliceai_FP__splam_FP_idices], s = 0.3, color="#8F118C", alpha=1.0)
                        spliceai_FN__splam_FN_main, spliceai_FP__splam_FP_main= plot_scatter_plot(d_score, a_score, spliceai_FN__splam_FN_idices, spliceai_FP__splam_FP_idices, axes["main"], 'Greys', 'Oranges')
                        spliceai_FN__splam_FN_len = len(d_score[spliceai_FN__splam_FN_idices])
                        spliceai_FP__splam_FP_len = len(d_score[spliceai_FP__splam_FP_idices])




                        # spliceai_TP__splam_FN = axes["TPFN"].scatter(d_score[spliceai_TP__splam_FN_idices], a_score[spliceai_TP__splam_FN_idices], s = 0.3, color="#13E8EC", alpha=1.0)
                        spliceai_TP__splam_FN, spliceai_FN__splam_TP= plot_scatter_plot(d_score, a_score, spliceai_TP__splam_FN_idices, spliceai_FN__splam_TP_idices, axes["TPFN"], 'Blues', 'Reds')
                        spliceai_TP__splam_FN_len = len(d_score[spliceai_TP__splam_FN_idices])
                        # spliceai_FN__splam_TP = plot_scatter_plot(d_score, a_score, spliceai_FN__splam_TP_idices, axes["TPFN"], 'Reds')
                        # spliceai_FN__splam_TP = axes["TPFN"].scatter(d_score[spliceai_FN__splam_TP_idices], a_score[spliceai_FN__splam_TP_idices], s = 0.3, color="#EC1713", alpha=1.0)
                        spliceai_FN__splam_TP_len = len(d_score[spliceai_FN__splam_TP_idices])


                        spliceai_TN__splam_FP, spliceai_FP__splam_TN= plot_scatter_plot(d_score, a_score, spliceai_TN__splam_FP_idices, spliceai_FP__splam_TN_idices, axes["FPTN"], 'Blues', 'Reds')
                        # spliceai_TN__splam_FP = axes["FPTN"].scatter(d_score[spliceai_TN__splam_FP_idices], a_score[spliceai_TN__splam_FP_idices], s = 0.3, color="#FFA300", alpha=1.0)
                        spliceai_TN__splam_FP_len = len(d_score[spliceai_TN__splam_FP_idices])
                        # spliceai_FP__splam_TN = axes["FPTN"].scatter(d_score[spliceai_FP__splam_TN_idices], a_score[spliceai_FP__splam_TN_idices], s = 0.3, color="#005CFF", alpha=1.0)
                        spliceai_FP__splam_TN_len = len(d_score[spliceai_FP__splam_TN_idices])

                        # spliceai_FN__splam_FN, spliceai_FP__splam_FP= plot_scatter_plot(d_score, a_score, spliceai_FN__splam_FN_idices, spliceai_FP__splam_FP_idices, axes["FPFN"], 'Blues', 'Reds')

                        # spliceai_FP__splam_FP = axes["FPFN"].scatter(d_score[spliceai_FP__splam_FP_idices], a_score[spliceai_FP__splam_FP_idices], s = 0.3, color="#8F118C", alpha=1.0)
                        spliceai_FP__splam_FP_len = len(d_score[spliceai_FP__splam_FP_idices])
                        # spliceai_FN__splam_FN = axes["FPFN"].scatter(d_score[spliceai_FN__splam_FN_idices], a_score[spliceai_FN__splam_FN_idices], s = 0.3, color="#118F14", alpha=1.0)
                        spliceai_FN__splam_FN_len = len(d_score[spliceai_FN__splam_FN_idices])


                        # non_junc_legend = ax.scatter(score_d[non_junc_prob], score_a[non_junc_prob], s = 0.3)
                        # lgd = fig.legend([spliceai_TP__splam_TP, spliceai_TP__splam_FN, spliceai_FN__splam_TP, spliceai_FN__splam_FN, spliceai_TN__splam_TN, spliceai_TN__splam_FP, spliceai_FP__splam_TN, spliceai_FP__splam_FP], ["spliceai_TP__splam_TP ("+str(spliceai_TP__splam_TP_len)+")", "spliceai_TP__splam_FN ("+str(spliceai_TP__splam_FN_len)+")", "spliceai_FN__splam_TP ("+str(spliceai_FN__splam_TP_len)+")", "spliceai_FN__splam_FN ("+str(spliceai_FN__splam_FN_len)+")", "spliceai_TN__splam_TN ("+str(spliceai_TN__splam_TN_len)+")", "spliceai_TN__splam_FP ("+str(spliceai_TN__splam_FP_len)+")", "spliceai_FP__splam_TN ("+str(spliceai_FP__splam_TN_len)+")", "spliceai_FP__splam_FP ("+str(spliceai_FP__splam_FP_len)+")"], loc='lower right', bbox_to_anchor=(0.9, 0.27))

                        lgd = fig.legend([spliceai_TP__splam_FN_main, spliceai_FN__splam_TP_main, spliceai_TN__splam_FP_main, spliceai_FP__splam_TN_main, spliceai_FN__splam_FN_main, spliceai_FP__splam_FP_main], ["spliceai_TP__splam_FN ("+str(spliceai_TP__splam_FN_len)+")", "spliceai_FN & splam_TP ("+str(spliceai_FN__splam_TP_len)+")", "spliceai_TN & splam_FP ("+str(spliceai_TN__splam_FP_len)+")", "spliceai_FP & splam_TN ("+str(spliceai_FP__splam_TN_len)+")", "spliceai_FN & splam_FN ("+str(spliceai_FN__splam_FN_len)+")", "spliceai_FP & splam_FP ("+str(spliceai_FP__splam_FP_len)+")"], loc='lower right', bbox_to_anchor=(0.9, 0.29), ncol=3)


                        # lgd = ax.legend([spliceai_TP__splam_TP, spliceai_TP__splam_FN, spliceai_FN__splam_TP], ["spliceai_TP__splam_TP ("+str(spliceai_TP__splam_TP_len)+")", "spliceai_TP__splam_FN ("+str(spliceai_TP__splam_FN_len)+")", "spliceai_FN__splam_TP ("+str(spliceai_FN__splam_TP_len)+")"], loc='center left', bbox_to_anchor=(1, 0.5))

                        # lgd = ax.legend([spliceai_FN__splam_FN, spliceai_FP__splam_FP], ["spliceai_FN__splam_FN ("+str(spliceai_FN__splam_FN_len)+")", "spliceai_FP__splam_FP ("+str(spliceai_FP__splam_FP_len)+")"], loc='center left', bbox_to_anchor=(1, 0.5))

                        # lgd = axes.legend([spliceai_TN__splam_FP, spliceai_FP__splam_TN], ["spliceai_TN__splam_FP ("+str(spliceai_TN__splam_FP_len)+")", "spliceai_FP__splam_TN ("+str(spliceai_FP__splam_TN_len)+")"], loc='center left', bbox_to_anchor=(1, 0.5))

                        axes["main"].set_aspect('equal', adjustable='box')
                        axes["TPFN"].set_aspect('equal', adjustable='box')
                        axes["FPTN"].set_aspect('equal', adjustable='box')
                        axes["FPFN"].set_aspect('equal', adjustable='box')

                        # handles, labels = .get_legend_handles_labels()
                        # fig.legend(handles, labels, loc='lower center')
                        plt.savefig(target_figure_root + TOOL + "_" + str(threshold) + "_" + MANE_OR_ALTS + ".png", dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
                        # fig.legend(labels, loc='lower right', bbox_to_anchor=(1,-0.1), ncol=len(labels), bbox_transform=fig.transFigure)
                        # plt.savefig(target_figure_root + TOOL + "_" + str(threshold) + ".png", dpi=300, bbox_inches='tight')
                        plt.close()


if __name__ == "__main__":
    main()


